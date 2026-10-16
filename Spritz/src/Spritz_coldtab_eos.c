#include <math.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

/*
 * Centralized matched-polytrope replacement for EOS_Omni_press/
 * DPressByDRho/DPressByDEps, used for EVERY cold_tabulated EOS call in
 * Spritz's active production path (Noble con2prim, prim2con_mhd,
 * Spritz_HLLE_MHD reconstruction) -- not just the con2prim-failure
 * atmosphere reset and initial-data floor, which already had their own
 * hand-written copies of this formula (Spritz_Con2PrimM_Noble_pt_EOSOmni.c's
 * atmosphere-reset branch, Spritz_Prim2ConInitial.c).
 *
 * EOS_Omni's own cold_tabulated extrapolation below the table's rho_min
 * hardcodes Gamma=2 (EOS_Omni_ColdEOSReadTable.F90), even though it already
 * reads the table's real local Gamma into coldeos_gamma(1) two lines
 * earlier. Rather than patching EOS_Omni itself (a shared, native ET
 * thorn), every Spritz call site that can see rho below the table's
 * rho_min is routed through here instead, using the real local Gamma
 * (coldtab_gamma_lo) calibrated once at startup by Spritz_ColdTabCalibration
 * (Spritz_Minima.cc) via a two-point secant against the actual table.
 *
 * Cold-part formula: a pure power law anchored at the table's own first row
 * (press_lo, eps_lo) -- continuous in P and in dP/drho (hence c_s^2) at
 * the matching point, by construction. This is the GRHayL-style choice
 * (see NRPyEOS_hybrid_coldtab_compute_P_cold_and_eps_cold.c), deliberately
 * NOT the WhiskyTHC-style thermodynamic-identity closure, which trades
 * that c_s^2 continuity for exact first-law consistency instead -- the
 * two are not simultaneously achievable with a single Gamma (verified
 * analytically and numerically earlier in this investigation).
 *
 * Thermal part: mirrors EOS_Omni's own case(5) table-INTERIOR formula
 * exactly (EOS_Omni_SingleVarCalls.F90's EOS_Omni_EOS_Press/DPressByDRho/
 * DPressByDEps, case(5), rho>=rhomin branch) --
 *   eps_th    = max(0, eps_input - eps_cold)   if keytemp==0 (eps is a
 *               genuine input, possibly with real thermal content)
 *             = 0, and eps forced to eps_cold  if keytemp==1 (atmosphere
 *               reset / pure-cold request)
 *   press_th  = coldtab_thfac*(coldtab_gammath-1)*rho*eps_th
 *   dP/drho   = Gamma_lo*press_cold/rho + press_th/rho
 *               (Gamma_lo*press_cold/rho is exactly dP_cold/drho for the
 *               pure power law above; EOS_Omni's own interior formula is
 *               cs2_cold*h_cold + press_th/rho, and cs2_cold*h_cold is by
 *               definition dP_cold/drho, so this is the same identity)
 *   dP/deps   = coldtab_thfac*(coldtab_gammath-1)*rho
 * coldtab_gammath and coldtab_thfac are calibrated once at startup by
 * Spritz_ColdTabCalibration (Spritz_Minima.cc), which parses the table
 * file's own "GammaTh = ..." header line directly (coldeos_gammath is a
 * private EOS_Omni_Module Fortran variable with no exposed CCTK parameter)
 * and reads the public EOS_Omni::coldeos_use_thermal_gamma_law parameter
 * for thfac -- so this is never previously computed by EOS_Omni's own
 * extrapolated branch, which drops the thermal term entirely below
 * rho_min (EOS_Omni_ColdEOSReadTable.F90's "cycle" statements skip past
 * it) -- a second, independent defect from the hardcoded Gamma=2, stacked
 * in the same region. Fixing it here is strictly new correct behavior,
 * not a restoration of something EOS_Omni itself ever did.
 *
 * For rho >= rho_min (or for any EOS table other than cold_tabulated,
 * signaled by coldtab_rhomin<=0 since Spritz_ColdTabCalibration only ever
 * sets it when Spritz_eos_table=="cold_tabulated"), these fall through to
 * the real EOS_Omni_* calls, unchanged.
 */

void Spritz_coldtab_press(
    CCTK_INT handle, CCTK_INT keytemp, CCTK_REAL rf_precision, CCTK_INT npoints,
    CCTK_REAL *rho, CCTK_REAL *eps, CCTK_REAL *temp, CCTK_REAL *ye,
    CCTK_REAL *press, CCTK_INT *keyerr, CCTK_INT *anyerr,
    CCTK_REAL coldtab_rhomin, CCTK_REAL coldtab_press_lo,
    CCTK_REAL coldtab_eps_lo, CCTK_REAL coldtab_gamma_lo,
    CCTK_REAL coldtab_gammath, CCTK_REAL coldtab_thfac)
{
  for (CCTK_INT i = 0; i < npoints; i++) {
    if (coldtab_rhomin > 0.0 && rho[i] < coldtab_rhomin) {
      const CCTK_REAL ratio = rho[i] / coldtab_rhomin;
      const CCTK_REAL eps_cold   = coldtab_eps_lo   * pow(ratio, coldtab_gamma_lo - 1.0);
      const CCTK_REAL press_cold = coldtab_press_lo * pow(ratio, coldtab_gamma_lo);
      CCTK_REAL eps_th;
      if (keytemp == 1) {
        eps_th = 0.0;
        eps[i] = eps_cold;
      } else {
        eps_th = fmax(0.0, eps[i] - eps_cold);
      }
      const CCTK_REAL press_th = coldtab_thfac * (coldtab_gammath - 1.0) * rho[i] * eps_th;
      press[i] = press_cold + press_th;
      keyerr[i] = 0;
    } else {
      CCTK_REAL rho1[1]   = {rho[i]};
      CCTK_REAL eps1[1]   = {eps[i]};
      CCTK_REAL temp1[1]  = {temp[i]};
      CCTK_REAL ye1[1]    = {ye[i]};
      CCTK_REAL press1[1] = {0.0};
      CCTK_INT  keyerr1[1] = {0};
      EOS_Omni_press(handle, keytemp, rf_precision, 1, rho1, eps1, temp1, ye1,
                      press1, keyerr1, anyerr);
      press[i]  = press1[0];
      eps[i]    = eps1[0];
      temp[i]   = temp1[0];
      keyerr[i] = keyerr1[0];
    }
  }
}

void Spritz_coldtab_dpdrho(
    CCTK_INT handle, CCTK_INT keytemp, CCTK_REAL rf_precision, CCTK_INT npoints,
    CCTK_REAL *rho, CCTK_REAL *eps, CCTK_REAL *temp, CCTK_REAL *ye,
    CCTK_REAL *dpdrho, CCTK_INT *keyerr, CCTK_INT *anyerr,
    CCTK_REAL coldtab_rhomin, CCTK_REAL coldtab_press_lo,
    CCTK_REAL coldtab_eps_lo, CCTK_REAL coldtab_gamma_lo,
    CCTK_REAL coldtab_gammath, CCTK_REAL coldtab_thfac)
{
  for (CCTK_INT i = 0; i < npoints; i++) {
    if (coldtab_rhomin > 0.0 && rho[i] < coldtab_rhomin) {
      const CCTK_REAL ratio = rho[i] / coldtab_rhomin;
      const CCTK_REAL eps_cold   = coldtab_eps_lo   * pow(ratio, coldtab_gamma_lo - 1.0);
      const CCTK_REAL press_cold = coldtab_press_lo * pow(ratio, coldtab_gamma_lo);
      CCTK_REAL eps_th;
      if (keytemp == 1) {
        eps_th = 0.0;
        eps[i] = eps_cold;
      } else {
        eps_th = fmax(0.0, eps[i] - eps_cold);
      }
      const CCTK_REAL press_th = coldtab_thfac * (coldtab_gammath - 1.0) * rho[i] * eps_th;
      /* dP/drho|_eps = dP_cold/drho + press_th/rho, matching EOS_Omni's own
       * case(5) formula (dpdrhoe = cs2_cold*h_cold + press_th/rho, where
       * cs2_cold*h_cold is dP_cold/drho by definition of cs2_cold); for the
       * pure power law P_cold=K*rho^Gamma_lo, dP_cold/drho=Gamma_lo*P_cold/rho
       * exactly. */
      dpdrho[i] = (coldtab_gamma_lo * press_cold + press_th) / rho[i];
      keyerr[i] = 0;
    } else {
      CCTK_REAL rho1[1]    = {rho[i]};
      CCTK_REAL eps1[1]    = {eps[i]};
      CCTK_REAL temp1[1]   = {temp[i]};
      CCTK_REAL ye1[1]     = {ye[i]};
      CCTK_REAL dpdrho1[1] = {0.0};
      CCTK_INT  keyerr1[1] = {0};
      EOS_Omni_DPressByDRho(handle, keytemp, rf_precision, 1, rho1, eps1, temp1, ye1,
                             dpdrho1, keyerr1, anyerr);
      dpdrho[i] = dpdrho1[0];
      eps[i]    = eps1[0];
      temp[i]   = temp1[0];
      keyerr[i] = keyerr1[0];
    }
  }
}

void Spritz_coldtab_dpdeps(
    CCTK_INT handle, CCTK_INT keytemp, CCTK_REAL rf_precision, CCTK_INT npoints,
    CCTK_REAL *rho, CCTK_REAL *eps, CCTK_REAL *temp, CCTK_REAL *ye,
    CCTK_REAL *dpdeps, CCTK_INT *keyerr, CCTK_INT *anyerr,
    CCTK_REAL coldtab_rhomin, CCTK_REAL coldtab_press_lo,
    CCTK_REAL coldtab_eps_lo, CCTK_REAL coldtab_gamma_lo,
    CCTK_REAL coldtab_gammath, CCTK_REAL coldtab_thfac)
{
  for (CCTK_INT i = 0; i < npoints; i++) {
    if (coldtab_rhomin > 0.0 && rho[i] < coldtab_rhomin) {
      if (keytemp == 1) {
        const CCTK_REAL ratio = rho[i] / coldtab_rhomin;
        eps[i] = coldtab_eps_lo * pow(ratio, coldtab_gamma_lo - 1.0);
      }
      /* dP/deps|_rho = thfac*(Gamma_th-1)*rho, matching EOS_Omni's own
       * case(5) formula exactly -- independent of keytemp/eps_th, same as
       * EOS_Omni's own interior-region implementation. */
      dpdeps[i] = coldtab_thfac * (coldtab_gammath - 1.0) * rho[i];
      keyerr[i] = 0;
    } else {
      CCTK_REAL rho1[1]    = {rho[i]};
      CCTK_REAL eps1[1]    = {eps[i]};
      CCTK_REAL temp1[1]   = {temp[i]};
      CCTK_REAL ye1[1]     = {ye[i]};
      CCTK_REAL dpdeps1[1] = {0.0};
      CCTK_INT  keyerr1[1] = {0};
      EOS_Omni_DPressByDEps(handle, keytemp, rf_precision, 1, rho1, eps1, temp1, ye1,
                             dpdeps1, keyerr1, anyerr);
      dpdeps[i] = dpdeps1[0];
      eps[i]    = eps1[0];
      temp[i]   = temp1[0];
      keyerr[i] = keyerr1[0];
    }
  }
}

/* Fortran-callable entry points (CCTK_FCALL/CCTK_FNAME convention, same as
 * Spritz_Con2Prim_Noble): Fortran passes every argument by reference, so
 * these just dereference the scalars and forward to the plain C functions
 * above, which both this file's own callers and Spritz_Con2PrimM_Noble_pt_
 * EOSOmni.c (a plain C call, no Fortran involved there) use directly. */

void CCTK_FCALL CCTK_FNAME(Spritz_coldtab_press)(
    CCTK_INT *handle, CCTK_INT *keytemp, CCTK_REAL *rf_precision, CCTK_INT *npoints,
    CCTK_REAL *rho, CCTK_REAL *eps, CCTK_REAL *temp, CCTK_REAL *ye,
    CCTK_REAL *press, CCTK_INT *keyerr, CCTK_INT *anyerr,
    CCTK_REAL *coldtab_rhomin, CCTK_REAL *coldtab_press_lo,
    CCTK_REAL *coldtab_eps_lo, CCTK_REAL *coldtab_gamma_lo,
    CCTK_REAL *coldtab_gammath, CCTK_REAL *coldtab_thfac)
{
  Spritz_coldtab_press(*handle, *keytemp, *rf_precision, *npoints,
                        rho, eps, temp, ye, press, keyerr, anyerr,
                        *coldtab_rhomin, *coldtab_press_lo,
                        *coldtab_eps_lo, *coldtab_gamma_lo,
                        *coldtab_gammath, *coldtab_thfac);
}

void CCTK_FCALL CCTK_FNAME(Spritz_coldtab_dpdrho)(
    CCTK_INT *handle, CCTK_INT *keytemp, CCTK_REAL *rf_precision, CCTK_INT *npoints,
    CCTK_REAL *rho, CCTK_REAL *eps, CCTK_REAL *temp, CCTK_REAL *ye,
    CCTK_REAL *dpdrho, CCTK_INT *keyerr, CCTK_INT *anyerr,
    CCTK_REAL *coldtab_rhomin, CCTK_REAL *coldtab_press_lo,
    CCTK_REAL *coldtab_eps_lo, CCTK_REAL *coldtab_gamma_lo,
    CCTK_REAL *coldtab_gammath, CCTK_REAL *coldtab_thfac)
{
  Spritz_coldtab_dpdrho(*handle, *keytemp, *rf_precision, *npoints,
                         rho, eps, temp, ye, dpdrho, keyerr, anyerr,
                         *coldtab_rhomin, *coldtab_press_lo,
                         *coldtab_eps_lo, *coldtab_gamma_lo,
                         *coldtab_gammath, *coldtab_thfac);
}

void CCTK_FCALL CCTK_FNAME(Spritz_coldtab_dpdeps)(
    CCTK_INT *handle, CCTK_INT *keytemp, CCTK_REAL *rf_precision, CCTK_INT *npoints,
    CCTK_REAL *rho, CCTK_REAL *eps, CCTK_REAL *temp, CCTK_REAL *ye,
    CCTK_REAL *dpdeps, CCTK_INT *keyerr, CCTK_INT *anyerr,
    CCTK_REAL *coldtab_rhomin, CCTK_REAL *coldtab_press_lo,
    CCTK_REAL *coldtab_eps_lo, CCTK_REAL *coldtab_gamma_lo,
    CCTK_REAL *coldtab_gammath, CCTK_REAL *coldtab_thfac)
{
  Spritz_coldtab_dpdeps(*handle, *keytemp, *rf_precision, *npoints,
                         rho, eps, temp, ye, dpdeps, keyerr, anyerr,
                         *coldtab_rhomin, *coldtab_press_lo,
                         *coldtab_eps_lo, *coldtab_gamma_lo,
                         *coldtab_gammath, *coldtab_thfac);
}
