 /*@@
   @file      Spritz_Minima.cc
   @date      Tue Aug 29 18:52:10 2006
   @author
   @desc
        Sets up the scalars used for the atmosphere, after initial data.
   @enddesc
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include <cstdio>
#include <cmath>

//#include "Carpet/Carpet/src/carpet.hh"
#include "carpet.hh"

#ifdef HAVE_CARPET
using namespace Carpet;
#endif

#ifdef __cplusplus
  extern "C" {
#endif

/*#include "EinsteinBase/EOS_Base/src/EOS_Base.h"*/

    /* Scheduled functions */
    void Spritz_Rho_Minima_Setup_Final(CCTK_ARGUMENTS);

    void Spritz_Rho_Minima_Setup_Final_PUGH(CCTK_ARGUMENTS);

    void Spritz_ColdTabCalibration(CCTK_ARGUMENTS);

#ifdef __cplusplus
  } /* extern "C" */
#endif

/* Parse "RhoMin = <val>" from the header of an EOS_Omni cold_tabulated ASCII
 * table file. Used to calibrate Spritz's own matched-polytrope extrapolation
 * below the table's rho_min (see Spritz_coldtab_calib in interface.ccl),
 * without re-implementing EOS_Omni's table reader and without depending on
 * EOS_Omni_Module's internal (non-public-API) coldeos_rhomin. */
static int Spritz_ReadColdTabRhoMin(const char *table_name, CCTK_REAL *rho_min_out)
{
  FILE *fp = fopen(table_name, "r");
  if (!fp) return 0;

  char line[512];
  double rmin = 0.0, rmax = 0.0;
  int found = 0;
  while (fgets(line, sizeof(line), fp)) {
    if (sscanf(line, "RhoMin = %lf RhoMax = %lf", &rmin, &rmax) == 2) {
      found = 1;
      break;
    }
  }
  fclose(fp);

  if (!found || !(rmin > 0.0)) return 0;
  *rho_min_out = rmin;
  return 1;
}

/* Parse "GammaTh = <val>" from the same table file, the same way
 * EOS_Omni_ColdEOSReadTable.F90 does internally (coldeos_gammath is a
 * private EOS_Omni_Module Fortran variable with no exposed CCTK parameter,
 * so this is the only way to get the table's own thermal Gamma without
 * modifying EOS_Omni). Reading it straight from the file guarantees it
 * always matches whatever EOS_Omni itself is using, with no separate
 * Spritz-level value that could drift out of sync. */
static int Spritz_ReadColdTabGammaTh(const char *table_name, CCTK_REAL *gammath_out)
{
  FILE *fp = fopen(table_name, "r");
  if (!fp) return 0;

  char line[512];
  double gammath = 0.0;
  int found = 0;
  while (fgets(line, sizeof(line), fp)) {
    if (sscanf(line, "GammaTh = %lf", &gammath) == 1) {
      found = 1;
      break;
    }
  }
  fclose(fp);

  if (!found) return 0;
  *gammath_out = gammath;
  return 1;
}

 /*@@
   @routine    Spritz_ColdTabCalibration
   @desc
        Calibrate Spritz's own matched-polytrope extrapolation below the
        cold_tabulated EOS table's own rho_min. EOS_Omni's built-in
        extrapolation for this table (EOS_Omni_ColdEOSReadTable.F90)
        hardcodes Gamma=2 regardless of the table's actual local slope;
        doing our own calibration here lets any Spritz atmosphere-setting
        code (initial data in Spritz_Prim2ConInitial.c, the con2prim-failure
        reset in Spritz_Con2PrimM_Noble_pt_EOSOmni.c) use the table's real
        behavior instead, with zero EOS_Omni calls at the point of use --
        this one-time calibration, run here before either of those, is the
        only place that still calls EOS_Omni_press for this purpose.

        Scheduled AT CCTK_INITIAL AFTER HydroBase_Initial BEFORE
        HydroBase_Prim2ConInitial -- i.e. before Spritz_Prim2ConInitial runs,
        and well before Spritz_Rho_Minima_Setup_Final (AT
        CCTK_PostPostInitial), which used to compute this inline and now
        just reads the already-set Spritz_coldtab_* scalars.
   @enddesc
@@*/

void Spritz_ColdTabCalibration(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (!CCTK_EQUALS(Spritz_eos_table,"cold_tabulated")) return;

  CCTK_REAL rho_table_min;
  if (!Spritz_ReadColdTabRhoMin(coldeos_table_name, &rho_table_min))
    {
      CCTK_WARN(CCTK_WARN_ABORT, "Spritz: could not read RhoMin from EOS_Omni::coldeos_table_name for the coldtab matched-polytrope calibration.");
    }

  CCTK_REAL gammath;
  if (!Spritz_ReadColdTabGammaTh(coldeos_table_name, &gammath))
    {
      CCTK_WARN(CCTK_WARN_ABORT, "Spritz: could not read GammaTh from EOS_Omni::coldeos_table_name for the coldtab matched-polytrope calibration.");
    }
  const CCTK_REAL thfac = coldeos_use_thermal_gamma_law ? 1.0 : 0.0;

  const CCTK_INT n_calib = 1;
  CCTK_INT calib_keytemp = 1, calib_keyerr[1], calib_anyerr;
  CCTK_REAL cxrho[1], cxeps[1], cxpress[1], cxtemp[1], cxye[1];
  cxtemp[0] = Spritz_hot_atmo_temp;
  cxye[0]   = Spritz_Y_e_max;

  // First point: the table's own rho_min (gives P_lo, eps_lo).
  cxrho[0] = rho_table_min;
  EOS_Omni_press(*Spritz_eos_handle,calib_keytemp,Eos_rf_prec,n_calib,cxrho,cxeps,cxtemp,cxye,cxpress,calib_keyerr,&calib_anyerr);
  const CCTK_REAL P_lo   = cxpress[0];
  const CCTK_REAL eps_lo = cxeps[0];

  // Second, nearby point (WhiskyTHC-style two-point secant, see
  // thcextra/EOS_Barotropic/src/eos_tabulated.cpp) to recover the
  // table's actual local Gamma there, via the public EOS_Omni_press
  // API only (EOS_Omni does not expose the per-row Gamma it
  // interpolates internally).
  cxrho[0] = 1.1 * rho_table_min;
  EOS_Omni_press(*Spritz_eos_handle,calib_keytemp,Eos_rf_prec,n_calib,cxrho,cxeps,cxtemp,cxye,cxpress,calib_keyerr,&calib_anyerr);
  const CCTK_REAL P2 = cxpress[0];

  const CCTK_REAL Gamma_lo = log(P2/P_lo) / log(1.1);

  *Spritz_coldtab_rhomin   = rho_table_min;
  *Spritz_coldtab_press_lo = P_lo;
  *Spritz_coldtab_eps_lo   = eps_lo;
  *Spritz_coldtab_gamma_lo = Gamma_lo;
  *Spritz_coldtab_gammath  = gammath;
  *Spritz_coldtab_thfac    = thfac;

  CCTK_VInfo(CCTK_THORNSTRING,
    "coldtab matched-polytrope calibration: rho_min=%.6e P_lo=%.6e eps_lo=%.6e Gamma_lo=%.6f GammaTh=%.6f thfac=%.1f",
    (double)rho_table_min, (double)P_lo, (double)eps_lo, (double)Gamma_lo,
    (double)gammath, (double)thfac);
}

 /*@@
   @routine    Spritz_Rho_Minima_Setup_Final
   @date       Tue Aug 29 18:54:05 2006
   @author     Luca Baiotti
   @desc
        After initial data, set up the scalar Spritz_rho_min used for checking the atmosphere.
        This is computed taking into account the actual maximum rest-mass density on the grid.
	Version for Carpet.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory
@@*/

void Spritz_Rho_Minima_Setup_Final(CCTK_ARGUMENTS)
{

#ifdef HAVE_CARPET

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL  max_rho;

  static int flag = true;

  if (flag) // Actually run this reduction only the first time the routine is called.
    {
      if (rho_abs_min > 0.0)
	{
	  *Spritz_rho_min = rho_abs_min;
	}
      else
	{
	  // Go to global mode

	  BEGIN_GLOBAL_MODE (cctkGH) {

            // Find the global maximum of rho

            const int Reduction_Handle = CCTK_ReductionHandle ("maximum");

            const CCTK_INT input_array_variable_indices= {CCTK_VarIndex("HydroBase::rho")};

            const int ierr = CCTK_Reduce(cctkGH,
					 -1, // target processors; -1 -> all
					 Reduction_Handle,
					 1,  // number of output variables
					 CCTK_VARIABLE_REAL,
					 &max_rho,
					 1,  // number of variables to be reduced
					 input_array_variable_indices);

            if (ierr != 0)
	      {
		CCTK_WARN(0, "Failed to compute the global maximum of rho");
	      }

            *Spritz_rho_min = max_rho * rho_rel_min;

            // Go back to local mode
	  } END_GLOBAL_MODE;

	}

      // Spritz_coldtab_calib (rho_min/press_lo/eps_lo/gamma_lo) is now
      // computed earlier by Spritz_ColdTabCalibration (AT CCTK_INITIAL,
      // before HydroBase_Prim2ConInitial, so both Spritz_Prim2ConInitial and
      // this routine's own Spritz_Cs2_min setup below can use it) --
      // nothing to do here anymore.

      // After this has run once, set the flag so that this does not run again
      flag = false;
    } // end if (flag)
  else
    {
      return;
    }

#define  eps_dummy   0.0
  CCTK_REAL rho_atm, press_atm, eps_atm, dpdrho_atm, dpdeps_atm;

  /*!!!!!!
  ! begin EOS Omni vars*/
  /*CCTK_REAL  Eos_rf_prec;*/
  CCTK_INT n,keytemp,anyerr,keyerr[1];  /* * */
  CCTK_REAL xrho[1],xpress[1],xpressinit[1],xeps[1],xtemp[1],xye[1];
  CCTK_REAL xdpdeps[1],xdpdrho[1],xcs2[1];

  n = 1;
  keytemp = 0;
  keyerr[0] = 0;
  xpress[0] = 0.0;
  xpressinit[0] = eps_dummy;
  xeps[0] = eps_dummy;
  xtemp[0] = 0.0;
  xye[0] = 0.0;
  xcs2[0] = 0.0;
  /*! end EOS Omni vars
  !!!!!!*/

  rho_atm      = *Spritz_rho_min;
  /*press_atm    = EOS_Pressure         ( *Spritz_polytrope_handle, rho_atm, eps_dummy );
  eps_atm      = EOS_SpecificIntEnergy( *Spritz_polytrope_handle, rho_atm, eps_dummy );
  dpdrho_atm   = EOS_DPressByDRho     ( *Spritz_polytrope_handle, rho_atm, eps_dummy );
  dpdeps_atm   = EOS_DPressByDEps     ( *Spritz_polytrope_handle, rho_atm, eps_dummy );*/

  /*!!!!!!*/
  xrho[0] = rho_atm;
  xtemp[0]= Spritz_hot_atmo_temp;
  xye[0]  = Spritz_Y_e_max;
  keytemp = 1;

  if(CCTK_EQUALS(Spritz_eos_table,"nuc_eos")) {
      EOS_Omni_press(*Spritz_eos_handle,keytemp,Eos_rf_prec,n,xrho,xeps,xtemp,xye,xpress,keyerr,&anyerr);
      press_atm = xpress[0];
      eps_atm = xeps[0];
      EOS_Omni_cs2(*Spritz_eos_handle,keytemp,Eos_rf_prec,n,xrho,xeps,xtemp,xye,xcs2,keyerr,&anyerr);
      *Spritz_Cs2_min = xcs2[0];
  }
  else if (CCTK_EQUALS(Spritz_eos_table,"cold_tabulated")) {
      /* Matched-polytrope extrapolation below the table's own rho_min,
       * calibrated once above into Spritz_coldtab_calib -- zero EOS_Omni
       * calls here. Uses the pure power-law ("multiplicative") form
       * eps=eps_lo*ratio^(Gamma_lo-1), P=P_lo*ratio^Gamma_lo, matching
       * GRHayL's NRPyEOS_hybrid_coldtab_compute_P_cold_and_eps_cold.c
       * exactly -- verified numerically that the alternative additive form
       * (integrating deps=-P d(1/rho) from the matching point) goes
       * NEGATIVE at Spritz_rho_min for this table (rho_min/Spritz_rho_min
       * ~200), consistent with that file's own comment explaining why it
       * was rejected there. The multiplicative form is a pure power law
       * (always positive) and is exactly continuous in both P and eps at
       * the matching point (ratio=1) by construction. dpdeps_atm=0 below
       * rho_min matches EOS_Omni's own EOS_Omni_EOS_DPressByDEps case(5)
       * (P=P(rho) only in this regime, confirmed against source);
       * dpdrho_atm is the analytic derivative of the same power law used
       * for press_atm. */
      const CCTK_REAL Gamma_lo = *Spritz_coldtab_gamma_lo;
      const CCTK_REAL ratio = rho_atm / (*Spritz_coldtab_rhomin);
      press_atm  = (*Spritz_coldtab_press_lo) * pow(ratio, Gamma_lo);
      eps_atm    = (*Spritz_coldtab_eps_lo) * pow(ratio, Gamma_lo - 1.0);
      dpdeps_atm = 0.0;
      dpdrho_atm = Gamma_lo * press_atm / rho_atm;
      *Spritz_Cs2_min = ( dpdrho_atm + press_atm * dpdeps_atm / ( rho_atm*rho_atm ) )
        / ( 1.0 + eps_atm + press_atm / rho_atm );
  }
  else {
      EOS_Omni_press(*Spritz_eos_handle,keytemp,Eos_rf_prec,n,xrho,xeps,xtemp,xye,xpress,keyerr,&anyerr);
      press_atm = xpress[0];
      /* NOTE: EOS_Omni_EpsFromPress used to be called here to set eps_atm,
       * but with keytemp=1 the EOS_Omni_press call above already writes the
       * cold-curve eps directly into xeps -- the EpsFromPress call was pure
       * redundancy (same handle, same keytemp, no new inputs). It also has
       * no case(5) (cold_tabulated) and hard-aborts the run if this EOS is
       * in use, so it's dropped entirely; xeps[0] already holds the right
       * value from the call above. */
      eps_atm = xeps[0];
      EOS_Omni_DPressByDEps(*Spritz_eos_handle,keytemp,Eos_rf_prec,n,xrho,xeps,xtemp,xye,xdpdeps,keyerr,&anyerr);
      dpdeps_atm = xdpdeps[0];
      EOS_Omni_DPressByDRho(*Spritz_eos_handle,keytemp,Eos_rf_prec,n,xrho,xeps,xtemp,xye,xdpdrho,keyerr,&anyerr);
      dpdrho_atm = xdpdrho[0];
      /*!!!!!!*/
      *Spritz_Cs2_min = ( dpdrho_atm + press_atm * dpdeps_atm / ( rho_atm*rho_atm ) )
        / ( 1.0 + eps_atm + press_atm / rho_atm );
  }
  if ( *Spritz_Cs2_min < 0.0 ) CCTK_WARN( CCTK_WARN_ABORT, "Spritz_Cs2_min is negative.");

#undef  eps_dummy

#endif

  //Debug stuff
  // char warnline;

  //CCTK_VInfo (CCTK_THORNSTRING, "STEP 2: compute rho max; rho min: %13.12e \n",*Spritz_rho_min);

 //  printf(warnline,"STEP 2: compute rho max; rho min: %13.12e \n",*Spritz_rho_min);
 //CCTK_WARN(1,warnline);

  //  printf(warnline,"STEP 3: recompute ID with new atmosphere; rho min: ', Spritz_rho_min, 'reflev: ',%i), Spritz_rho_min, Spritz_reflevel);



}


 /*@@
   @routine    Spritz_Rho_Minima_Setup_Final_PUGH
   @date       Tue Aug 29 18:57:42 2006
   @author     Luca Baiotti
   @desc
         As above, but for PUGH.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory
@@*/

void Spritz_Rho_Minima_Setup_Final_PUGH(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL  max_rho;

  if (rho_abs_min > 0.0)
    {
      *Spritz_rho_min = rho_abs_min;
    }
  else
    {
      // Find the global maximum of rho

      const int Reduction_Handle = CCTK_ReductionHandle("maximum");

      const CCTK_INT input_array_variable_indices={CCTK_VarIndex("HydroBase::rho")};

      const int ierr = CCTK_Reduce(cctkGH,
			 -1, // target processors; -1 -> all
			 Reduction_Handle,
			 1,  // number of output variables
			 CCTK_VARIABLE_REAL,
			 &max_rho,
			 1,  // number of variables to be reduced
			 input_array_variable_indices);

      if (ierr != 0)
	{
	  CCTK_WARN(0, "Failed to compute the global maximum of rho");
	}

      *Spritz_rho_min = max_rho * rho_rel_min;
    }

#define  eps_dummy   0.0
  CCTK_REAL rho_atm, press_atm, eps_atm, dpdrho_atm, dpdeps_atm;

  /*!!!!!!
  ! begin EOS Omni vars*/
  /*CCTK_REAL  Eos_rf_prec;*/
  CCTK_INT n,keytemp,anyerr,keyerr[1]; /* * */
  CCTK_REAL xrho[1],xpress[1],xpressinit[1],xeps[1],xtemp[1],xye[1];
  CCTK_REAL xdpdeps[1],xdpdrho[1],xcs2[1];

  n = 1;
  keytemp = 0;
  keyerr[0] = 0;
  xpress[0] = 0.0;
  xpressinit[0] = eps_dummy;
  xeps[0] = eps_dummy;
  xtemp[0] = 0.0;
  xye[0] = 0.0;
  /*! end EOS Omni vars
  !!!!!!*/

  rho_atm      = *Spritz_rho_min;
  /*press_atm    = EOS_Pressure         ( *Spritz_polytrope_handle, rho_atm, eps_dummy );
  eps_atm      = EOS_SpecificIntEnergy( *Spritz_polytrope_handle, rho_atm, eps_dummy );
  dpdrho_atm   = EOS_DPressByDRho     ( *Spritz_polytrope_handle, rho_atm, eps_dummy );
  dpdeps_atm   = EOS_DPressByDEps     ( *Spritz_polytrope_handle, rho_atm, eps_dummy );*/

  /*!!!!!!*/
  xrho[0] = rho_atm;
  xtemp[0]= Spritz_hot_atmo_temp;
  xye[0]  = Spritz_Y_e_max;
  keytemp = 1;
  EOS_Omni_press(*Spritz_eos_handle,keytemp,Eos_rf_prec,n,xrho,xeps,xtemp,xye,xpress,keyerr,&anyerr);
  press_atm = xpress[0];
  /* NOTE: see the identical case above -- keytemp=1 already makes the
   * EOS_Omni_press call write the cold-curve eps into xeps directly, so the
   * EpsFromPress call that used to sit here was pure redundancy (and has no
   * case(5)/cold_tabulated support, hard-aborting the run). Dropped. */
  eps_atm = xeps[0];
  EOS_Omni_DPressByDEps(*Spritz_eos_handle,keytemp,Eos_rf_prec,n,xrho,xeps,xtemp,xye,xdpdeps,keyerr,&anyerr);
  dpdeps_atm = xdpdeps[0];
  EOS_Omni_DPressByDRho(*Spritz_eos_handle,keytemp,Eos_rf_prec,n,xrho,xeps,xtemp,xye,xdpdrho,keyerr,&anyerr);
  dpdrho_atm = xdpdrho[0];
  EOS_Omni_DPressByDRho(*Spritz_eos_handle,keytemp,Eos_rf_prec,n,xrho,xeps,xtemp,xye,xcs2,keyerr,&anyerr);
  *Spritz_Cs2_min = xcs2[0];

  if ( *Spritz_Cs2_min < 0.0 )
    { CCTK_WARN( CCTK_WARN_ABORT, "Spritz_Cs2_min is negative."); }
#undef  eps_dummy

}
