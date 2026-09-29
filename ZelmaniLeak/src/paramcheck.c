#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <assert.h>

double min(double a, double b) {
  return a < b ? a : b;
}

void ZelmaniLeak_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *in_merger = 0; //at the very beginning, no merger. LSMOD

  if(force_postmerger_at_time > 0.0 && force_postmerger == 0) { //LSMOD
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Parameter force_postmerger_at_time won't be taken into account, since force_postmerger is set to 'no'. Please, set force_postmerger to yes if you want to force leakage.");
  }

  // rad_max_outer = min(rad_max_outer,min(xmax,min(ymax,zmax)))

  if(do_pnu) {
    // check that Spritz RHS are available when adding pressure terms
    int varindex, rhsindex, numtl;

    varindex = CCTK_VarIndex("Spritz::tau");
    rhsindex = MoLQueryEvolvedRHS (varindex);
    numtl = CCTK_ActiveTimeLevelsVI(cctkGH, rhsindex);
    if(numtl <= 0) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "No storage for %s. varindex = %d, rhsindex = %d, numtl = %d",
                  "taurhs", varindex, rhsindex, numtl);
    }

    varindex = CCTK_VarIndex("Spritz::scon[0]");
    rhsindex = MoLQueryEvolvedRHS (varindex);
    numtl = CCTK_ActiveTimeLevelsVI(cctkGH, rhsindex);
    if(numtl <= 0) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "No storage for %s. varindex = %d, rhsindex = %d, numtl = %d",
                  "sxrhs", varindex, rhsindex, numtl);
    }

    varindex = CCTK_VarIndex("Spritz::scon[1]");
    rhsindex = MoLQueryEvolvedRHS (varindex);
    numtl = CCTK_ActiveTimeLevelsVI(cctkGH, rhsindex);
    if(numtl <= 0) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "No storage for %s. varindex = %d, rhsindex = %d, numtl = %d",
                  "syrhs", varindex, rhsindex, numtl);
    }

    varindex = CCTK_VarIndex("Spritz::scon[2]");
    rhsindex = MoLQueryEvolvedRHS (varindex);
    numtl = CCTK_ActiveTimeLevelsVI(cctkGH, rhsindex);
    if(numtl <= 0) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "No storage for %s. varindex = %d, rhsindex = %d, numtl = %d",
                  "szrhs", varindex, rhsindex, numtl);
    }
  }

  if(do_tau) {
    if(!CCTK_IsImplementationActive("Spritz")) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Please enable a thorn implementing Spritz when using do_tau. While this *might* work without Spritz if another prim2con is provided, I am not yet aware of sucha thorn");
    }
  }

  if(CCTK_IsImplementationActive("Spritz")) {
    if(!compute_DCandMaxima){
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Please set Spritz::compute_DCandMaxima = yes, in order to correctly compute the global maximum of rho");
    }
  }

  // TODO: instead consider adding STORAGE statements to schedule.ccl?
  if(CCTK_ActiveTimeLevels(cctkGH, "HydroBase::entropy") == 0 ||
     CCTK_ActiveTimeLevels(cctkGH, "HydroBase::temperature") == 0) {
    CCTK_ERROR("Require storage for entropy and temperature. Please make sure to set HydroBase::initial_entropy and HydroBase::initial_temperature to the correct values");
  }

}
