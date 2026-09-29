#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <fstream>
#include <string.h>

extern "C"

void ID_NrotStar_check_parameters (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (not CCTK_EQUALS (initial_data,    "ID_NrotStar") or
      not CCTK_EQUALS (initial_lapse,   "ID_NrotStar") or
      not CCTK_EQUALS (initial_shift,   "ID_NrotStar") or
      not (CCTK_EQUALS (initial_dtlapse, "ID_NrotStar") or
           CCTK_EQUALS (initial_dtlapse, "none")) or
      not (CCTK_EQUALS (initial_dtshift, "ID_NrotStar") or
           CCTK_EQUALS (initial_dtshift, "none")) or
      not CCTK_EQUALS (initial_hydro,   "ID_NrotStar"))
  {
    CCTK_PARAMWARN ("The parameters ADMBase::initial_data, ADMBase::initial_lapse, ADMBase::initial_shift, ADMBase::initial_dtlapse, ADMBase::initial_dtshift, HydroBase::initial_hydro, and HydroBase::initial_Bvec must all be set to the value \"ID_NrotStar\" or \"none\"");
  }

  if(not CCTK_EQUALS (temperature_evolution_method, "none") and
     not CCTK_IsThornActive("Spritz_SetBeta")){
    CCTK_PARAMWARN ("If you wanted to evolve temperature, you should activate Spritz_SetBeta to set up initial beta-equilibrated data");
  }
}
