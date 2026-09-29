#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

// LSADD

void ZelmaniLeakLS_OutputSettings(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  CCTK_VInfo(CCTK_THORNSTRING,"*********** Leakage Settings ************");

  if(force_postmerger == 1) {
    CCTK_VInfo(CCTK_THORNSTRING,"force_postmerger has been selected: leakage will be forced at time %g.", force_postmerger_at_time);
  }

  return;
}
