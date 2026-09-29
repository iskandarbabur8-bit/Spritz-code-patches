#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "carpet.hh"
#define RHOGF 6.1755e17

// LSADD from CCC

void ZelmaniLeakLS_OutputControl(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  double inv_rho_gf = 1.0/RHOGF;
  if(force_postmerger && (!(*in_merger))) {
    if(cctk_time >= force_postmerger_at_time) {
        *in_merger = 1;
        CCTK_VInfo(CCTK_THORNSTRING,"ZelmaniLeak activating at time %g, as requested by parfile (actually, the requested time was %g)!",double(cctk_time),double(force_postmerger_at_time));
    }
  }

  return;
}
