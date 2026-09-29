// register.cc -- register variables with various thorns that need-to-know
// $Header$
//
// Spritz_register_GZPatchSystem - register with GZPatchSystem
//
// Cut 'n paste job from BackgroundWaveToy, author J Thornburg...

#include <cstdio>
#include <string>

#include "cctk.h"
#include "cctk_Arguments.h"

using namespace std;

//******************************************************************************

//
// This function is called by the Cactus scheduler to (maybe) register
// our to-be-interpatch-synchronized variables with GZPatchSystem.  It
// checks if the GZPatchSystem registration function has been provided
// (which it should be if and only if GZPatchSystem is active), and if so,
// does the registration.
//
// If we're using Carpet, this function must be called in meta mode.
//
extern "C"void Spritz_register_GZPatchSystem(CCTK_ARGUMENTS)
{

  if (CCTK_IsFunctionAliased("GZPatchSystem_register_sync"))
  {
    CCTK_VInfo(CCTK_THORNSTRING,
               "registering to-be-interpatch-synchronized variables "
               "with GZPatchSystem");

    string var[8] = {"HydroBase::rho", "HydroBase::press", "HydroBase::eps",
                     "HydroBase::vel",
                     "Spritz::dens", "Spritz::tau", "HydroBase::w_lorentz",
                     "Spritz::scon"};
    for (int i = 0; i < 8; i++)
    {
      int status = 
        GZPatchSystem_register_sync(var[i].c_str());
      if (status < 0)
      {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "***** Spritz_register_GZPatchSystem():\n"
                   "        error registering var group %s to be "
                   "interpatch-synchronized!\n"
                   "        (GZPatchSystem_register_sync() error code %d)\n",
                   var[i].c_str(), status);
      }
    }
    
  } 
  else
  {
    CCTK_WARN(1, "Function GZPatchSystem_register_sync not registered!");
  }

  if (CCTK_IsFunctionAliased("GZPatchSystem_register_cxform"))
  {
    CCTK_VInfo(CCTK_THORNSTRING,
               "registering to-be-cxformed variables with GZPatchSystem");

    string var[11] = {"HydroBase::rho", "HydroBase::press", "HydroBase::eps",
                      "HydroBase::vel",
                      "Spritz::dens", "Spritz::tau", "HydroBase::w_lorentz",
                      "Spritz::scon",
                      "ADMBase::metric", "ADMBase::curv", "ADMBase::shift"};
    for (int i = 0; i < 11; i++)
     for (int j = 0; j < 3; j++)
     {
      int ps_status = 
        GZPatchSystem_register_cxform(j, var[i].c_str());
      if (ps_status < 0)
      {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "***** Spritz_register_GZPatchSystem():\n"
                   "        error registering var group %s to be "
                   "interpatch-cxformhronized!\n"
                   "        (GZPatchSystem_register_cxform() error code %d)\n",
                   var[i].c_str(), ps_status);
      }
     }
  }
  else
  {
    CCTK_WARN(1, "Function GZPatchSystem_register_cxform not registered!");
  }
  
}
