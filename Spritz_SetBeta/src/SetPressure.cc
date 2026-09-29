/*@@
  @file      SetPressure.cc
  @date
  @author    Lorenzo Sala
  @desc
  Sets up initial pressure, employing the EoS.
  @enddesc
  @@*/

#include <cassert>
#include <iostream>
#include <cstdio>
#include <vector>
#include <ios>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

// headers C
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sstream>

// headers Cactus
#include "cctk_Functions.h"

using namespace std;

extern "C"

void Spritz_SetPressure (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO ("Setting up initial pressure.");
  try {
    int const local_grid_size = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
    CCTK_INT  anyerr;
    CCTK_INT nn,xkeytemp,keyerr[1];
    CCTK_REAL xrho[1],xpress[1],xeps[1],xtemp[1],xye[1];
    nn = 1;
    for (int k=0; k < cctk_lsh[2]; k++) {
      for (int j=0; j < cctk_lsh[1]; j++) {
        for (int i=0; i < cctk_lsh[0]; i++) {
    	  int i3D = CCTK_GFINDEX3D(cctkGH, i, j, k);
          xkeytemp = 0;
          anyerr = 0;
          keyerr[0] = 0;
          xrho[0] = rho[i3D];
          xeps[0] = eps[i3D];
          xtemp[0] = temperature[i3D];
          xye[0] = Y_e[i3D];

          xkeytemp = 1; //temperature has already been set, there is no need of root-finding it again.
          // keytemp=1 means that press is found interpolating in rho,Ye,T table.
          EOS_Omni_press(*Spritz_eos_handle,xkeytemp,Eos_rf_prec,nn,xrho,xeps,xtemp,xye,xpress,keyerr,&anyerr);
          press[i3D] = xpress[0];
        }
      }
    }

    // if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::press") > 1) {
    //   memcpy(press_p, press, local_grid_size*sizeof(CCTK_REAL));
    //   if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::press") > 2) {
    //     memcpy(press_p_p, press, local_grid_size*sizeof(CCTK_REAL));
	//   }
    // }
  CCTK_VINFO ("Pressure set on on reflevel %d.",int(log10(double(cctk_levfac[0]))/log10(2.0)));
  } catch (ios::failure e) {
    CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Could not read initial pressure data.");
  }
}
