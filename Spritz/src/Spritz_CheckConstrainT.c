 /*@@
   @file      Spritz_Check_ConstrainT.c
   @date      06-08-2019 16:10
   @author    Lorenzo Sala
   @desc
   Check if temperature is constrained.
   @enddesc
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Spritz_CheckConstrainT(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(*Constrain_Temp_now == 1) {
    if(Spritz_EvolveT_at_time >= 0.0 &&
         cctk_time >= Spritz_EvolveT_at_time &&
         *Spritz_reflevel == 0) {
      CCTK_VINFO("Not constraining temperature anymore, time %f, it %d", cctk_time, cctk_iteration);
      *Constrain_Temp_now = 0;
      int chk = CCTK_ParameterSet ("Spritz_Constrain_Temperature",CCTK_THORNSTRING,"no");
      if(chk != 0) CCTK_INFO("Error in setting Spritz_hot_atmo_Y_e");
    }
  }
  return;
}
