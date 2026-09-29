 /*@@
   @file      Spritz_EOS.c
   @date      Wed Feb  6 18:25:33 2002
   @author    
   @desc 
   Sets the EOS handle number from CactusEOS for use by
   all the Spritz routines
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include <stdio.h>
#include <stdlib.h>

/*#include "EOS_Base.h"*/

 /*@@
   @routine    Spritz_EOSHandle
   @date       Wed Feb  6 18:28:01 2002
   @author     Ian Hawke
   @desc 
   Sets the EOS handle number
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void Spritz_EOSHandle(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int ierr = 0;
  char *infoline;

  /* Check that the parfile does not activate two EoSs which would conflict */
  if ((CCTK_IsThornActive("EOS_Polytrope"))&&(CCTK_IsThornActive("EOS_HybridPolytrope")))
    {
      CCTK_WARN(0, "Do not activate in the parfile both EOS_Polytrope and EOS_HybridPolytrope.");
    }
  
  /**Spritz_eos_handle = EOS_Handle(Spritz_eos_table);*/
  *Spritz_eos_handle = EOS_Omni_GetHandle(Spritz_eos_table);
  
  if (*Spritz_eos_handle < 0) {
    ierr = *Spritz_eos_handle;
    CCTK_WARN(0, "Table handle is incorrect.");
  }

  /**Spritz_polytrope_handle = EOS_Handle("2D_Polytrope");*/
  *Spritz_polytrope_handle = EOS_Omni_GetHandle("2D_Polytrope");
  
  if (*Spritz_polytrope_handle < 0) {

    /**Spritz_polytrope_handle = EOS_Handle("2D_HybridPolytrope");*/
    *Spritz_polytrope_handle = EOS_Omni_GetHandle("2D_HybridPolytrope");

    if (*Spritz_polytrope_handle < 0) {
      ierr = *Spritz_polytrope_handle;
      CCTK_WARN(0, "Spritz requires the 2D_Polytrope or the 2D_HybridPolytrope EoSs to be registered, for the atmosphere.");
    }
  }

  infoline = (char *)malloc(100*sizeof(char));
  if (!infoline) {
    CCTK_WARN(0, "Failed to malloc an array of 100 chars in Spritz_EOS.");
  }
  else {
    sprintf(infoline,"Spritz will use the %s equation of state.",
            Spritz_eos_table);
    CCTK_INFO(infoline);
  }
  free(infoline);
  infoline = NULL;

}
