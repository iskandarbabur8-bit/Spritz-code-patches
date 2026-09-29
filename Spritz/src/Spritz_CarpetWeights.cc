 /*@@
   @file     Spritz_CarpetWeights.cc
   @date     Wed Aug 20 23:08:28 2008
   @author   Luca Baiotti
   @desc 
   Copy the variable CarpetReduce::weight in a Spritz variable.
   Written with the help of Thomas.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "SpaceMask.h"

#include <stdio.h>
#include <stdlib.h>


#ifdef __cplusplus
  extern "C" {
#endif
    
    /* Scheduled functions */
    void Spritz_CarpetWeights(CCTK_ARGUMENTS);

#ifdef __cplusplus
  } /* extern "C" */
#endif


void Spritz_CarpetWeights(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int i;

  // Define the variable we will use in Con2Prim (F90) as a pointer to itself
  Spritz_CarpetWeights = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,0,"Spritz::Spritz_CarpetWeights"));


  if (CCTK_IsThornActive("CarpetReduce")!=0){
    // If we use Carpet then obtain a pointer to rhe carpet variable. This will be CarpetWeights
    const CCTK_REAL *CarpetWeights = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,0,"CarpetReduce::weight"));
    // Then assign the values of the carpet variable to the Spritz variable
    for (i=0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
      {
	Spritz_CarpetWeights[i] = CarpetWeights[i];
	//printf("%13.12e\n",Spritz_CarpetWeights[i]);
      }
  }
  else{
    // Then assign the values of the carpet variable to the Spritz variable
    CCTK_INFO("We are not using Carpet, so I set Spritz_CarpetWeights equal to 1");
    for (i=0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
      {
	Spritz_CarpetWeights[i] = 1; //if we do not use Carpet then we set this variable to be equal to 1
	//printf("%13.12e\n",Spritz_CarpetWeights[i]);
      }
  }

  return;
}
