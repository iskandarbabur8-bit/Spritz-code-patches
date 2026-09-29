 /*@@
   @file      Spritz_Maxpress.cc
   @date      Tue Aug 29 18:52:10 2006
   @author    
   @desc 
        Find the maximum of the pressure
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "carpet.hh"

#ifdef __cplusplus
  extern "C" {
#endif
    
    /* Scheduled functions */
    void Spritz_Maxpress_Setup_Final(CCTK_ARGUMENTS);

    void Spritz_Maxpress_Setup_Final_PUGH(CCTK_ARGUMENTS);
    
#ifdef __cplusplus
  } /* extern "C" */
#endif

 /*@@
   @routine    Spritz_MaxPress_Setup_Final
   @date       Tue Aug 29 18:54:05 2006
   @author     Bruno Giacomazzo
   @desc 
        Find the maximum of the pressure to be used to compute the pcut, i.e. the value of the pressure at which the vector potential is zero.
	Version for Carpet.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
@@*/

void Spritz_Maxpress_Setup_Final(CCTK_ARGUMENTS)
{

#ifdef HAVE_CARPET

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL maxpress_tmp;

  static int flag = true;
  
  if (flag) // Actually run this reduction only the first time the routine is called.
    {
      // Go to global mode
      
      //      BEGIN_GLOBAL_MODE (cctkGH) {

	char* infoline;

	infoline = (char *)malloc(1000*sizeof(char));

	snprintf(infoline,1000,"Finding the maximum of the pressure to be used in the definition of the vector potential");
	CCTK_INFO(infoline);

	// Find the global maximum of press
	
	int Reduction_Handle = CCTK_ReductionHandle ("maximum");
	
	CCTK_INT input_array_variable_indices= {CCTK_VarIndex("HydroBase::press")};
	
	int ierr = CCTK_Reduce(cctkGH,
			       -1, // target processors; -1 -> all
			       Reduction_Handle,
			       1,  // number of output variables
			       CCTK_VARIABLE_REAL,
			       &maxpress_tmp,
			       1,  // number of variables to be reduced
			       input_array_variable_indices);
	
	if (ierr != 0)
	  {
	    CCTK_WARN(0, "Failed to compute the global maximum of press");
	  }
	
	*max_press = maxpress_tmp;


	//---------------------------------------------------------------------------------------------------------
	// Find the maximum of press on x>0
	
	maxpress_tmp=0;

	Reduction_Handle = CCTK_ReductionHandle ("maximum");
	
	input_array_variable_indices= CCTK_VarIndex("Spritz_InitData::press_sign");
	
	ierr = CCTK_Reduce(cctkGH,
			   -1, // target processors; -1 -> all
			   Reduction_Handle,
			   1,  // number of output variables
			   CCTK_VARIABLE_REAL,
			   &maxpress_tmp,
			   1,  // number of variables to be reduced
			   input_array_variable_indices);
	
	if (ierr != 0)
	  {
	    CCTK_WARN(0, "Failed to compute the maximum of press_sign");
	  }
	
	*max_press_plus = maxpress_tmp;


	//---------------------------------------------------------------------------------------------------------
	// Find the maximum of press on x<0

	maxpress_tmp=0;
	
	Reduction_Handle = CCTK_ReductionHandle ("minimum");
	
	input_array_variable_indices= CCTK_VarIndex("Spritz_InitData::press_sign");
	
	ierr = CCTK_Reduce(cctkGH,
			   -1, // target processors; -1 -> all
			   Reduction_Handle,
			   1,  // number of output variables
			   CCTK_VARIABLE_REAL,
			   &maxpress_tmp,
			   1,  // number of variables to be reduced
			   input_array_variable_indices);
	
	if (ierr != 0)
	  {
	    CCTK_WARN(0, "Failed to compute the maximum of press_sign");
	  }
	
	*max_press_minus = -1.0*maxpress_tmp;




	// Go back to local mode
	//      } END_GLOBAL_MODE;
      
      // After this has run once, set the flag so that this does not run again
      //flag = false;
    } // end if (flag)
  else
    {
      return;
    }

#endif

}


 /*@@
   @routine    Spritz_MaxPress_Setup_Final_PUGH
   @date       Tue Aug 29 18:57:42 2006
   @author     Bruno Giacomazzo
   @desc 
         As above, but for PUGH.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
@@*/

void Spritz_Maxpress_Setup_Final_PUGH(CCTK_ARGUMENTS)
{ 
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_REAL maxpress_tmp;
  
  // Find the global maximum of press
  
  int Reduction_Handle = CCTK_ReductionHandle("maximum");
  
  CCTK_INT input_array_variable_indices={CCTK_VarIndex("HydroBase::press")};
  
  int ierr = CCTK_Reduce(cctkGH,
			 -1, // target processors; -1 -> all
			 Reduction_Handle,
			 1,  // number of output variables   
			 CCTK_VARIABLE_REAL,
			 &maxpress_tmp,
			 1,  // number of variables to be reduced
			 input_array_variable_indices);
  
  if (ierr != 0)
    {
      CCTK_WARN(0, "Failed to compute the global maximum of press");
    }
  
  *max_press = maxpress_tmp;
  
  
  
  
  // Find the maximum of press on x>0
  
  maxpress_tmp=0;
  
  Reduction_Handle = CCTK_ReductionHandle("maximum");
  
  input_array_variable_indices=CCTK_VarIndex("Spritz_InitData::press_sign");
  
  ierr = CCTK_Reduce(cctkGH,
		     -1, // target processors; -1 -> all
		     Reduction_Handle,
		     1,  // number of output variables   
		     CCTK_VARIABLE_REAL,
		     &maxpress_tmp,
		     1,  // number of variables to be reduced
		     input_array_variable_indices);
  
  if (ierr != 0)
    {
      CCTK_WARN(0, "Failed to compute the global maximum of press");
    }
  
  *max_press_plus = maxpress_tmp;
  
  
  // Find the maximum of press on x<0
  
  maxpress_tmp=0;
  
  Reduction_Handle = CCTK_ReductionHandle("minimum");
  
  input_array_variable_indices=CCTK_VarIndex("Spritz_InitData::press_sign");
  
  ierr = CCTK_Reduce(cctkGH,
		     -1, // target processors; -1 -> all
		     Reduction_Handle,
		     1,  // number of output variables   
		     CCTK_VARIABLE_REAL,
		     &maxpress_tmp,
		     1,  // number of variables to be reduced
		     input_array_variable_indices);
  
  if (ierr != 0)
    {
      CCTK_WARN(0, "Failed to compute the global maximum of press");
    }

  *max_press_minus = -1.0*maxpress_tmp;
    
}
