 /*@@
   @file      Spritz_Startup.F90
   @date      Sun Feb 10 00:02:52 2002
   @author    Ian Hawke
   @desc 
   Startup banner.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

 /*@@
   @routine    Spritz_Startup
   @date       Sun Feb 10 00:03:09 2002
   @author     Ian Hawke
   @desc 
   Startup banner.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

integer function Spritz_Startup()

  implicit none

  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: ierr

  call CCTK_RegisterBanner(ierr, "Spritz: magnetized flows on mountains")
  
  Spritz_Startup = 0

end function Spritz_Startup
