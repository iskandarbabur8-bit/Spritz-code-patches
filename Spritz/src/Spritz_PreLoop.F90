 /*@@
   @file      Spritz_PreLoop.F90
   @date      Mon Feb 25 11:43:36 2002
   @author    
   @desc 
   Sets up various scalars used for efficiency reasons.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    Spritz_Scalar_Setup
   @date       Mon Feb 25 11:25:27 2002
   @author     
   @desc 
   Sets up the logical scalars from the parameters. These are
   solely used for efficiency.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Spritz_Scalar_Setup(CCTK_ARGUMENTS)

  USE Spritz_Scalars

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  MINMOD = .false.
  MINMOD2 = .false.
  MINMOD3 = .false.
  MCB = .false.
  MC1 = .false.
  MC2 = .false.
  SUPERBEE = .false.
  VANLEERBG = .false.


  PPM3=.false.
  PPM4=.false.

  if (CCTK_EQUALS(ppm_flatten,"stencil_3")) then
    PPM3=.true.
  else if (CCTK_EQUALS(ppm_flatten,"stencil_4")) then
    PPM4=.true.
  else
    call CCTK_WARN(0, "PPM Flattening Procedure not recognized!")
  end if

  ANALYTICAL = .false.


  FAST = .false.


end subroutine Spritz_Scalar_Setup
