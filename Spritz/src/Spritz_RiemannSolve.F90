 /*@@
   @file      Spritz_RiemannSolve.F90
   @date      Sat Jan 26 02:20:25 2002
   @author    
   @desc 
   A wrapper routine to call the correct Riemann solver.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    RiemannSolve
   @date       Sat Jan 26 02:20:48 2002
   @author     Pedro Montero, Ian Hawke
   @desc 
   A wrapper routine to switch between the different Riemann solvers.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine RiemannSolve(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS


  call Spritz_HLLE_MHD(CCTK_PASS_FTOF)
  

    
  
end subroutine RiemannSolve

