 /*@@
   @file      Spritz_MINMOD.F90
   @date      Sat Jan 26 02:11:44 2002
   @author    Luca Baiotti
   @desc 
   The Minmod reconstruction routine.
   Based on Whisky_TVDReconstruct.F90. Modified by Federico Cipolletta and Jay Kalinani. 
   Cleaned different methods and left only minmod. 
   @enddesc 
 @@*/
 

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

  /*@@
   @routine    minmod_reconstruct_mhd
   @date       
   @author     
   @desc 
   Performs minmod reconstruction on the magnetic filed components using linear extrapolation
   @enddesc 
   @calls     
   @calledby   
@@*/


subroutine minmod_reconstruct_mhd(nx, ny, nz, xoffset, yoffset, zoffset, &
     orig, bextp, bextm, trivial_rp) 

  USE Spritz_Scalars

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: i, j, k, xoffset, yoffset, zoffset, nx, ny, nz
  CCTK_REAL, dimension(nx, ny, nz) :: orig, bextp, bextm
  CCTK_REAL :: dupw, dloc, delta, ratio, hdelta
  logical, dimension(nx,ny,nz) :: trivial_rp

  bextp = 0.d0
  bextm = 0.d0

  trivial_rp = .false.

  do k = Spritz_stencil-1, nz-Spritz_stencil+2
     do j = Spritz_stencil-1, ny-Spritz_stencil+2
        do i = Spritz_stencil-1, nx-Spritz_stencil+2

!!$        Calculate the appropriate slopes.

              dupw = orig(i, j, k) - orig(i-xoffset, j-yoffset, k-zoffset)
              dloc = orig(i+xoffset, j+yoffset, k+zoffset) - orig(i, j, k)

!! 	   Call to min mod directly

              delta = minmod_func(dupw,dloc)

              hdelta = 0.5d0 * delta 
              bextm(i, j, k) = orig(i, j, k) - hdelta
              bextp(i, j, k) = orig(i, j, k) + hdelta

        end do
     end do
  end do

contains
  function minmod_func(a_in,b_in) result(minmod_result)
    implicit none
    CCTK_REAL,intent(IN)::a_in,b_in
    CCTK_REAL::minmod_result
    
    minmod_result=0.5D0*(sign(1.0D0,a_in)+sign(1.0D0,b_in))*min(abs(a_in),abs(b_in))
  end function minmod_func

end subroutine minmod_reconstruct_mhd


