  /*@@
   @file      Spritz_CalcUpdate.F90
   @date      Thu Jan  11 11:03:32 2002
   @author    Ian Hawke
   @desc
   Calculates the update terms given the fluxes. Moved to here so that
   @enddesc
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"


 /*@@
   @routine    UpdateCalculation
   @date       Wed Feb 13 11:03:32 2002
   @author     Ian Hawke
   @desc
   Calculates the update terms from the fluxes.
   @enddesc
   @calls
   @calledby
   @history
   Moved out of the Riemann solver routines to make the FishEye /
   weighted flux calculation easier.
   But we commented out the FishEye code, since it was not used.
   @endhistory

@@*/


subroutine UpdateCalculation(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i,j,k
  CCTK_REAL :: idx, alp_l, alp_r, psi4_loc, det, sqrtdet,  uxx, uxy, uxz, uyy, uyz, uzz


  idx = 1.d0 / CCTK_DELTA_SPACE(flux_direction)


        !$OMP PARALLEL DO PRIVATE(k, j, i, alp_l, alp_r, psi4_loc, &
        !$OMP det, uxx, uxy, uxz, uyy, uyz, uzz)
        do k = Spritz_stencil + 1, cctk_lsh(3) - Spritz_stencil
           do j = Spritz_stencil + 1, cctk_lsh(2) - Spritz_stencil
              do i = Spritz_stencil + 1, cctk_lsh(1) - Spritz_stencil
						
                 if (high_order .ne. 0) then
!!$ If high-order correction is used, the fluxes have been already multiplied by alpha in 
!!$ Spritz_HLLE_MHD.F90, because the correction must take into account the complete values of the fluxes
                 		alp_l = 1.0d0
                 		alp_r = 1.0d0
                 else
                 		alp_l = 0.5d0 * (alp(i,j,k) + &
		                   alp(i-xoffset,j-yoffset,k-zoffset))
		               alp_r = 0.5d0 * (alp(i,j,k) + &
		                   alp(i+xoffset,j+yoffset,k+zoffset))
                 end if
		              	

                 densrhs(i,j,k) = densrhs(i,j,k) + &
                      (alp_l * densflux(i-xoffset,j-yoffset,k-zoffset) - &
                      alp_r * densflux(i,j,k)) * idx
                 srhs(i,j,k,1) = srhs(i,j,k,1) + &
                      (alp_l * sxflux(i-xoffset,j-yoffset,k-zoffset) - &
                      alp_r * sxflux(i,j,k)) * idx
                 srhs(i,j,k,2) = srhs(i,j,k,2) + &
                      (alp_l * syflux(i-xoffset,j-yoffset,k-zoffset) - &
                      alp_r * syflux(i,j,k)) * idx
                 srhs(i,j,k,3) = srhs(i,j,k,3) + &
                      (alp_l * szflux(i-xoffset,j-yoffset,k-zoffset) - &
                      alp_r * szflux(i,j,k)) * idx
                 taurhs(i,j,k) = taurhs(i,j,k) + &
                      (alp_l * tauflux(i-xoffset,j-yoffset,k-zoffset) - &
                      alp_r * tauflux(i,j,k)) * idx

                 Bxrhs(i,j,k) = Bxrhs(i,j,k) + &
                      (alp_l * Bxflux(i-xoffset,j-yoffset,k-zoffset) - &
                      alp_r * Bxflux(i,j,k)) * idx
                 Byrhs(i,j,k) = Byrhs(i,j,k) + &
                      (alp_l * Byflux(i-xoffset,j-yoffset,k-zoffset) - &
                      alp_r * Byflux(i,j,k)) * idx
                 Bzrhs(i,j,k) = Bzrhs(i,j,k) + &
                      (alp_l * Bzflux(i-xoffset,j-yoffset,k-zoffset) - &
                      alp_r * Bzflux(i,j,k)) * idx

                 if (evolve_Y_e .ne. 0) then
                 !!$ Compute alp_l and alp_r if they were not computed before
!                 	  if (high_order .ne. 0) then
!                 	  		alp_l = 0.5d0 * (alp(i,j,k) + &
!		                   alp(i-xoffset,j-yoffset,k-zoffset))
!		               	alp_r = 0.5d0 * (alp(i,j,k) + &
!		                   alp(i+xoffset,j+yoffset,k+zoffset))
!                 	  end if
                    Y_e_con_rhs(i,j,k) = Y_e_con_rhs(i,j,k) + &
                         (alp_l * Y_e_con_flux(i-xoffset,j-yoffset,k-zoffset) - &
                         alp_r * Y_e_con_flux(i,j,k)) * idx
                 end if


              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO

  return

end subroutine UpdateCalculation
