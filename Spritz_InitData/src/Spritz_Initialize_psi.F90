 /*@@
   @file      Spritz_Initialize_psi
   @date
   @author    Bruno Giacomazzo
   @desc 
   Initialize psi_mhd to zero
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

 /*@@
   @file      Spritz_Initialize_psi
   @date
   @author    Bruno Giacomazzo
   @desc 
   Initialize psi_mhd to zero
   @enddesc 
@@*/

subroutine Spritz_Initialize_psi(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS


  CCTK_INT::nx,ny,nz,i,j,k


  call CCTK_INFO("Initialize psi_mhd to zero")

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  do i=1,nx
     do j=1,ny
        do k=1,nz
           psi_mhd(i,j,k)=0.0d0
           psi_mhd_p(i,j,k)=0.0d0
           if(timelevels==3) then
              psi_mhd_p_p(i,j,k)=0.0d0
           end if
        end do
     end do
  end do

end subroutine Spritz_Initialize_psi
