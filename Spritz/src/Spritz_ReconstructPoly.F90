 /*@@
   @file      Spritz_ReconstructPoly.F90
   @date      Sat Jan 26 02:13:25 2002
   @author
   @desc
   Wrapper routine to perform the reconstruction for polytropes.
   @enddesc
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#define velx(i,j,k) vel(i,j,k,1)
#define vely(i,j,k) vel(i,j,k,2)
#define velz(i,j,k) vel(i,j,k,3)
#define sx(i,j,k) scon(i,j,k,1)
#define sy(i,j,k) scon(i,j,k,2)
#define sz(i,j,k) scon(i,j,k,3)

 /*@@
   @routine    ReconstructionPolytype
   @date       Tue Mar 19 11:36:55 2002
   @author     Ian Hawke
   @desc
   If using a polytropic  EOS, we do not have to reconstruct all the 
   variables.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/

subroutine ReconstructionPolytype(CCTK_ARGUMENTS)
  use Interfaces, only: SimplePPM_1d_MHD, wenoz_reconstruct_1D
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: nx, ny, nz, i, j, k

  logical, dimension(:,:,:), allocatable :: trivial_rp
!!$  logical, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: trivial_rp

  CCTK_REAL, dimension(:,:,:),allocatable :: &
       &psi4, lbetax, lbetay, lbetaz
!!$  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: &
!!$       &psi4, lbetax, lbetay, lbetaz

  CCTK_INT :: ierr, mhd_stc

!!$  CCTK_REAL :: v2p, v2m
!!$  character(len=200) warnline
!!$  CCTK_REAL :: gxxl,gxyl,gxzl,gyyl,gyzl,gzzl
!!$  CCTK_REAL :: gxxr,gxyr,gxzr,gyyr,gyzr,gzzr


  allocate(trivial_rp(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)),STAT=ierr)
  if (ierr .ne. 0) then
    call CCTK_WARN(0, "Allocation problems with trivial_rp")
  end if

  allocate(psi4(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)),&
       lbetax(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)),&
       lbetay(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)),&
       lbetaz(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)),STAT=ierr)
  if (ierr .ne. 0) then
    call CCTK_WARN(0, "Allocation problems with lbeta")
  end if


!we need a larger stencil, actually because of Constrained transport and high-order methods
  if (high_order .ne. 0) then
  		if (approximation_order == 4) mhd_stc = 1
  		if (approximation_order == 6) mhd_stc = 2
  else
  		mhd_stc = 1
  end if


  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  trivial_rp = .false.

  if (conformal_state .eq. 0) then
    psi4 = 1.d0
  else
    psi4 = psi**4
  end if
  if (shift_state .ne. 0) then
    lbetax = betax
    lbetay = betay
    lbetaz = betaz
  else
    lbetax = 0.d0
    lbetay = 0.d0
    lbetaz = 0.d0
  end if

	if (CCTK_EQUALS(recon_method, "ppm")) then
	  if (flux_direction == 1) then
		  !$OMP PARALLEL DO PRIVATE(k, j, i)
		  do k = Spritz_stencil, nz - Spritz_stencil + 1
			  do j = Spritz_stencil, ny - Spritz_stencil + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! The routine SimplePPM_1d() expects to get passed 1D arrays (which are laid out
! contiguous in memory) as parameters.
! Thus passing a Fortran90 array section from a 3D grid function rho
! into SimplePPM_1d() -- as for example in rho(j,:,k) --
! requires the compiler for the y and z dimensions to
!
!   (1) allocate a temporary 1D array and copy the values from rho into it
!   (2) call SimplePPM_1d() with that contiguous 1D array
!   (3) copy the values back into rho and release the 1D array
!
! While it is convenient for the programmer to use Fortran90 array sections,
! they become expensive when frequently used.
! It is more efficient to avoid them and let the called routine handle
! the 1D slicing within a 3D grid function itself.
! For this purpose we pass here the address of the first array element
! of the corresponding slice into SimplePPM_1d(), as well as the number
! of points in that slice (np) and the offset between points (di=1 for
! an x slice, di=nx for a y slice, di=nx*ny for a z slice).
! In SimplePPM_1d() we then declare a 1D slice array parameter as
!
!   CCTK_REAL, dimension(di, np) :: rho
!
! and use it as rho(1,i) where i runs from 1..np.
! Note that using the constant 1 as the first subscript entirely leaves the
! index calculations to the compiler (which we hope is proabably more efficient
! than any handwritten indexing scheme).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			     call SimplePPM_1d_MHD(Spritz_eos_handle,1,nx,CCTK_DELTA_SPACE(1),&
			          rho(:,j,k),velx(:,j,k),vely(:,j,k),velz(:,j,k),eps(:,j,k),&
			          press(:,j,k),Bx(:,j,k),By(:,j,k),Bz(:,j,k),&
			          rhominus(:,j,k),velxminus(:,j,k),velyminus(:,j,k),&
			          velzminus(:,j,k),epsminus(:,j,k),&
			          Bxminus(:,j,k),Byminus(:,j,k),Bzminus(:,j,k),rhoplus(:,j,k),&
			          velxplus(:,j,k),velyplus(:,j,k),velzplus(:,j,k),epsplus(:,j,k),&
			          Bxplus(:,j,k),Byplus(:,j,k),Bzplus(:,j,k),&
			          trivial_rp(:,j,k), &
			          gxx(:,j,k), gxy(:,j,k), gxz(:,j,k), gyy(:,j,k), gyz(:,j,k), &
			          gzz(:,j,k), psi4(:,j,k), lbetax(:,j,k), alp(:,j,k),&
			          w_lorentz(:,j,k), &
			          1, j, k, nx, ny, nz, Spritz_Cs2_min )

			  end do
		  end do
		  !$OMP END PARALLEL DO
	  else if (flux_direction == 2) then
		  !$OMP PARALLEL DO PRIVATE(k, j, i)
		  do k = Spritz_stencil, nz - Spritz_stencil + 1
			  do j = Spritz_stencil, nx - Spritz_stencil + 1
			     call SimplePPM_1d_MHD(Spritz_eos_handle,1,ny,CCTK_DELTA_SPACE(2),&
			          rho(j,:,k),vely(j,:,k),velz(j,:,k),velx(j,:,k),eps(j,:,k),&
			          press(j,:,k),By(j,:,k),Bz(j,:,k),Bx(j,:,k),&
			          rhominus(j,:,k),velyminus(j,:,k),velzminus(j,:,k),&
			          velxminus(j,:,k),epsminus(j,:,k),&
			          Byminus(j,:,k),Bzminus(j,:,k),Bxminus(j,:,k),rhoplus(j,:,k),&
			          velyplus(j,:,k),velzplus(j,:,k),velxplus(j,:,k),epsplus(j,:,k),&
			          Byplus(j,:,k),Bzplus(j,:,k),Bxplus(j,:,k),&
			          trivial_rp(j,:,k), &
			          gyy(j,:,k), gyz(j,:,k), gxy(j,:,k), gzz(j,:,k), gxz(j,:,k), &
			          gxx(j,:,k), psi4(j,:,k), lbetay(j,:,k), alp(j,:,k),&
			          w_lorentz(j,:,k), &
			          2, j, k, nx, ny, nz, Spritz_Cs2_min )

			  end do
		  end do
		  !$OMP END PARALLEL DO
	  else if (flux_direction == 3) then
		  !$OMP PARALLEL DO PRIVATE(k, j, i)
		  do k = Spritz_stencil, ny - Spritz_stencil + 1
			  do j = Spritz_stencil, nx - Spritz_stencil + 1
			     call SimplePPM_1d_MHD(Spritz_eos_handle,1,nz,CCTK_DELTA_SPACE(3),&
			          rho(j,k,:),velz(j,k,:),velx(j,k,:),vely(j,k,:),eps(j,k,:),&
			          press(j,k,:),Bz(j,k,:),Bx(j,k,:),By(j,k,:),&
			          rhominus(j,k,:),velzminus(j,k,:),velxminus(j,k,:),&
			          velyminus(j,k,:),epsminus(j,k,:),&
			          Bzminus(j,k,:),Bxminus(j,k,:),Byminus(j,k,:),rhoplus(j,k,:),&
			          velzplus(j,k,:),velxplus(j,k,:),velyplus(j,k,:),epsplus(j,k,:),&
			          Bzplus(j,k,:),Bxplus(j,k,:),Byplus(j,k,:),&
			          trivial_rp(j,k,:), &
			          gzz(j,k,:), gxz(j,k,:), gyz(j,k,:), gxx(j,k,:), gxy(j,k,:), &
			          gyy(j,k,:), psi4(j,k,:), lbetaz(j,k,:), alp(j,k,:),&
			          w_lorentz(j,k,:), &
			          3, j, k, nx, ny, nz, Spritz_Cs2_min )


			  end do
		  end do
		  !$OMP END PARALLEL DO
	  else
		  !$OMP CRITICAL
		  call CCTK_WARN(0, "Flux direction not x,y,z")
		  !$OMP END CRITICAL
	  end if
	
	else if (CCTK_EQUALS(recon_method, "weno-z")) then
		
		if (flux_direction == 1) then
       !$OMP PARALLEL DO PRIVATE(i, j, k)
       do k = Spritz_stencil - mhd_stc, nz - Spritz_stencil + 1 + mhd_stc 
          do j = Spritz_stencil - mhd_stc, ny - Spritz_stencil + 1 + mhd_stc 
             call wenoz_reconstruct_1D(nx, Bx(:,j,k), Bxminus(:,j,k), Bxplus(:,j,k), trivial_rp(:,j,k))
             call wenoz_reconstruct_1D(nx, By(:,j,k), Byminus(:,j,k), Byplus(:,j,k), trivial_rp(:,j,k))
             call wenoz_reconstruct_1D(nx, Bz(:,j,k), Bzminus(:,j,k), Bzplus(:,j,k), trivial_rp(:,j,k))
             call wenoz_reconstruct_1D(nx, rho(:,j,k), rhominus(:,j,k), rhoplus(:,j,k), trivial_rp(:,j,k))
!             call wenoz_reconstruct_1D(nx, eps(:,j,k), epsminus(:,j,k), epsplus(:,j,k), trivial_rp(:,j,k))
             call wenoz_reconstruct_1D(nx, velx(:,j,k), velxminus(:,j,k), velxplus(:,j,k), trivial_rp(:,j,k))
             call wenoz_reconstruct_1D(nx, vely(:,j,k), velyminus(:,j,k), velyplus(:,j,k), trivial_rp(:,j,k))
             call wenoz_reconstruct_1D(nx, velz(:,j,k), velzminus(:,j,k), velzplus(:,j,k), trivial_rp(:,j,k))
          end do
       end do
       !$OMP END PARALLEL DO
  
       
     else if (flux_direction == 2) then
       !$OMP PARALLEL DO PRIVATE(i, j, k)
       do k = Spritz_stencil - mhd_stc, nz - Spritz_stencil + 1 + mhd_stc
          do j = Spritz_stencil - mhd_stc, nx - Spritz_stencil + 1 + mhd_stc
             call wenoz_reconstruct_1D(ny, Bx(j,:,k), Bxminus(j,:,k), Bxplus(j,:,k), trivial_rp(j,:,k))
             call wenoz_reconstruct_1D(ny, By(j,:,k), Byminus(j,:,k), Byplus(j,:,k), trivial_rp(j,:,k))
             call wenoz_reconstruct_1D(ny, Bz(j,:,k), Bzminus(j,:,k), Bzplus(j,:,k), trivial_rp(j,:,k))
             call wenoz_reconstruct_1D(ny, rho(j,:,k), rhominus(j,:,k), rhoplus(j,:,k), trivial_rp(j,:,k))
!             call wenoz_reconstruct_1D(ny, eps(j,:,k), epsminus(j,:,k), epsplus(j,:,k), trivial_rp(j,:,k))
             call wenoz_reconstruct_1D(ny, velx(j,:,k), velxminus(j,:,k), velxplus(j,:,k), trivial_rp(j,:,k))
             call wenoz_reconstruct_1D(ny, vely(j,:,k), velyminus(j,:,k), velyplus(j,:,k), trivial_rp(j,:,k))
             call wenoz_reconstruct_1D(ny, velz(j,:,k), velzminus(j,:,k), velzplus(j,:,k), trivial_rp(j,:,k))
          end do
       end do
       !$OMP END PARALLEL DO
          
       
     else if (flux_direction == 3) then
       !$OMP PARALLEL DO PRIVATE(i, j, k)
       do k = Spritz_stencil - mhd_stc, ny - Spritz_stencil + 1 + mhd_stc
          do j = Spritz_stencil - mhd_stc, nx - Spritz_stencil + 1 + mhd_stc
             call wenoz_reconstruct_1D(nz, Bx(j,k,:), Bxminus(j,k,:), Bxplus(j,k,:), trivial_rp(j,k,:))
             call wenoz_reconstruct_1D(nz, By(j,k,:), Byminus(j,k,:), Byplus(j,k,:), trivial_rp(j,k,:))
             call wenoz_reconstruct_1D(nz, Bz(j,k,:), Bzminus(j,k,:), Bzplus(j,k,:), trivial_rp(j,k,:))
             call wenoz_reconstruct_1D(nz, rho(j,k,:), rhominus(j,k,:), rhoplus(j,k,:), trivial_rp(j,k,:))
!             call wenoz_reconstruct_1D(nz, eps(j,k,:), epsminus(j,k,:), epsplus(j,k,:), trivial_rp(j,k,:))
             call wenoz_reconstruct_1D(nz, velx(j,k,:), velxminus(j,k,:), velxplus(j,k,:), trivial_rp(j,k,:))
             call wenoz_reconstruct_1D(nz, vely(j,k,:), velyminus(j,k,:), velyplus(j,k,:), trivial_rp(j,k,:))
             call wenoz_reconstruct_1D(nz, velz(j,k,:), velzminus(j,k,:), velzplus(j,k,:), trivial_rp(j,k,:))
          end do
       end do
       !$OMP END PARALLEL DO
       
       
     else
        !$OMP CRITICAL
        call CCTK_WARN(0, "Flux direction not x,y,z")
        !$OMP END CRITICAL
     end if

  else
     !$OMP CRITICAL
     call CCTK_WARN(0, "Reconstruction method not recognised.")
     !$OMP END CRITICAL

  end if

  deallocate(psi4, lbetax, lbetay, lbetaz)

  call Prim2ConservativePolytype(CCTK_PASS_FTOF)

  return

end subroutine ReconstructionPolytype
