 /*@@
   @file      Spritz_Reconstruct.F90
   @date      Sat Jan 26 02:13:25 2002
   @author
   @desc
   Wrapper routine to perform the reconstruction.
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
   @routine    Reconstruction
   @date       Sat Jan 26 02:13:47 2002
   @author     Luca Baiotti, Ian Hawke
   @desc
   A wrapper routine to do reconstruction. Currently just does
   TVD on the primitive variables.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/

subroutine Reconstruction(CCTK_ARGUMENTS)

  use Interfaces, only: SimplePPM_1d_MHD, wenoz_reconstruct_1D
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: nx, ny, nz, i, j, k, itracer

  logical, dimension(:,:,:), allocatable :: trivial_rp

  CCTK_INT :: mhd_stc
  CCTK_REAL, dimension(:,:,:),allocatable :: &
       &psi4, lbetax, lbetay, lbetaz
  CCTK_INT :: ierr
  CCTK_REAL :: local_min_tracer

  !we need a larger stencil, actually because of Constrained transport and high-order methods
  if (high_order .ne. 0) then
  		if (approximation_order == 4) mhd_stc = 1
  		if (approximation_order == 6) mhd_stc = 2
  else
  		mhd_stc = 1
  end if
  

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




  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

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

!!$ Initialize variables that store reconstructed quantities

  rhoplus = 0.0d0
  rhominus = 0.0d0
  epsplus = 0.0d0
  epsminus = 0.0d0
  velxplus = 0.0d0
  velxminus = 0.0d0
  velyplus = 0.0d0
  velyminus = 0.0d0
  velzplus = 0.0d0
  velzminus = 0.0d0
  Bxplus = 0.0d0
  Bxminus = 0.0d0
  Byplus = 0.0d0
  Byminus = 0.0d0
  Bzplus = 0.0d0
  Bzminus = 0.0d0


  if (CCTK_EQUALS(recon_method, "ppm")) then

     if (flux_direction == 1) then
        !$OMP PARALLEL DO PRIVATE(i, j, k)
        do k = Spritz_stencil - mhd_stc, nz - Spritz_stencil + 1 + mhd_stc
           do j = Spritz_stencil - mhd_stc, ny - Spritz_stencil + 1 + mhd_stc
              call SimplePPM_1d_MHD(Spritz_eos_handle,0,nx,CCTK_DELTA_SPACE(1),&
                   rho(:,j,k),velx(:,j,k),vely(:,j,k),velz(:,j,k),eps(:,j,k),&
                   press(:,j,k),Bx(:,j,k),By(:,j,k),Bz(:,j,k),&
                   rhominus(:,j,k),velxminus(:,j,k),velyminus(:,j,k),&
                   velzminus(:,j,k),epsminus(:,j,k),&
                   Bxminus(:,j,k),Byminus(:,j,k),Bzminus(:,j,k),rhoplus(:,j,k),&
                   velxplus(:,j,k),velyplus(:,j,k),velzplus(:,j,k),epsplus(:,j,k),&
                   Bxplus(:,j,k),Byplus(:,j,k),Bzplus(:,j,k),&
                   trivial_rp(:,j,k),&
                   gxx(:,j,k), gxy(:,j,k), gxz(:,j,k), gyy(:,j,k), gyz(:,j,k), &
                   gzz(:,j,k), psi4(:,j,k), lbetax(:,j,k), alp(:,j,k),&
                   w_lorentz(:,j,k), &
                   1, j, k, nx, ny, nz, Spritz_Cs2_min )

           end do
        end do
        !$OMP END PARALLEL DO

    if(evolve_Y_e.ne.0) then
        !$OMP PARALLEL DO PRIVATE(j, k)
         do k = Spritz_stencil, nz - Spritz_stencil + 1
           do j = Spritz_stencil, ny - Spritz_stencil + 1
                call SimplePPM_ye_1d(nx,CCTK_DELTA_SPACE(1),rho(:,j,k), &
                     vel(:,j,k,1),vel(:,j,k,2),vel(:,j,k,3), &
                     Y_e(:,j,k),Y_e_minus(:,j,k),Y_e_plus(:,j,k), &
                     press(:,j,k))
           end do
        end do
    end if

     else if (flux_direction == 2) then
        !$OMP PARALLEL DO PRIVATE(i, j, k)
        do k = Spritz_stencil - mhd_stc, nz - Spritz_stencil + 1 + mhd_stc
           do j = Spritz_stencil - mhd_stc, nx - Spritz_stencil + 1 + mhd_stc
              call SimplePPM_1d_MHD(Spritz_eos_handle,0,ny,CCTK_DELTA_SPACE(2),&
                   rho(j,:,k),vely(j,:,k),velz(j,:,k),velx(j,:,k),eps(j,:,k),&
                   press(j,:,k),By(j,:,k),Bz(j,:,k),Bx(j,:,k),&
                   rhominus(j,:,k),velyminus(j,:,k),velzminus(j,:,k),&
                   velxminus(j,:,k),epsminus(j,:,k),&
                   Byminus(j,:,k),Bzminus(j,:,k),Bxminus(j,:,k),rhoplus(j,:,k),&
                   velyplus(j,:,k),velzplus(j,:,k),velxplus(j,:,k),epsplus(j,:,k),&
                   Byplus(j,:,k),Bzplus(j,:,k),Bxplus(j,:,k),&
                   trivial_rp(j,:,k),&
                   gyy(j,:,k), gyz(j,:,k), gxy(j,:,k), gzz(j,:,k), gxz(j,:,k), &
                   gxx(j,:,k), psi4(j,:,k), lbetay(j,:,k), alp(j,:,k),&
                   w_lorentz(j,:,k), &
                   2, j, k, nx, ny, nz, Spritz_Cs2_min )
           end do
        end do
        !$OMP END PARALLEL DO

      if(evolve_Y_e.ne.0) then
         !$OMP PARALLEL DO PRIVATE(j, k)
          do k = Spritz_stencil, nz - Spritz_stencil + 1
            do j = Spritz_stencil, nx - Spritz_stencil + 1
                 call SimplePPM_ye_1d(ny,CCTK_DELTA_SPACE(2),rho(j,:,k), &
                      vel(j,:,k,2),vel(j,:,k,3),vel(j,:,k,1), &
                      Y_e(j,:,k),Y_e_minus(j,:,k),Y_e_plus(j,:,k), &
                      press(j,:,k))
            end do
         end do
         !$OMP END PARALLEL DO
      end if

     else if (flux_direction == 3) then
        !$OMP PARALLEL DO PRIVATE(i, j, k)
        do k = Spritz_stencil - mhd_stc, ny - Spritz_stencil + 1 + mhd_stc
           do j = Spritz_stencil - mhd_stc, nx - Spritz_stencil + 1 + mhd_stc
              call SimplePPM_1d_MHD(Spritz_eos_handle,0,nz,CCTK_DELTA_SPACE(3),&
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

      if(evolve_Y_e.ne.0) then
         !$OMP PARALLEL DO PRIVATE(j, k)
          do k = Spritz_stencil, ny - Spritz_stencil + 1
            do j = Spritz_stencil, nx - Spritz_stencil + 1
                 call SimplePPM_ye_1d(nz,CCTK_DELTA_SPACE(3),rho(j,k,:), &
                      vel(j,k,:,3),vel(j,k,:,1),vel(j,k,:,2), &
                      Y_e(j,k,:),Y_e_minus(j,k,:),Y_e_plus(j,k,:), &
                      press(j,k,:))
            end do
         end do
         !$OMP END PARALLEL DO
      end if

     else
        !$OMP CRITICAL
        call CCTK_WARN(0, "Flux direction not x,y,z")
        !$OMP END CRITICAL
     end if

  else if (CCTK_EQUALS(recon_method, "minmod")) then

     call minmod_reconstruct_mhd(nx, ny, nz, xoffset, yoffset, zoffset, &
          rho, rhoplus, rhominus, trivial_rp)
     call minmod_reconstruct_mhd(nx, ny, nz, xoffset, yoffset, zoffset, &
          vel(:,:,:,1), velxplus, velxminus, trivial_rp)
     call minmod_reconstruct_mhd(nx, ny, nz, xoffset, yoffset, zoffset, &
          vel(:,:,:,2), velyplus, velyminus, trivial_rp)
     call minmod_reconstruct_mhd(nx, ny, nz, xoffset, yoffset, zoffset, &
          vel(:,:,:,3), velzplus, velzminus, trivial_rp)
     call minmod_reconstruct_mhd(nx, ny, nz, xoffset, yoffset, zoffset, &
          eps, epsplus, epsminus, trivial_rp)
     call minmod_reconstruct_mhd(nx, ny, nz, xoffset, yoffset, zoffset, &
          Bx, Bxplus, Bxminus, trivial_rp)
     call minmod_reconstruct_mhd(nx, ny, nz, xoffset, yoffset, zoffset, &
          By, Byplus, Byminus, trivial_rp)
     call minmod_reconstruct_mhd(nx, ny, nz, xoffset, yoffset, zoffset, &
          Bz, Bzplus, Bzminus, trivial_rp)
     if(evolve_Y_e.ne.0) then
       call minmod_reconstruct_mhd(nx, ny, nz, xoffset, yoffset, zoffset, &
            Y_e, Y_e_plus, Y_e_minus, trivial_rp)
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
             call wenoz_reconstruct_1D(nx, eps(:,j,k), epsminus(:,j,k), epsplus(:,j,k), trivial_rp(:,j,k))
             call wenoz_reconstruct_1D(nx, velx(:,j,k), velxminus(:,j,k), velxplus(:,j,k), trivial_rp(:,j,k))
             call wenoz_reconstruct_1D(nx, vely(:,j,k), velyminus(:,j,k), velyplus(:,j,k), trivial_rp(:,j,k))
             call wenoz_reconstruct_1D(nx, velz(:,j,k), velzminus(:,j,k), velzplus(:,j,k), trivial_rp(:,j,k))
             if(evolve_Y_e.ne.0) then
       			call wenoz_reconstruct_1D(nx, Y_e(:,j,k), Y_e_minus(:,j,k), Y_e_plus(:,j,k), trivial_rp(:,j,k))
       		 end if
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
             call wenoz_reconstruct_1D(ny, eps(j,:,k), epsminus(j,:,k), epsplus(j,:,k), trivial_rp(j,:,k))
             call wenoz_reconstruct_1D(ny, velx(j,:,k), velxminus(j,:,k), velxplus(j,:,k), trivial_rp(j,:,k))
             call wenoz_reconstruct_1D(ny, vely(j,:,k), velyminus(j,:,k), velyplus(j,:,k), trivial_rp(j,:,k))
             call wenoz_reconstruct_1D(ny, velz(j,:,k), velzminus(j,:,k), velzplus(j,:,k), trivial_rp(j,:,k))
             if(evolve_Y_e.ne.0) then
       			call wenoz_reconstruct_1D(ny, Y_e(j,:,k), Y_e_minus(j,:,k), Y_e_plus(j,:,k), trivial_rp(j,:,k))  		
       		 end if
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
             call wenoz_reconstruct_1D(nz, eps(j,k,:), epsminus(j,k,:), epsplus(j,k,:), trivial_rp(j,k,:))
             call wenoz_reconstruct_1D(nz, velx(j,k,:), velxminus(j,k,:), velxplus(j,k,:), trivial_rp(j,k,:))
             call wenoz_reconstruct_1D(nz, vely(j,k,:), velyminus(j,k,:), velyplus(j,k,:), trivial_rp(j,k,:))
             call wenoz_reconstruct_1D(nz, velz(j,k,:), velzminus(j,k,:), velzplus(j,k,:), trivial_rp(j,k,:))
             if(evolve_Y_e.ne.0) then
             	call wenoz_reconstruct_1D(nz, Y_e(j,k,:), Y_e_minus(j,k,:), Y_e_plus(j,k,:), trivial_rp(j,k,:))
          	 end if
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


  !$OMP WORKSHARE
  where ( (rhoplus < Spritz_rho_min).or.(rhominus < Spritz_rho_min).or.(epsplus < 0.d0).or.(epsminus < 0.d0) )
     rhoplus = rho
     rhominus = rho
     velxplus = vel(:,:,:,1)
     velxminus = vel(:,:,:,1)
     velyplus = vel(:,:,:,2)
     velyminus = vel(:,:,:,2)
     velzplus = vel(:,:,:,3)
     velzminus = vel(:,:,:,3)
     epsplus = eps
     epsminus = eps
     Bxplus = Bx
     Bxminus = Bx
     Byplus = By
     Byminus = By
     Bzplus = Bz
     Bzminus = Bz
  end where
  !$OMP END WORKSHARE



  call primitive2conservative(CCTK_PASS_FTOF)

  return

end subroutine Reconstruction
