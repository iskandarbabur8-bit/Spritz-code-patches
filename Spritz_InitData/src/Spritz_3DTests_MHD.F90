 /*@@
   @file      Spritz_3DTests_MHD.F90
   @date      Jun 2019
   @author    Federico Cipolletta
   @desc 
   Initial data for MHD tests in 3D.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#define velx(i,j,k) vel(i,j,k,1)
#define vely(i,j,k) vel(i,j,k,2)
#define velz(i,j,k) vel(i,j,k,3)
#define sx(i,j,k) scon(i,j,k,1)
#define sy(i,j,k) scon(i,j,k,2)
#define sz(i,j,k) scon(i,j,k,3)

subroutine Spritz_3DTests_MHD(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k, nx, ny, nz
  CCTK_REAL:: dx, dy, dz
  CCTK_REAL :: det !direction, 
  CCTK_REAL :: rhoin, rhoout, pressin, pressout, epsin, epsout, &
               Bxin, Bxout, Byin, Byout, Bzin, Bzout, & 
               velxin, velxout, velyin, velyout, velzin, velzout, &
               radiusin, radiusout

  CCTK_REAL :: sphexp_prof
  CCTK_REAL :: rtemp

  !!!!!!
  ! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1),one(1)
  !CCTK_REAL :: xeos_gamma(1)
  character(len=500):: warnline
  character(len=500):: infoline
  logical:: condition
  
  n = 1;keytemp = 0;anyerr = 0;keyerr = 0;xrho=0.0d0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
  one = 1.0d0
  ! end EOS Omni vars
  !!!!!!

  psi_mhd = 0.0d0
  psi_mhd_p = 0.0d0
  if (timelevels == 3) then
     psi_mhd_p_p = 0.0d0
  end if

  if (CCTK_EQUALS(initial_shift,"none")) then
     call CCTK_WARN(0,"You have to set ADMBase::initial_shift=zero")
  end if
  
  ! Assign initial parameter, in order to compute initial
  ! data on gridpoints
  if (CCTK_EQUALS(threedtest_case,"SpherExp")) then
     ! See GRHydro paper, Beckwith, Stone (2011) and Komissarov (1999)

     ! Values IN     
     radiusin = R2D_in !0.8d0
     rhoin = rho2D_in !0.01d0
     pressin = P2D_in !1.0d0
     xrho = rhoin
     xpress = pressin
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsin = xeps(1)
     Bxin = Bx2D_init
     Byin = By2D_init
     Bzin = Bz2D_init
     velxin = 0.0d0
     velyin = 0.0d0
     velzin = 0.0d0
     
     !Values OUT
     radiusout = R2D_out !1.0d0 
     rhoout = rho2D_out !0.001d0 !0.0001d0
     pressout = P2D_out !0.0003d0 !0.00003d0
     xrho = rhoout
     xpress = pressout
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsout = xeps(1)
     Bxout = Bx2D_init
     Byout = By2D_init
     Bzout = Bz2D_init
     velxout = 0.0d0
     velyout = 0.0d0
     velzout = 0.0d0
     
     write(infoline,'(a20,g26.15,a11,g26.15,a11,g26.15,a20,g26.15,a11,g26.15,a11,g26.15)') &
                    'CylExp CASE: pin  = ',pressin, ', rhoin  = ',rhoin, ', epsin  = ',epsin,&
                    '             pout = ',pressout,', rhoout = ',rhoout,', epsout = ',epsout
     call CCTK_INFO(infoline)

  else
     call CCTK_WARN(0,"Shock case not recognized")
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Starting Loop for setting initial data at all !
  ! the gridpoints.                               !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  dx = CCTK_DELTA_SPACE(1)	
  dy = CCTK_DELTA_SPACE(2)
  dz = CCTK_DELTA_SPACE(3)

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  do i=1,nx
     do j=1,ny
        do k=1,nz
  
           gxx(i,j,k) = 1.0d0
           gxy(i,j,k) = 0.0d0
           gxz(i,j,k) = 0.0d0
           gyy(i,j,k) = 1.0d0
           gyz(i,j,k) = 0.0d0
           gzz(i,j,k) = 1.0d0
           alp(i,j,k) = alp_value

           if (CCTK_EQUALS(initial_shift,"zero")) then
              betax(i,j,k) = 0.0d0
              betay(i,j,k) = 0.0d0
              betaz(i,j,k) = 0.0d0
           end if

           ! Computing initial data on gridpoints
           if (CCTK_EQUALS(threedtest_case,"SpherExp")) then
              ! See Cerda-Duran et al. (2008) A&A.

              rtemp = sqrt( ( (x(i,j,k) - xc) * (x(i,j,k) - xc) + &
                              (y(i,j,k) - yc) * (y(i,j,k) - yc) + &
                              (z(i,j,k) - zc) * (z(i,j,k) - zc) ) )
              
              rho(i,j,k) = sphexp_prof(rtemp,radiusin,radiusout,rhoin,rhoout)                             
              press(i,j,k) = sphexp_prof(rtemp,radiusin,radiusout,pressin,pressout)

              n = 1
              keytemp = 0
              anyerr = 0 
              keyerr = 0
              xrho = rho(i,j,k)
              xpress = press(i,j,k)
              xtemp = 0.0d0
              xye = 0.0d0
              call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                         xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
              eps(i,j,k) = xeps(1)
              velx(i,j,k) = 0.0d0
              vely(i,j,k) = 0.0d0
              velz(i,j,k) = 0.0d0
              Bx(i,j,k) = Bx2D_init
              By(i,j,k) = By2D_init
              Bz(i,j,k) = Bz2D_init
              dBx(i,j,k) = Bx(i,j,k)
              dBy(i,j,k) = By(i,j,k)
              dBz(i,j,k) = Bz(i,j,k)
              ! The following components of the magnetic vector potential, have been found by trials
              ! but they fulfill the definition of B = curl(Avec)
              ! NOTE: They have not been tested for axis direction different from zdir...

              ! NOTE: there is an addition to half of gridspacing because of staggering!!!
           
              Avec_x(i,j,k) = 0.0d0 !0.5d0 * ( (By(i,j,k)*z(i,j,k)) - (Bz(i,j,k)*y(i,j,k)) )
              Avec_y(i,j,k) = ( Bz(i,j,k)*(x(i,j,k) + (dx/2.0d0)) ) !0.5d0 * ( (Bz(i,j,k)*x(i,j,k)) - (Bx(i,j,k)*z(i,j,k)) )
              Avec_z(i,j,k) = 0.0d0 !0.5d0 * ( (Bx(i,j,k)*y(i,j,k)) - (By(i,j,k)*x(i,j,k)) )

              !if ((abs(x(i,j,k) - 2.0d0) .lt. 1.0d0)) then 
              !   write(warnline,'(a20,g26.15,a11,g26.15,a11,g26.15,a20,g26.15,a11,g26.15,a11,g26.15)') &
              !         'CylExp CASE: p  = ',press(i,j,k), ', rho  = ',rho(i,j,k), ', eps  = ',eps(i,j,k),&
              !         '             Bx = ',Bx(i,j,k),', By = ',By(i,j,k),', Bz = ',Bz(i,j,k)
              !   call CCTK_WARN(1,warnline)
              !end if
           
           else
              call CCTK_WARN(0,"Shock case not recognized")
           end if
           ! End of computation of initial data on gridpoints
        
           call SpatialDet(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),&
                gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det)

           !if ((i .eq. 200) .and. (j .eq. 200) .and. (k .eq. 30)) then
           !   write(infoline,'(a20,g26.15,a11,g26.15,a11,g26.15,a20,g26.15,a11,g26.15,a11,g26.15)') &
           !                  'Before P2C : pres = ',press(i,j,k), ', rho    = ',rho(i,j,k), ', eps    = ',eps(i,j,k),&
           !                  '             x    = ',x(i,j,k),', y      = ',y(i,j,k),', z      = ',z(i,j,k)
           !   call CCTK_INFO(infoline)
           !end if

           if (CCTK_EQUALS(Spritz_eos_type,"Polytype")) then
              call prim2conpolytype_mhd(Spritz_eos_handle,gxx(i,j,k),gxy(i,j,k),&
                   gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
                   det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
                   tau(i,j,k),rho(i,j,k),&
                   velx(i,j,k),vely(i,j,k),velz(i,j,k),&
                   eps(i,j,k),press(i,j,k),w_lorentz(i,j,k), &
                   Bx(i,j,k), By(i,j,k), Bz(i,j,k), &
                   dBx(i,j,k), dBy(i,j,k), dBz(i,j,k), alp(i,j,k),&
                   betax(i,j,k),betay(i,j,k),betaz(i,j,k))
           else
              call prim2con_mhd(Spritz_eos_handle,gxx(i,j,k),gxy(i,j,k),&
                   gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
                   det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
                   tau(i,j,k),rho(i,j,k),&
                   velx(i,j,k),vely(i,j,k),velz(i,j,k),&
                   eps(i,j,k),press(i,j,k),w_lorentz(i,j,k), &
                   Bx(i,j,k), By(i,j,k), Bz(i,j,k), &
                   dBx(i,j,k), dBy(i,j,k), dBz(i,j,k),alp(i,j,k),&
                   betax(i,j,k),betay(i,j,k),betaz(i,j,k), &
                   Spritz_coldtab_rhomin,Spritz_coldtab_press_lo, &
                   Spritz_coldtab_eps_lo,Spritz_coldtab_gamma_lo, &
                   Spritz_coldtab_gammath,Spritz_coldtab_thfac)
           end if

           !if ((i .eq. 200) .and. (j .eq. 200) .and. (k .eq. 30)) then
           !   write(infoline,'(a20,g26.15,a11,g26.15,a11,g26.15,a20,g26.15,a11,g26.15,a11,g26.15)') &
           !                  'After P2C  : pres = ',press(i,j,k), ', rho    = ',rho(i,j,k), ', eps    = ',eps(i,j,k),&
           !                  '             x    = ',x(i,j,k),', y      = ',y(i,j,k),', z      = ',z(i,j,k)
           !   call CCTK_INFO(infoline)
           !end if

           rho_p(i,j,k)      = rho(i,j,k)
           dens_p(i,j,k)     = dens(i,j,k)

           scon_p(i,j,k,1)   = sx(i,j,k)
           scon_p(i,j,k,2)   = sy(i,j,k)
           scon_p(i,j,k,3)   = sz(i,j,k)
           vel_p(i,j,k,1)    = vel(i,j,k,1)
           vel_p(i,j,k,2)    = vel(i,j,k,2)
           vel_p(i,j,k,3)    = vel(i,j,k,3)

           tau_p(i,j,k)      = tau(i,j,k)
           press_p(i,j,k)    = press(i,j,k)
           eps_p(i,j,k)      = eps(i,j,k)

           Bx_p(i,j,k)       = Bx(i,j,k)
           By_p(i,j,k)       = By(i,j,k)
           Bz_p(i,j,k)       = Bz(i,j,k)
           dBx_p(i,j,k)      = dBx(i,j,k)
           dBy_p(i,j,k)      = dBy(i,j,k)
           dBz_p(i,j,k)      = dBz(i,j,k)
           Avec_x_p(i,j,k)   = Avec_x(i,j,k)
           Avec_y_p(i,j,k)   = Avec_y(i,j,k)
           Avec_z_p(i,j,k)   = Avec_z(i,j,k)

           if (timelevels==3) then
              rho_p_p(i,j,k)      = rho(i,j,k)
              dens_p_p(i,j,k)     = dens(i,j,k)

              scon_p_p(i,j,k,1)   = sx(i,j,k)
              scon_p_p(i,j,k,2)   = sy(i,j,k)
              scon_p_p(i,j,k,3)   = sz(i,j,k)
              vel_p_p(i,j,k,1)    = vel(i,j,k,1)
              vel_p_p(i,j,k,2)    = vel(i,j,k,2)
              vel_p_p(i,j,k,3)    = vel(i,j,k,3)

              tau_p_p(i,j,k)      = tau(i,j,k)
              press_p_p(i,j,k)    = press(i,j,k)
              eps_p_p(i,j,k)      = eps(i,j,k)

              Bx_p_p(i,j,k)       = Bx(i,j,k)
              By_p_p(i,j,k)       = By(i,j,k)
              Bz_p_p(i,j,k)       = Bz(i,j,k)
              dBx_p_p(i,j,k)      = dBx(i,j,k)
              dBy_p_p(i,j,k)      = dBy(i,j,k)
              dBz_p_p(i,j,k)      = dBz(i,j,k)
              Avec_x_p_p(i,j,k)   = Avec_x(i,j,k)
              Avec_y_p_p(i,j,k)   = Avec_y(i,j,k)
              Avec_z_p_p(i,j,k)   = Avec_z(i,j,k)
           end if

           !write(infoline,'(a25,g20.7,a8,g20.7,a8,g20.7)') '     input values:  Bx = ',Bx(i,j,k), &
           !                                                ',  By = ',By(i,j,k),',  Bz = ',Bz(i,j,k)
           !call CCTK_INFO(infoline)


        end do
     end do
  end do

  densrhs = 0.d0
  srhs = 0.d0
  taurhs = 0.d0
  Bxrhs = 0.0d0
  Byrhs = 0.0d0
  Bzrhs = 0.0d0
  Avec_x_rhs = 0.0d0
  Avec_y_rhs = 0.0d0
  Avec_z_rhs = 0.0d0

  Ex=0.0d0
  Ey=0.0d0
  Ez=0.0d0

  return

end subroutine Spritz_3DTests_MHD

function sphexp_prof(rr,rin,rout,fin,fout)
  ! Exponential decay in given region between given values
  implicit none

  CCTK_REAL :: sphexp_prof
  CCTK_REAL :: rr,rin,rout,fin,fout

  if(rr .lt. rin) then 
     sphexp_prof = fin
  else if( rr .gt. rout) then 
     sphexp_prof = fout
  else 
     sphexp_prof = exp( ( (log(fin)*(rout-rr) + log(fout)*(rr-rin)) / (rout - rin) ) ) 
  end if

  return 

end function sphexp_prof
