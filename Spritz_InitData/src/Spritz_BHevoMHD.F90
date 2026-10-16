 /*@@
   @file      Spritz_BHevoMHD.F90
   @date      Nov 2019
   @author    Federico Cipolletta
   @desc 
   Initial data for MHD evolution of BHs.
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

subroutine Spritz_BHevoMHD(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k, nx, ny, nz
  CCTK_REAL :: dx, dy, dz
  CCTK_REAL :: det 
  CCTK_REAL :: rhoin, pressin, epsin, &
               Bxin, Byin, Bzin, & 
               velxin, velyin, velzin

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

  write(infoline,'(a52,a57)') 'Setting ID for BH test. Add Bx2D_init, By2D_init and', &
                              ' Bz2D_init to parfile if you want non-zero magnetic field'
  call CCTK_INFO(infoline)
  
  ! Assign initial parameter, in order to compute initial
  ! data on gridpoints
  if (CCTK_EQUALS(Spritz_initial_data,"BBH")) then
     ! Setting magneto-hydro values for BH 

     ! Values for input     
     rhoin = rho_abs_min !Set everything to atmosphere
     !Since I don't know what to do in this case for eps and press
     !I set them using a polytropic EOS which means that the atmosphere is always cold
     xeps = 0.0d0 !This is a dummy value, not used in case of polytropic EOS
     call EOS_Omni_press(Spritz_polytrope_handle,keytemp,Eos_rf_prec,n,&
                         xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
     pressin = xpress(1)
     !Following doc for EOS_Omni, I reset the keytemp flag in order to eval
     !eps from the temperature
     keytemp = 1
     call EOS_Omni_EpsFromPress(Spritz_polytrope_handle,keytemp,Eos_rf_prec,n,&
                                xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)         
     epsin = xeps(1)
     velxin = 0.0d0
     velyin = 0.0d0
     velzin = 0.0d0
     Bxin = Bx2D_init
     Byin = By2D_init
     Bzin = Bz2D_init

     write(infoline,'(a7,g20.7,a9,g20.7,a9,g20.7)') 'Bxin = ',Bxin,', Byin = ',Byin,', Bzin = ',Bzin
     call CCTK_INFO(infoline)
     
  else
     call CCTK_WARN(0,"Should not be here...")
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
  
           !Assuming metric already assigned via other thorns
           !e.g. TwoPunctures

           call SpatialDet(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),&
                           gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det)

           ! Computing initial data on gridpoints
           if (CCTK_EQUALS(BHevo_case,"Single Uniform")) then
              
              rho(i,j,k) = rhoin                             
              press(i,j,k) = pressin
              eps(i,j,k) = epsin
              velx(i,j,k) = velxin
              vely(i,j,k) = velyin
              velz(i,j,k) = velzin

              Bx(i,j,k) = Bxin/sqrt(det)
              By(i,j,k) = Byin/sqrt(det)
              Bz(i,j,k) = Bzin/sqrt(det)
              dBx(i,j,k) = Bxin
              dBy(i,j,k) = Byin
              dBz(i,j,k) = Bzin

              condition = (abs(Bxin) > 1.0d-10) .or. &
                          (abs(Byin) > 1.0d-10) .or. &
                          (abs(Bzin) > 1.0d-10)

              !write(infoline,'(a12,l1)') 'condition = ',condition
              !call CCTK_INFO(infoline)

              if (condition) then
                 ! The following components of the magnetic vector potential, should fulfill the definition of B = curl(Avec)
                 ! for a staggered vector potential.           
                 Avec_x(i,j,k) = ( ( dBy(i,j,k)*(z(i,j,k) + (dz/2.0d0)) ) )
                 Avec_y(i,j,k) = ( ( dBz(i,j,k)*(x(i,j,k) + (dx/2.0d0)) ) )
                 Avec_z(i,j,k) = ( ( dBx(i,j,k)*(y(i,j,k) + (dy/2.0d0)) ) )
              else
                 Avec_x(i,j,k) = 0.0d0
                 Avec_y(i,j,k) = 0.0d0
                 Avec_z(i,j,k) = 0.0d0
              end if

           else
              call CCTK_WARN(0,"BH case not recognized")
           end if
           ! End of computation of initial data on gridpoints
        
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

end subroutine Spritz_BHevoMHD
