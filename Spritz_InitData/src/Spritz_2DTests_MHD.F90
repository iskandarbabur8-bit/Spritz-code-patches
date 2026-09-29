 /*@@
   @file      Spritz_2DTests_MHD.F90
   @date      Dec 2018
   @author    Federico Cipolletta
   @desc
   Initial data for MHD tests in 2D. Based on general structure of Spritz_ShockTube.F90,
   by Bruno Giacomazzo.
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

subroutine Spritz_2DTests_MHD(CCTK_ARGUMENTS)

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
               vmaxperp, &
               Omegain, Omegaout, &
               radiusin, radiusout

  CCTK_REAL :: cylblw_prof
  CCTK_REAL :: Aloop, rtemp

  CCTK_REAL :: hval, valf, t1, t2, t3, & !AA Defined as parameter
               pi, range_x, range_y, range_z, &
               wnbr, diaglength, range_d, &
               cos_theta, sin_theta, xnew, &
               vparallel, vperp, &
               Bparallel, Bperp

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
  if (CCTK_EQUALS(twodtest_case,"CylBlastExp")) then
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
     Omegain = OM2D_in !0.0d0

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
     Omegaout = OM2D_out !0.0d0

     !write(infoline,'(a20,g26.15,a11,g26.15,a11,g26.15,a20,g26.15,a11,g26.15,a11,g26.15)') &
     !               'CylExp CASE: pin  = ',pressin, ', rhoin  = ',rhoin, ', epsin  = ',epsin,&
     !               '             pout = ',pressout,', rhoout = ',rhoout,', epsout = ',epsout
     !call CCTK_INFO(infoline)

  else if (CCTK_EQUALS(twodtest_case,"MagRot")) then
     ! See GRHydro paper and Del Zanna, Bucciantini, Londrillo (2003)

     ! Values IN
     radiusin = R2D_in !0.1d0
     rhoin = rho2D_in !10.0d0
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
     Omegain = OM2D_in !9.95 ! From literature - In principle one could give any rotation profile

     !Values OUT
     radiusout = R2D_out !0.0d0 ! *** NOT USED ***
     rhoout = rho2D_out !1.0d0
     pressout = 1.0d0
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
     Omegaout = OM2D_out !0.0d0 ! Outside the rotor the fluid is at rest

  else if (CCTK_EQUALS(twodtest_case,"Alfven")) then
     ! Adapted from GRHydro public code.
     ! See GRHydro paper,
     !     Beckwith, Stone (2011),
     !     Del Zanna et al (2007).
     ! NOTE: a direction of the shock should be defined in
     !       order to place the wave at an oblique angle to
     !       the grid. 

     !!$pressure, density, B^x, specific internal energy and enthalpy:
     rhoin   = rho2D_in    !Should be set to 1.0
     pressin = P2D_in      !Should be set to 1.0
     Bxin    = Bx2D_init   !Should be set to 1.0
     
     !xrho = rhoin
     !xpress = pressin
     !call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
     !                           xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     !epsin = xeps(1)

     epsin = pressin/((5.0d0/3.0d0)-1.0d0)/rhoin
     hval   = 1.0d0 + epsin + pressin/rhoin

     !!$ DZ: rho=P=Bxval=AA = 1
     !!$ Using DZ parameters, epsval=1.5, hval=3.5
  
     !!$ Alfven Wave Amplitude: (CAN BE PASSED AS PARAMETER) !Default is 1.0
     !AA=1.0d0

     !!$ Alfven wave speed:
     t1 = rhoin*hval+Bxin**2*(1.0d0+AA**2)
     t2 = 2.0d0*AA*Bxin**2/t1
     t3 = 0.5d0*(1.0d0+sqrt(1.0d0-t2**2))
     valf = sqrt(Bxin**2/t1/t3)

     write(infoline,'(a17,g20.7)')'Alfven velocity: ',valf
     call CCTK_INFO(infoline)
     write(infoline,'(a17,g20.7)')'rhoin = ',rhoin
     call CCTK_INFO(infoline)
     write(infoline,'(a17,g20.7)')'pressin = ',pressin
     call CCTK_INFO(infoline)
     write(infoline,'(a17,g20.7)')'epsin = ',epsin
     call CCTK_INFO(infoline)
     write(infoline,'(a17,g20.7)')'Bxin = ',Bxin
     call CCTK_INFO(infoline)
     write(infoline,'(a17,g20.7)')'hval = ',hval
     call CCTK_INFO(infoline)
     !!$ Using DZ parameters with P=1.0, we have:
     !!$ t1=5.5
     !!$ t2=4/11
     !!$ t3=0.96577
     !!$ valf=0.43389

     !!$ Vx value:
     velxin=velx_init !Default is 0

  else if (CCTK_EQUALS(twodtest_case,"LoopAdv")) then
     ! See GRHydro paper and Beckwith, Stone (2011)

     Aloop = Aloop_init

     ! Values IN
     radiusin =  R2D_in   ! 0.3d0	! R2D_in
     rhoin = 	 rho2D_in ! 1.0d0	! rho2D_in
     pressin = 	 P2D_in	  ! 3.0d0	! P2D_in
     xrho = rhoin
     xpress = pressin
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsin = xeps(1)
     Bxin = 0.0d0 ! Bx2D_init
     Byin = 0.0d0 ! By2D_init
     Bzin = 0.0d0 ! By2D_init
     velxin = velx_init         !		1.0d0/12.0d0
     velyin = vely_init		!		1.0d0/24.0d0
     velzin = velz_init		! 		0.0d0

     !Values OUT
     radiusout = 0.0d0 !R2D_out ! *** NOT USED ***
     rhoout =    rho2D_out
     pressout =  P2D_out
     xrho = rhoout
     xpress = pressout
     call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     epsout = xeps(1)
     Bxout = 0.0d0
     Byout = 0.0d0
     Bzout = 0.0d0
     velxout = 0.0d0
     velyout = 0.0d0
     velzout = 0.0d0

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

  ! The following 4 lines are only used in Alven Wave Test
  ! Taken from GRHydro code
  !!$ Note that the 3D test wasn't deviced to be used with AMR!
  pi=4.0d0*atan(1.0d0)
  range_x = 1.0d0 !(cctk_gsh(1)-2*cctk_nghostzones(1))*dx
  range_y = (cctk_gsh(2)-2*cctk_nghostzones(2))*dy
  range_z = (cctk_gsh(3)-2*cctk_nghostzones(3))*dz

  !range_d = sqrt(range_x**2+range_y**2)
  !diaglength=range_x*range_y/range_d
  !cos_theta = range_y/range_d
  !sin_theta = range_x/range_d

  write(infoline,'(a17,g20.7)')'pi = ',pi
  call CCTK_INFO(infoline)
  write(infoline,'(a17,g20.7)')'range_x = ',range_x
  call CCTK_INFO(infoline)
  write(infoline,'(a17,g20.7)')'range_y = ',range_y
  call CCTK_INFO(infoline)
  write(infoline,'(a17,g20.7)')'range_z = ',range_z
  call CCTK_INFO(infoline)
  !write(infoline,'(a17,g20.7)')'range_d = ',range_d
  !call CCTK_INFO(infoline)
  !write(infoline,'(a17,g20.7)')'diaglength = ',diaglength
  !call CCTK_INFO(infoline)
  !write(infoline,'(a17,g20.7)')'cos_theta = ',cos_theta
  !call CCTK_INFO(infoline)
  !write(infoline,'(a17,g20.7)')'sin_theta = ',sin_theta
  !call CCTK_INFO(infoline)

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
           if (CCTK_EQUALS(twodtest_case,"CylBlastExp")) then
              ! See GRHydro paper, Beckwith, Stone (2011) and Komissarov (1999)
              ! Assuming the center of the cylinder to be (xc,yc,zc) given as input
              ! and its axis directed on cylexp_axis_dir

              if (CCTK_EQUALS(cylexp_axis_dir,"x_dir")) then
                 rtemp = sqrt( ( (y(i,j,k) - yc) * (y(i,j,k) - yc) + &
                                 (z(i,j,k) - zc) * (z(i,j,k) - zc) ) )
              else if (CCTK_EQUALS(cylexp_axis_dir,"y_dir")) then
                 rtemp = sqrt( ( (z(i,j,k) - zc) * (z(i,j,k) - zc) + &
                                 (x(i,j,k) - xc) * (x(i,j,k) - xc) ) )
              else if (CCTK_EQUALS(cylexp_axis_dir,"z_dir")) then
                 rtemp = sqrt( ( (x(i,j,k) - xc) * (x(i,j,k) - xc) + &
                                 (y(i,j,k) - yc) * (y(i,j,k) - yc) ) )
              else
                 call CCTK_WARN(0,"Cylindrical Explosion case not recognized")
              end if

              rho(i,j,k) = cylblw_prof(rtemp,radiusin,radiusout,rhoin,rhoout)
              press(i,j,k) = cylblw_prof(rtemp,radiusin,radiusout,pressin,pressout)           !!!!!
              !eps(i,j,k) = cylblw_prof(rtemp,radiusin,radiusout,epsin,epsout)
              !eps(i,j,k) = cylblw_prof(rtemp,radiusin,radiusout,pressin,pressout)/((gl_gamma-1.0d0)*rho(i,j,k))
              n = 1
              keytemp = 0
              anyerr = 0
              keyerr = 0
              xrho = rho(i,j,k)
              xpress = press(i,j,k)
              !xeps = eps(i,j,k)
              xtemp = 0.0d0
              xye = 0.0d0
              !call EOS_Omni_press(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
	      !                    xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
              !press(i,j,k) = xpress(1)
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

              ! NOTE: there is an adiition to half of gridspacing because of staggering!!!

              if (CCTK_EQUALS(cylexp_axis_dir,"x_dir")) then
                 Avec_x(i,j,k) = ( By(i,j,k)*(z(i,j,k) + (dz/2.0d0)) ) !0.5d0 * ( (By(i,j,k)*z(i,j,k)) - (Bz(i,j,k)*y(i,j,k)) )
                 Avec_y(i,j,k) = 0.0d0 !0.5d0 * ( (Bz(i,j,k)*x(i,j,k)) - (Bx(i,j,k)*z(i,j,k)) )
                 Avec_z(i,j,k) = 0.0d0 !0.5d0 * ( (Bx(i,j,k)*y(i,j,k)) - (By(i,j,k)*x(i,j,k)) )
              else if (CCTK_EQUALS(cylexp_axis_dir,"y_dir")) then
                 Avec_x(i,j,k) = 0.0d0 !0.5d0 * ( (By(i,j,k)*z(i,j,k)) - (Bz(i,j,k)*y(i,j,k)) )
                 Avec_y(i,j,k) = ( Bz(i,j,k)*(x(i,j,k) + (dx/2.0d0)) ) !0.5d0 * ( (Bz(i,j,k)*x(i,j,k)) - (Bx(i,j,k)*z(i,j,k)) )
                 Avec_z(i,j,k) = 0.0d0 !0.5d0 * ( (Bx(i,j,k)*y(i,j,k)) - (By(i,j,k)*x(i,j,k)) )
              else if (CCTK_EQUALS(cylexp_axis_dir,"z_dir")) then
                 Avec_x(i,j,k) = 0.0d0 !0.5d0 * ( (By(i,j,k)*z(i,j,k)) - (Bz(i,j,k)*y(i,j,k)) )
                 Avec_y(i,j,k) = 0.0d0 !0.5d0 * ( (Bz(i,j,k)*x(i,j,k)) - (Bx(i,j,k)*z(i,j,k)) )
                 Avec_z(i,j,k) = ( Bx(i,j,k)*(y(i,j,k) + (dy/2.0d0)) ) !0.5d0 * ( (Bx(i,j,k)*y(i,j,k)) - (By(i,j,k)*x(i,j,k)) )
              else
                 call CCTK_WARN(0,"Cylindrical Explosion case not recognized")
              end if

              !if ((abs(x(i,j,k) - 2.0d0) .lt. 1.0d0)) then
              !   write(warnline,'(a20,g26.15,a11,g26.15,a11,g26.15,a20,g26.15,a11,g26.15,a11,g26.15)') &
              !         'CylExp CASE: p  = ',press(i,j,k), ', rho  = ',rho(i,j,k), ', eps  = ',eps(i,j,k),&
              !         '             Bx = ',Bx(i,j,k),', By = ',By(i,j,k),', Bz = ',Bz(i,j,k)
              !   call CCTK_WARN(1,warnline)
              !end if

           else if (CCTK_EQUALS(twodtest_case,"MagRot")) then
              ! See GRHydro paper and Del Zanna, Bucciantini, Londrillo (2003)
              ! Assuming the center of the cylinder to be (xc,yc,0) given as input

              rtemp = sqrt( ( (x(i,j,k) - xc) * (x(i,j,k) - xc) ) + &
                            ( (y(i,j,k) - yc) * (y(i,j,k) - yc) ) )
              vmaxperp = Omegain * radiusin
              condition = ( (rtemp .gt. radiusin) )
              !write(warnline,'(a8,g26.15,a5,g26.15,a9,g26.15,a7,l2)') &
              !               '!!! x = ',x(i,j,k),' y = ',y(i,j,k),' rtemp = ',rtemp,' cond =',condition
              !call CCTK_WARN(1,warnline)

              if( condition ) then
                rho(i,j,k) = rhoout
                press(i,j,k) = pressout
                velx(i,j,k) = 0.0d0
                vely(i,j,k) = 0.0d0
                velz(i,j,k) = 0.0d0
              else
              	rho(i,j,k) = rhoin
                press(i,j,k) = pressin
                velx(i,j,k) = ( -1.0d0 * vmaxperp*(y(i,j,k) - yc) )/radiusin
                vely(i,j,k) = ( 1.0d0 * vmaxperp*(x(i,j,k) - xc) )/radiusin
                velz(i,j,k) = 0.0d0
              end if

              keyerr = 0
              anyerr = 0
              call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                         rho(i:i,j:j,k:k),eps(i:i,j:j,k:k),xtemp,xye,&
                                         press(i:i,j:j,k:k),eps(i:i,j:j,k:k),&
                                         keyerr,anyerr)
              Bx(i,j,k) = Bx2D_init
              By(i,j,k) = By2D_init
              Bz(i,j,k) = Bz2D_init
              ! The following components of the magnetic vector potential, have been found by trials
              ! but they fulfill the definition of B = curl(Avec)
              Avec_x(i,j,k) = 0.0d0
              Avec_y(i,j,k) = 0.0d0
              Avec_z(i,j,k) = (Bx(i,j,k)*y(i,j,k))

           else if (CCTK_EQUALS(twodtest_case,"Alfven")) then
              ! Adapted from GRHydro public code
              ! See GRHydro paper,
              !     Beckwith, Stone (2011),
              !     Del Zanna et al (2007).
              ! NOTE: There should be used the direction defined
              !       in **** (or define it in here...)

              rho(i,j,k)=rhoin
              press(i,j,k)=pressin
              eps(i,j,k)=epsin
              !call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
              !                           rho(i:i,j:j,k:k),eps(i:i,j:j,k:k),xtemp,xye,&
              !                           press(i:i,j:j,k:k),eps(i:i,j:j,k:k),&
              !                           keyerr,anyerr)

              if (CCTK_EQUALS(alfvenwave_type,"1D")) then

                wnbr = 2.0d0*pi/range_x

                velx(i,j,k)=velxin
                vely(i,j,k)=-valf*AA*cos(wnbr*x(i,j,k))
                velz(i,j,k)=-valf*AA*sin(wnbr*x(i,j,k))
                Bx(i,j,k)=Bxin
                By(i,j,k)=AA*Bxin*cos(wnbr*x(i,j,k))
                Bz(i,j,k)=AA*Bxin*sin(wnbr*x(i,j,k))
                dBx(i,j,k)=Bx(i,j,k)
                dBy(i,j,k)=By(i,j,k)
                dBz(i,j,k)=Bz(i,j,k)

                !Assuming staggered vector potential
                !Found analytically
                
                !Avec_x(i,j,k) = - ( AA*Bxin*sin(wnbr*x(i,j,k)) )
                !Avec_y(i,j,k) = 0.0d0
                !Avec_z(i,j,k) = - (AA*Bxin*sin(wnbr*x(i,j,k)) / wnbr ) + ( Bxin*y(i,j,k) )

                Avec_x(i,j,k) = 0.0d0 !- ( AA*Bxin*sin(wnbr*x(i,j,k)) )
                Avec_y(i,j,k) = - (AA*Bxin*cos(wnbr*x(i,j,k)) / wnbr )
                Avec_z(i,j,k) = - (AA*Bxin*sin(wnbr*x(i,j,k)) / wnbr ) + ( Bxin*y(i,j,k) )

              else if (CCTK_EQUALS(alfvenwave_type,"2D")) then

                range_d = sqrt(range_x**2+range_y**2)
                diaglength=range_x*range_y/range_d
                cos_theta = range_y/range_d
                sin_theta = range_x/range_d
                wnbr = 2.0d0*pi/diaglength

                xnew = cos_theta*x(i,j,k)+sin_theta*y(i,j,k)

                vparallel=velxin
                vperp=-valf*AA*cos(wnbr*xnew)
                velx(i,j,k)=vparallel*cos_theta-vperp*sin_theta
                vely(i,j,k)=vparallel*sin_theta+vperp*cos_theta
                velz(i,j,k)=-valf*AA*sin(wnbr*xnew)
                Bparallel=Bxin
                Bperp=AA*Bxin*cos(wnbr*xnew)
                Bx(i,j,k)=Bparallel*cos_theta-Bperp*sin_theta
                By(i,j,k)=Bparallel*sin_theta+Bperp*cos_theta
                Bz(i,j,k)=AA*Bxin*sin(wnbr*xnew)
                dBx(i,j,k)=Bx(i,j,k)
                dBy(i,j,k)=By(i,j,k)
                dBz(i,j,k)=Bz(i,j,k)

                !Vector Potential to be defined
                Avec_x(i,j,k) = 0.0d0
                Avec_y(i,j,k) = 0.0d0
                Avec_z(i,j,k) = 0.0d0

              else
                call CCTK_WARN(0,"Alfven wave case not recognized!")
              end if

           else if (CCTK_EQUALS(twodtest_case,"LoopAdv")) then
              ! See GRHydro paper and Beckwith, Stone (2011)

              rtemp = sqrt( ( (x(i,j,k) - xc) * (x(i,j,k) - xc) ) + &
                            ( (y(i,j,k) - yc) * (y(i,j,k) - yc) ) )

              velx(i,j,k) = velxin
              vely(i,j,k) = velyin
              velz(i,j,k) = velzin

              condition = ( (rtemp .gt. radiusin) )

              if( condition ) then
                rho(i,j,k) = rhoout
                press(i,j,k) = pressout
                Avec_x(i,j,k) = 0.0d0
                Avec_y(i,j,k) = 0.0d0
                Avec_z(i,j,k) = 0.0d0
                Bx(i,j,k) = 0.0d0
                By(i,j,k) = 0.0d0
                Bz(i,j,k) = 0.0d0
                Y_e(i,j,k)= 0.0d0
              else
              	rho(i,j,k) = rhoin
                press(i,j,k) = pressin
                Avec_x(i,j,k) = 0.0d0
                Avec_y(i,j,k) = 0.0d0
                Avec_z(i,j,k) = Aloop*(radiusin-rtemp)
                if(x(i,j,k).eq.0.0d0) then
                  Bx(i,j,k) = 0.0d0
                else
                  Bx(i,j,k) = -1.0d0*Aloop*y(i,j,k)/rtemp    !Bx2D_init
                end if
                if(y(i,j,k).eq.0.0d0) then
                  By(i,j,k) = 0.0d0
                else
                  By(i,j,k) = Aloop*x(i,j,k)/rtemp           !By2D_init
                end if
                Bz(i,j,k) = 0.0d0			   !Bz2D_init
                Y_e(i,j,k) = Y_e_in_loop*(radiusin-rtemp)/radiusin !default is 0 everywhere
              end if

              keyerr = 0
              anyerr = 0
              call EOS_Omni_EpsFromPress(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                                         rho(i:i,j:j,k:k),eps(i:i,j:j,k:k),xtemp,xye,&
                                         press(i:i,j:j,k:k),eps(i:i,j:j,k:k),&
                                         keyerr,anyerr)

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
                   betax(i,j,k),betay(i,j,k),betaz(i,j,k))
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
           Y_e_p(i,j,k)      = Y_e(i,j,k)

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
              Y_e_p_p(i,j,k)      = Y_e(i,j,k)

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

end subroutine Spritz_2DTests_MHD

function cylblw_prof(rr,rin,rout,fin,fout)
  ! From GRHydro paper - Note a typo in eq (96) of GRHydro paper wrt to its code
  ! This corresponds to the expression used in the GRHydro code
  implicit none

  CCTK_REAL :: cylblw_prof
  CCTK_REAL :: rr,rin,rout,fin,fout

  if(rr .lt. rin) then
     cylblw_prof = fin
  else if( rr .gt. rout) then
     cylblw_prof = fout
  else
     cylblw_prof = exp( ( (log(fin)*(rout-rr) + log(fout)*(rr-rin)) / (rout - rin) ) )
  end if

  return

end function cylblw_prof
