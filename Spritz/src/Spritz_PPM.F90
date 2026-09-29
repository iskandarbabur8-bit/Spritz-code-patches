 /*@@
   @file      Spritz_PPM.F90
   @date      Sun Feb 10 16:53:29 2002
   @author    Ian Hawke, Toni Font, Luca Baiotti, Frank Loeffler
   @desc
   Routines to do PPM reconstruction.
   @enddesc
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine PPM_TVD(origm, orig, origp, bextm, bextp)
  CCTK_REAL :: origm, orig, origp, bextm, bextp
  CCTK_REAL :: dloc, dupw, delta

  dupw = orig - origm
  dloc = origp - orig
  if (dupw*dloc < 0.d0) then
    delta=0.d0
  else if (abs(dupw) < abs(dloc)) then
    delta=dupw
  else
    delta=dloc
  end if
  bextm = orig - 0.5d0 * delta
  bextp = orig + 0.5d0 * delta
end subroutine PPM_TVD



 /*@@
   @routine    SimplePPM_1d_MHD
   @date
   @author     Bruno Giacomazzo
   @desc
   It is simply the PPM scheme extended in a naive way
   to MHD. No information on the hyperbolic structure of the equations
   is used.
   @enddesc
   @calls
   @calledby
   @history
   @endhistory

@@*/

  subroutine SimplePPM_1d_MHD(handle,poly,&
       & nx,dx,rho,velx,vely,velz,eps,press,Bx,By,Bz,rhominus,&
       & velxminus,velyminus,velzminus,epsminus,Bxminus,Byminus,Bzminus,&
       & rhoplus,velxplus,velyplus,&
       & velzplus,epsplus,Bxplus,Byplus,Bzplus,trivial_rp,&
       & gxx, gxy, gxz, gyy, gyz, gzz, psi4, beta, alp, w_lorentz, &
       & dir, ni, nj, nrx, nry, nrz, Spritz_Cs2_min, &
       & psi_mhd, psi_mhd_minus, psi_mhd_plus )

    USE Spritz_Scalars

    implicit none

    DECLARE_CCTK_PARAMETERS

!#ifdef _EOS_BASE_INC_
!#undef _EOS_BASE_INC_
!#endif
!#include "EOS_Base.inc"

    CCTK_INT :: handle,poly,nx
    CCTK_REAL :: dx
    CCTK_REAL, dimension(nx) :: rho,velx,vely,velz,eps,Bx,By,Bz
    CCTK_REAL, dimension(nx) :: rhominus,velxminus,velyminus,velzminus,epsminus,Bxminus,Byminus,Bzminus
    CCTK_REAL, dimension(nx) :: rhoplus,velxplus,velyplus,velzplus,epsplus,Bxplus,Byplus,Bzplus
    CCTK_REAL, dimension(nx) :: press
    CCTK_REAL :: Spritz_Cs2_min

    CCTK_REAL, dimension(nx),OPTIONAL :: psi_mhd, psi_mhd_minus, psi_mhd_plus

    CCTK_INT :: i,s
    CCTK_REAL, dimension(nx) :: drho,dvelx,dvely,dvelz,deps,dBx,dBy,dBz, dpsi_mhd
    CCTK_REAL, dimension(nx) :: dmrho,dmvelx,dmvely,dmvelz,dmeps,dmBx,dmBy,dmBz, dmpsi_mhd
    CCTK_REAL, dimension(nx) :: dpress,d2rho,tilde_flatten
    CCTK_REAL :: dpress2,dvel,w,flatten,eta,etatilde

    logical, dimension( nx) :: trivial_rp



    CCTK_REAL, dimension( nx) :: gxx, gxy, gxz, gyy, gyz, gzz, &
         psi4, beta, alp, w_lorentz
    CCTK_INT :: dir, ni, nj, nrx, nry, nrz
    CCTK_REAL, dimension(nrx, nry, nrz) :: ev_l, ev_r, xw

    CCTK_REAL :: uxx, uxy, uxz, uyy, uyz, uzz, det
    CCTK_REAL, dimension(5) :: lam
    CCTK_REAL :: agxx, agxy, agxz, agyy, agyz, agzz
    CCTK_REAL, dimension(nx) :: xwind, l_ev_l, l_ev_r
    CCTK_REAL, parameter :: one_sixth = 1.d0/6.d0, one_twelfth = 1.d0/12.d0

    logical :: cond


!!$  Initially, all the Riemann problems will be trivial

    trivial_rp(:) = .true.

!!$  Average slopes delta_m(a). See (1.7) of Colella and Woodward, p.178
!!$  This is the expression for an even grid.

    do i = 2, nx - 1
       drho(i) = 0.5d0 * (rho(i+1) - rho(i-1))
       dvelx(i) = 0.5d0 * (velx(i+1) - velx(i-1))
       dvely(i) = 0.5d0 * (vely(i+1) - vely(i-1))
       dvelz(i) = 0.5d0 * (velz(i+1) - velz(i-1))
       dBx(i) = 0.5d0 * (Bx(i+1) - Bx(i-1))
       dBy(i) = 0.5d0 * (By(i+1) - By(i-1))
       dBz(i) = 0.5d0 * (Bz(i+1) - Bz(i-1))
       if (present(psi_mhd)) then
          dpsi_mhd(i) = 0.5d0 * (psi_mhd(i+1) - psi_mhd(i-1))
       end if
       dpress(i) = press(i+1) - press(i-1)
       d2rho(i) = (rho(i+1) - 2.d0 * rho(i) + rho(i-1))! / 6.d0 / dx / dx
       ! since we use d2rho only for the condition d2rho(i+1)*d2rhoi(i-1)<0
       ! the denominator is not necessary
    end do
    if (poly .eq. 0) then
       do i = 2, nx - 1
          deps(i) = 0.5d0 * (eps(i+1) - eps(i-1))
       end do
    end if

!!$  Steepened slope. See (1.8) of Colella and Woodward, p.178

    do i = 2, nx - 1
#define STEEP(x,dx,dmx)                                              \
       if ( (x(i+1) - x(i)) * (x(i) - x(i-1)) > 0.d0 ) then           &&\
          dmx(i) = sign(1.d0, dx(i)) *                                   \
          min(abs(dx(i)), 2.d0 * abs(x(i) - x(i-1)),                \
        2.d0 * abs(x(i+1) - x(i)))              &&\
       else                                                           &&\
          dmx(i) = 0.d0                                                &&\
       end if
       STEEP(rho, drho, dmrho)
       STEEP(velx, dvelx, dmvelx)
       STEEP(vely, dvely, dmvely)
       STEEP(velz, dvelz, dmvelz)
       STEEP(Bx, dBx, dmBx)
       STEEP(By, dBy, dmBy)
       STEEP(Bz, dBz, dmBz)
       if (present(psi_mhd)) then
          STEEP(psi_mhd, dpsi_mhd, dmpsi_mhd)
       end if
    end do
    if (poly .eq. 0) then
       do i = 2, nx - 1
          STEEP(eps, deps, dmeps)
       end do
    end if

!!$  Initial boundary states. See (1.9) of Colella and Woodward, p.178

    do i = 2, nx-2
       rhoplus(i) = 0.5d0 * (rho(i) + rho(i+1)) + &
            (dmrho(i) - dmrho(i+1)) * one_sixth
       rhominus(i+1) = rhoplus(i)
       velxplus(i) = 0.5d0 * (velx(i) + velx(i+1)) + &
            (dmvelx(i) - dmvelx(i+1)) * one_sixth
       velxminus(i+1) = velxplus(i)
       velyplus(i) = 0.5d0 * (vely(i) + vely(i+1)) + &
            (dmvely(i) - dmvely(i+1)) * one_sixth
       velyminus(i+1) = velyplus(i)
       velzplus(i) = 0.5d0 * (velz(i) + velz(i+1)) + &
            (dmvelz(i) - dmvelz(i+1)) * one_sixth
       velzminus(i+1) = velzplus(i)
       Bxplus(i) = 0.5d0 * (Bx(i) + Bx(i+1)) + &
            (dmBx(i) - dmBx(i+1)) * one_sixth
       Bxminus(i+1) = Bxplus(i)
       Byplus(i) = 0.5d0 * (By(i) + By(i+1)) + &
            (dmBy(i) - dmBy(i+1)) * one_sixth
       Byminus(i+1) = Byplus(i)
       Bzplus(i) = 0.5d0 * (Bz(i) + Bz(i+1)) + &
            (dmBz(i) - dmBz(i+1)) * one_sixth
       Bzminus(i+1) = Bzplus(i)
       if (present(psi_mhd)) then
          psi_mhd_plus(i) = 0.5d0 * (psi_mhd(i) + psi_mhd(i+1)) + &
               (dmpsi_mhd(i) - dmpsi_mhd(i+1)) * one_sixth
          psi_mhd_minus(i+1) = psi_mhd_plus(i)
       end if
    end do
    if (poly .eq. 0) then
       do i = 2, nx-2
          epsplus(i) = 0.5d0 * (eps(i) + eps(i+1)) + &
               (dmeps(i) - dmeps(i+1)) * one_sixth
          epsminus(i+1) = epsplus(i)
       end do
    end if

!!$Discontinuity steepening. See (1.14-17) of C&W.
!!$This is the detect routine which mat be activated with the ppm_detect parameter
!!$Note that this part really also depends on the grid being even.
!!$Note also that we don''t have access to the gas constant gamma.
!!$So this is just dropped from eq. (3.2) of C&W.
!!$We can get around this by just rescaling the constant k0 (ppm_k0 here).

    if (ppm_detect .ne. 0) then

       do i = 3, nx - 2
          if ( (d2rho(i+1)*d2rho(i-1) < 0.d0).and.(abs(rho(i+1)-rho(i-1)) - &
               ppm_epsilon_shock * min(abs(rho(i+1)), abs(rho(i-1))) > 0.d0) ) then
             etatilde = (rho(i-2) - rho(i+2) + 4.d0 * drho(i)) / (drho(i) * 12.d0)
          else
             etatilde = 0.d0
          end if
          eta = max(0.d0, min(1.d0, ppm_eta1 * (etatilde - ppm_eta2)))
          if (ppm_k0 * abs(drho(i)) * min(press(i-1),press(i+1)) < &
               abs(dpress(i)) * min(rho(i-1), rho(i+1))) then
             eta = 0.d0
          end if
          if (eta > 0.d0) then
             trivial_rp(i-1) = .false.
             trivial_rp(i) = .false.
          end if
          rhominus(i) = rhominus(i) * (1.d0 - eta) + &
               (rho(i-1) + 0.5d0 * dmrho(i-1)) * eta
          rhoplus(i) = rhoplus(i) * (1.d0 - eta) + &
               (rho(i+1) - 0.5d0 * dmrho(i+1)) * eta
       end do

    end if


!!$  Zone flattening. See appendix of C&W, p. 197-8.

    do i = 3, nx - 2
       dpress2 = press(i+2) - press(i-2)
       dvel = velx(i+1) - velx(i-1)
       if ( (abs(dpress(i)) >  ppm_epsilon * min(press(i-1),press(i+1))) .and. &
            (dvel < 0.d0) ) then
          w = 1.d0
       else
          w = 0.d0
       end if
       if (abs(dpress2) < ppm_small) then
          tilde_flatten(i) = 1.d0
       else
          tilde_flatten(i) = max(0.d0, 1.d0 - w * max(0.d0, ppm_omega2 * &
               (dpress(i) / dpress2 - ppm_omega1)))
       end if
    end do



    if (PPM3) then !!$ Implement C&W, page 197, but with a workaround which allows to use stencil=3.
       do i = 3, nx - 2
          flatten = tilde_flatten(i)
          if (abs(1.d0 - flatten) > 0.d0) then
             trivial_rp(i-1) = .false.
             trivial_rp(i) = .false.
          end if
          rhoplus(i) = flatten * rhoplus(i) + (1.d0 - flatten) * rho(i)
          rhominus(i) = flatten * rhominus(i) + (1.d0 - flatten) * rho(i)
          velxplus(i) = flatten * velxplus(i) + (1.d0 - flatten) * velx(i)
          velxminus(i) = flatten * velxminus(i) + (1.d0 - flatten) * velx(i)
          velyplus(i) = flatten * velyplus(i) + (1.d0 - flatten) * vely(i)
          velyminus(i) = flatten * velyminus(i) + (1.d0 - flatten) * vely(i)
          velzplus(i) = flatten * velzplus(i) + (1.d0 - flatten) * velz(i)
          velzminus(i) = flatten * velzminus(i) + (1.d0 - flatten) * velz(i)
          Bxplus(i) = flatten * Bxplus(i) + (1.d0 - flatten) * Bx(i)
          Bxminus(i) = flatten * Bxminus(i) + (1.d0 - flatten) * Bx(i)
          Byplus(i) = flatten * Byplus(i) + (1.d0 - flatten) * By(i)
          Byminus(i) = flatten * Byminus(i) + (1.d0 - flatten) * By(i)
          Bzplus(i) = flatten * Bzplus(i) + (1.d0 - flatten) * Bz(i)
          Bzminus(i) = flatten * Bzminus(i) + (1.d0 - flatten) * Bz(i)
          if (present(psi_mhd)) then
             psi_mhd_plus(i) = flatten * psi_mhd_plus(i) + (1.d0 - flatten) * psi_mhd(i)
             psi_mhd_minus(i) = flatten * psi_mhd_minus(i) + (1.d0 - flatten) * psi_mhd(i)
          end if
          if (poly .eq. 0) then
             epsplus(i) = flatten * epsplus(i) + (1.d0 - flatten) * eps(i)
             epsminus(i) = flatten * epsminus(i) + (1.d0 - flatten) * eps(i)
          end if
       end do
    else  !!$ Really implement C&W, page 197; which requires stencil 4.
       do i = 4, nx - 3
          s=sign(1.d0, -dpress(i))
          flatten = min(tilde_flatten(i), tilde_flatten(i+s))  ! it must be min, because of our definition of flatten.
          ! in other words, flattening is applied to a point if the point or any adjacent point is found to be in a shock.
          if (abs(1.d0 - flatten) > 0.d0) then
             trivial_rp(i-1) = .false.
             trivial_rp(i) = .false.
          end if
          rhoplus(i) = flatten * rhoplus(i) + (1.d0 - flatten) * rho(i)
          rhominus(i) = flatten * rhominus(i) + (1.d0 - flatten) * rho(i)
          velxplus(i) = flatten * velxplus(i) + (1.d0 - flatten) * velx(i)
          velxminus(i) = flatten * velxminus(i) + (1.d0 - flatten) * velx(i)
          velyplus(i) = flatten * velyplus(i) + (1.d0 - flatten) * vely(i)
          velyminus(i) = flatten * velyminus(i) + (1.d0 - flatten) * vely(i)
          velzplus(i) = flatten * velzplus(i) + (1.d0 - flatten) * velz(i)
          velzminus(i) = flatten * velzminus(i) + (1.d0 - flatten) * velz(i)
          Bxplus(i) = flatten * Bxplus(i) + (1.d0 - flatten) * Bx(i)
          Bxminus(i) = flatten * Bxminus(i) + (1.d0 - flatten) * Bx(i)
          Byplus(i) = flatten * Byplus(i) + (1.d0 - flatten) * By(i)
          Byminus(i) = flatten * Byminus(i) + (1.d0 - flatten) * By(i)
          Bzplus(i) = flatten * Bzplus(i) + (1.d0 - flatten) * Bz(i)
          Bzminus(i) = flatten * Bzminus(i) + (1.d0 - flatten) * Bz(i)
          if (present(psi_mhd)) then
             psi_mhd_plus(i) = flatten * psi_mhd_plus(i) + (1.d0 - flatten) * psi_mhd(i)
             psi_mhd_minus(i) = flatten * psi_mhd_minus(i) + (1.d0 - flatten) * psi_mhd(i)
          end if
          if (poly .eq. 0) then
             epsplus(i) = flatten * epsplus(i) + (1.d0 - flatten) * eps(i)
             epsminus(i) = flatten * epsminus(i) + (1.d0 - flatten) * eps(i)
          end if
       end do
    end if


!!$ Monotonicity. See (1.10) of C&W.

    do i = Spritz_stencil, nx - Spritz_stencil + 1
#define MONOTONICITY(xminus,x,xplus)                                       \
       if (.not.( (xplus(i).eq.x(i)) .and. (x(i).eq.xminus(i)) )     \
          .and. ((xplus(i)-x(i))*(x(i)-xminus(i)) .le. 0.d0)) then&&\
          trivial_rp(i-1) = .false.                                 &&\
          trivial_rp(i) = .false.                                   &&\
          xminus(i) = x(i)                                          &&\
          xplus(i) = x(i)                                           &&\
       else if (6.d0 * (xplus(i) - xminus(i)) * (x(i) - 0.5d0 *      \
          (xplus(i) + xminus(i))) >                         \
          (xplus(i) - xminus(i))**2) then                 &&\
          xminus(i) = 3.d0 * x(i) - 2.d0 * xplus(i)                 &&\
          trivial_rp(i-1) = .false.                                 &&\
          trivial_rp(i) = .false.                                   &&\
       else if (6.d0 * (xplus(i) - xminus(i)) * (x(i) - 0.5d0 *      \
          (xplus(i) + xminus(i))) <                         \
          -(xplus(i) - xminus(i))**2) then                 &&\
          xplus(i) = 3.d0 * x(i) - 2.d0 * xminus(i)                 &&\
          trivial_rp(i-1) = .false.                                 &&\
          trivial_rp(i) = .false.                                   &&\
       end if&&\
       if (.not.( (xplus(i).eq.x(i)) .and. (x(i).eq.xminus(i)) ) ) then     &&\
          trivial_rp(i-1) = .false.                                 &&\
          trivial_rp(i) = .false.                                   &&\
       end if

       MONOTONICITY(rhominus,rho,rhoplus)
       MONOTONICITY(velxminus,velx,velxplus)
       MONOTONICITY(velyminus,vely,velyplus)
       MONOTONICITY(velzminus,velz,velzplus)
       MONOTONICITY(Bxminus,Bx,Bxplus)
       MONOTONICITY(Byminus,By,Byplus)
       MONOTONICITY(Bzminus,Bz,Bzplus)
       if (present(psi_mhd)) then
          MONOTONICITY(psi_mhd_minus,psi_mhd,psi_mhd_plus)
       end if
    end do
    if (poly .eq. 0) then
       do i = Spritz_stencil, nx - Spritz_stencil + 1
          MONOTONICITY(epsminus,eps,epsplus)
       end do
    end if

    !if (check_for_trivial_rp .eq. 0) then
    trivial_rp(:) = .false.
    !end if


    return

  end subroutine SimplePPM_1d_MHD


  /*@@
    @routine    SimplePPM_ye_1d
    @date
    @author     Added by Lorenzo Sala
    @desc
    Written modifying the homonymous routine in GRHydro.
    Note that apply_enhanced_ppm has been suppressed.
    @enddesc
    @calls
    @calledby
    @history
    @endhistory

 @@*/

subroutine SimplePPM_ye_1d(&
     nx,dx,rho,velx,vely,velz, &
     Y_e,Y_e_minus,Y_e_plus,press)

  USE Spritz_Scalars

  implicit none

  DECLARE_CCTK_PARAMETERS
  ! DECLARE_CCTK_FUNCTIONS


  CCTK_REAL, parameter :: one = 1

  CCTK_INT :: nx
  CCTK_REAL :: dx
  CCTK_REAL, dimension(nx) :: rho,velx,vely,velz
  CCTK_REAL, dimension(nx) :: Y_e,Y_e_minus,Y_e_plus


  CCTK_INT :: i,s
  CCTK_REAL, dimension(nx) :: press,dpress,tilde_flatten
  CCTK_REAL, dimension(nx) :: dmY_e, dY_e
  CCTK_REAL :: dpress2,w,flatten,dvel

  CCTK_REAL :: D2a, D2aL, D2aR, D2aC, D2aLim, rhi, daplus, daminus, D3a, D3aLL, D3aL, D3aR, D2aLL, D2aRR, D3aMin, D3aMax

  do i = 2, nx - 1
       dpress(i) = press(i+1) - press(i-1)
  end do



!!$  Average slopes delta_m(a). See (1.7) of Colella and Woodward, p.178
!!$  This is the expression for an even grid.


  do i = 2, nx - 1
        dY_e(i) = 0.5d0 * (Y_e(i+1) - Y_e(i-1))
  !        d2Y_e(i,iY_e) = (Y_e(i+1) - 2.d0 * Y_e(i) + Y_e(i-1))! / 6.d0 / dx / dx
  !    ! since we use d2Y_e only for the condition d2Y_e(i+1)*d2Y_e(i-1)<0
  !    ! the denominator is not necessary
  enddo


  !!$  Steepened slope. See (1.8) of Colella and Woodward, p.178

  do i = 2, nx - 1
     if( (Y_e(i+1) - Y_e(i)) * &
           (Y_e(i) - Y_e(i-1)) > 0.0d0 ) then
        dmY_e(i) = sign(one,dY_e(i)) * &
              min(abs(dY_e(i)), 2.0d0 * &
              abs(Y_e(i) - Y_e(i-1)), &
              2.0d0 * abs(Y_e(i+1) - Y_e(i)))
     else
        dmY_e(i) = 0.0d0
     endif
  end do

  !!$  Initial boundary states. See (1.9) of Colella and Woodward, p.178

  do i = 2, nx - 2
     Y_e_plus(i) = 0.5d0 * (Y_e(i) + Y_e(i+1)) + &
           (dmY_e(i) - dmY_e(i+1)) / 6.d0
     Y_e_minus(i+1) = Y_e_plus(i)
  enddo



!!$Discontinuity steepening. See (1.14-17) of C&W.
!!$This is the detect routine which mat be activated with the ppm_detect parameter
!!$Note that this part really also depends on the grid being even.
!!$Note also that we do not have access to the gas constant gamma.
!!$So this is just dropped from eq. (3.2) of C&W.
!!$We can get around this by just rescaling the constant k0 (ppm_k0 here).


!!$  Zone flattening. See appendix of C&W, p. 197-8.

  do i = 3, nx - 2
    dpress2 = press(i+2) - press(i-2)
    dvel = velx(i+1) - velx(i-1)
    if ( (abs(dpress(i)) >  ppm_epsilon * min(press(i-1),press(i+1))) .and. &
         (dvel < 0.d0) ) then
      w = 1.d0
    else
      w = 0.d0
    end if
    if (abs(dpress2) < ppm_small) then
      tilde_flatten(i) = 1.d0
    else
      tilde_flatten(i) = max(0.d0, 1.d0 - w * max(0.d0, ppm_omega2 * &
           (dpress(i) / dpress2 - ppm_omega1)))
    end if
  end do

  if (PPM3) then
     do i = 3, nx - 2
        flatten = tilde_flatten(i)
        Y_e_plus(i) = flatten * Y_e_plus(i) &
              + (1.d0 - flatten) * Y_e(i)
        Y_e_minus(i) = flatten * Y_e_minus(i) &
              + (1.d0 - flatten) * Y_e(i)
     end do
  else  !!$ Really implement C&W, page 197; which requires stencil 4.
     do i = 4, nx - 3
        s=int(sign(one, -dpress(i)))
        flatten = max(tilde_flatten(i), tilde_flatten(i+s))
        Y_e_plus(i) = flatten * Y_e_plus(i) + &
              (1.d0 - flatten) * Y_e(i)
        Y_e_minus(i) = flatten * Y_e_minus(i) &
              + (1.d0 - flatten) * Y_e(i)
     end do
  end if


#undef MON_WITH_LOCAL_EXTREMUM_STENCIL4
!! Monotonicity of PPM 2011 (McCorquodale & Colella 2011), Sec. 2.4.1, Eq. 23-34
!! This requires 4 stencil points.
#define MON_WITH_LOCAL_EXTREMUM_STENCIL4(aminus, a, aplus, C, C3) \
      daplus = aplus(i)-a(i)   &&\
      daminus = a(i)-aminus(i) &&\
      if (daplus*daminus .le. 0 .or. (a(i-2)-a(i))*(a(i)-a(i+2)) .le. 0) then &&\
         D2a  = - (12.0d0*a(i) - 6.0d0*(aminus(i)+aplus(i)))                  &&\
         D2aC  = (a(i-1) + a(i+1)) - 2.0d0*a(i)                               &&\
         D2aL  = (a(i-2) + a(i)  ) - 2.0d0*a(i-1)                             &&\
         D2aLL = (a(i-3) + a(i-1)) - 2.0d0*a(i-2)                             &&\
         D2aR  = (a(i)   + a(i+2)) - 2.0d0*a(i+1)                             &&\
         D2aRR = (a(i+1) + a(i+3)) - 2.0d0*a(i+2)                             &&\
         D3a = D2aR - D2aC   &&\
         D3aL = D2aC - D2aL  &&\
         D3aR = D2aRR - D2aR &&\
         D3aLL = D2aL - D2aLL &&\
         if (sign(one, D2a) .eq. sign(one, D2aC) .and. sign(one, D2a) .eq. sign(one, D2aL) .and. sign(one, D2a) .eq. sign(one, D2aR)) then &&\
            D2aLim = sign(one, D2a) * min(C*abs(D2aL), C*abs(D2aR), C*abs(D2aC), abs(D2a))  &&\
         else                                                      &&\
            D2aLim = 0  &&\
         end if                                                    &&\
         if (abs(D2a) .le. 1.d-12*max(abs(a(i-2)), abs(a(i-1)), abs(a(i)), abs(a(i+1)), abs(a(i+2)))) then  &&\
            rhi = 0  &&\
         else        &&\
            rhi = D2aLim / D2a   &&\
         endif  &&\
         if (.not. (rhi .ge. 1.0d0 - 1.d-12)) then   &&\
            D3aMin = min(D3aLL, D3aL, D3a, D3aR)   &&\
            D3aMax = max(D3aLL, D3aL, D3a, D3aR)   &&\
            if (C3 * max(abs(D3aMin), abs(D3aMax)) .le. D3aMax-D3aMin) then &&\
                  if (daplus*daminus .lt. 0) then        &&\
                     aplus(i)  = a(i) + daplus * rhi     &&\
                     aminus(i) = a(i) - daminus * rhi    &&\
                  else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then   &&\
                     aminus(i)  = a(i) - (2.0d0*(1.0d0-rhi)*daplus + rhi*daminus)   &&\
                  else if (abs(daplus) .ge. 2.0d0*abs(daminus)) then &&\
                     aplus(i)  = a(i) + (2.0d0*(1.0d0-rhi)*daminus + rhi*daplus)   &&\
                  endif &&\
            endif &&\
         endif &&\
      else                                                         &&\
            if (abs(daplus) .ge. 2.0d0*abs(daminus)) then     &&\
               aplus(i) = a(i) + 2.0d0*daminus               &&\
            else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then &&\
               aminus(i) = a(i) - 2.0d0*daplus               &&\
            end if                                                    &&\
      endif

#undef MON_WITH_LOCAL_EXTREMUM
!! Monotonicity of PPM 2011 (McCorquodale & Colella 2011), Sec. 2.4.1, Eq. 23-34.
!! This does not use the check for deviations from a cubic. Thus, it gets away with only 3 stencil points!
#define MON_WITH_LOCAL_EXTREMUM(aminus, a, aplus, C) \
      daplus = aplus(i)-a(i)   &&\
      daminus = a(i)-aminus(i) &&\
      if (daplus*daminus .le. 0 .or. (a(i-2)-a(i))*(a(i)-a(i+2)) .le. 0) then &&\
         D2a  = - (12.0d0*a(i) - 6.0d0*(aminus(i)+aplus(i)))                  &&\
         D2aC = (a(i-1) + a(i+1)) - 2.0d0*a(i)                                &&\
         D2aL = (a(i-2) + a(i)  ) - 2.0d0*a(i-1)                              &&\
         D2aR = (a(i)   + a(i+2)) - 2.0d0*a(i+1)                              &&\
         if (sign(one, D2a) .eq. sign(one, D2aC) .and. sign(one, D2a) .eq. sign(one, D2aL) .and. sign(one, D2a) .eq. sign(one, D2aR)) then &&\
            D2aLim = sign(one, D2a) * min(C*abs(D2aL), C*abs(D2aR), C*abs(D2aC), abs(D2a))  &&\
         else                                                      &&\
            D2aLim = 0  &&\
         end if                                                    &&\
         if (abs(D2a) .le. 1.d-12*max(abs(a(i-2)), abs(a(i-1)), abs(a(i)), abs(a(i+1)), abs(a(i+2)))) then  &&\
            rhi = 0  &&\
         else        &&\
            rhi = D2aLim / D2a   &&\
         endif  &&\
         if (.not. (rhi .ge. 1.0d0 - 1.d-12)) then   &&\
               if (daplus*daminus .lt. 0) then        &&\
                  aplus(i)  = a(i) + daplus * rhi     &&\
                  aminus(i) = a(i) - daminus * rhi    &&\
               else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then   &&\
                  aminus(i)  = a(i) - (2.0d0*(1.0d0-rhi)*daplus + rhi*daminus)   &&\
               else if (abs(daplus) .ge. 2.0d0*abs(daminus)) then &&\
                  aplus(i)  = a(i) + (2.0d0*(1.0d0-rhi)*daminus + rhi*daplus)   &&\
               endif &&\
         endif      &&\
      else                                                         &&\
            if (abs(daplus) .ge. 2.0d0*abs(daminus)) then     &&\
               aplus(i) = a(i) - 2.0d0*(aminus(i) - a(i))               &&\
            else if (abs(daminus) .ge. 2.0d0*abs(daplus)) then &&\
               aminus(i) = a(i) - 2.0d0*(aplus(i) - a(i))               &&\
            end if  &&\
      endif



!!$ Monotonicity. See (1.10) of C&W.
do i = Spritz_stencil, nx - Spritz_stencil + 1
   if (((Y_e_plus(i)-Y_e(i))*      &
         (Y_e(i)-Y_e_minus(i)) .le. 0.d0)) then
      Y_e_minus(i) = Y_e(i)
      Y_e_plus(i) = Y_e(i)
   else if ((Y_e_plus(i) - Y_e_minus(i)) * (Y_e(i) - 0.5d0 * &
         (Y_e_plus(i) + Y_e_minus(i))) > &
         (Y_e_plus(i) - Y_e_minus(i))**2 / 6.d0) then
      Y_e_minus(i) = 3.d0 * Y_e(i) - 2.0d0 * Y_e_plus(i)
   else if ((Y_e_plus(i) - Y_e_minus(i)) * (Y_e(i) - 0.5d0 * &
         (Y_e_plus(i) + Y_e_minus(i))) <  &
         -(Y_e_plus(i) - Y_e_minus(i))**2 / 6.0d0 ) then
      Y_e_plus(i) = 3.d0 * Y_e(i) - 2.d0 * Y_e_minus(i)
   end if
end do

return

end subroutine SimplePPM_ye_1d
