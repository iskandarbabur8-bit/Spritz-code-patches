  /*@@
   @file      primitive2conservative
   @date      Thu Jan  11 11:03:32 2002
   @author    Pedro Montero, Ian Hawke
   @desc
   Primitive to conservative routine
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
   @routine   primitive2conservative.f90
   @date       Thu Jan 11 11:03:32 2002
   @author     Pedro Montero, Ian Hawke
   @desc
   Converts primitive to conserved variables for the boundary extended data.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/

subroutine primitive2conservative(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  integer :: i, j, k, mhd_stc
  CCTK_REAL :: gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,avg_detl,&
       gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,avg_detr,psi4l,psi4r, &
       alpl,alpr,betaxl,betaxr,betayl,betayr,betazl,betazr,v2p,v2m,v2

  CCTK_REAL :: xtemp(1)

  character(len=200) warnline

  !we need a larger stencil, actually because of Constrained transport and high-order methods
  if (high_order .ne. 0) then
  		if (approximation_order == 4) mhd_stc = 1
  		if (approximation_order == 6) mhd_stc = 2
  else
  		mhd_stc = 1
  end if


  !$OMP PARALLEL DO PRIVATE(gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,avg_detl,&
  !$OMP gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,avg_detr,psi4l,psi4r, &
  !$OMP alpl,alpr,betaxl,betaxr,betayl,betayr,betazl,betazr,v2p,v2m, v2, xtemp,&
  !$OMP k, i, j, warnline)
  do k = Spritz_stencil-mhd_stc,cctk_lsh(3)-Spritz_stencil+1+mhd_stc
    do j = Spritz_stencil-mhd_stc,cctk_lsh(2)-Spritz_stencil+1+mhd_stc
      do i = Spritz_stencil-mhd_stc,cctk_lsh(1)-Spritz_stencil+1+mhd_stc

        gxxl = 0.5d0 * (gxx(i,j,k) + gxx(i-xoffset,j-yoffset,k-zoffset))
        gxyl = 0.5d0 * (gxy(i,j,k) + gxy(i-xoffset,j-yoffset,k-zoffset))
        gxzl = 0.5d0 * (gxz(i,j,k) + gxz(i-xoffset,j-yoffset,k-zoffset))
        gyyl = 0.5d0 * (gyy(i,j,k) + gyy(i-xoffset,j-yoffset,k-zoffset))
        gyzl = 0.5d0 * (gyz(i,j,k) + gyz(i-xoffset,j-yoffset,k-zoffset))
        gzzl = 0.5d0 * (gzz(i,j,k) + gzz(i-xoffset,j-yoffset,k-zoffset))
        gxxr = 0.5d0 * (gxx(i,j,k) + gxx(i+xoffset,j+yoffset,k+zoffset))
        gxyr = 0.5d0 * (gxy(i,j,k) + gxy(i+xoffset,j+yoffset,k+zoffset))
        gxzr = 0.5d0 * (gxz(i,j,k) + gxz(i+xoffset,j+yoffset,k+zoffset))
        gyyr = 0.5d0 * (gyy(i,j,k) + gyy(i+xoffset,j+yoffset,k+zoffset))
        gyzr = 0.5d0 * (gyz(i,j,k) + gyz(i+xoffset,j+yoffset,k+zoffset))
        gzzr = 0.5d0 * (gzz(i,j,k) + gzz(i+xoffset,j+yoffset,k+zoffset))

        alpl = 0.5d0 * (alp(i,j,k) + alp(i-xoffset,j-yoffset,k-zoffset))
        alpr = 0.5d0 * (alp(i,j,k) + alp(i+xoffset,j+yoffset,k+zoffset))
        betaxl = 0.5d0 * (betax(i,j,k) + betax(i-xoffset,j-yoffset,k-zoffset))
        betaxr = 0.5d0 * (betax(i,j,k) + betax(i+xoffset,j+yoffset,k+zoffset))
        betayl = 0.5d0 * (betay(i,j,k) + betay(i-xoffset,j-yoffset,k-zoffset))
        betayr = 0.5d0 * (betay(i,j,k) + betay(i+xoffset,j+yoffset,k+zoffset))
        betazl = 0.5d0 * (betaz(i,j,k) + betaz(i-xoffset,j-yoffset,k-zoffset))
        betazr = 0.5d0 * (betaz(i,j,k) + betaz(i+xoffset,j+yoffset,k+zoffset))


        if (.not.(conformal_state .eq. 0)) then
          psi4r = (0.5d0*(psi(i,j,k)+psi(i+xoffset,j+yoffset,k+zoffset)))**4
          psi4l = (0.5d0*(psi(i,j,k)+psi(i-xoffset,j-yoffset,k-zoffset)))**4
          gxxr = gxxr * psi4r
          gxyr = gxyr * psi4r
          gxzr = gxzr * psi4r
          gyyr = gyyr * psi4r
          gyzr = gyzr * psi4r
          gzzr = gzzr * psi4r
          gxxl = gxxl * psi4l
          gxyl = gxyl * psi4l
          gxzl = gxzl * psi4l
          gyyl = gyyl * psi4l
          gyzl = gyzl * psi4l
          gzzl = gzzl * psi4l
        end if

        if ( zero_slope_fix /= 0 ) then !start of the check on v^2 (used only if zero_slope_fix is equal to "yes")

           !Check the values of v^2_plus and v^2_minus to see if they exceed the velocity of light

           v2p =          gxxr*velxplus(i,j,k)*velxplus(i,j,k) &
                & +       gyyr*velyplus(i,j,k)*velyplus(i,j,k) &
                & +       gzzr*velzplus(i,j,k)*velzplus(i,j,k) &
                & + 2.0d0*gxyr*velxplus(i,j,k)*velyplus(i,j,k) &
                & + 2.0d0*gxzr*velxplus(i,j,k)*velzplus(i,j,k) &
                & + 2.0d0*gyzr*velyplus(i,j,k)*velzplus(i,j,k)

           v2m =          gxxl*velxminus(i,j,k)*velxminus(i,j,k) &
                & +       gyyl*velyminus(i,j,k)*velyminus(i,j,k) &
                & +       gzzl*velzminus(i,j,k)*velzminus(i,j,k) &
                & + 2.0d0*gxyl*velxminus(i,j,k)*velyminus(i,j,k) &
                & + 2.0d0*gxzl*velxminus(i,j,k)*velzminus(i,j,k) &
                & + 2.0d0*gyzl*velyminus(i,j,k)*velzminus(i,j,k)

           ! if ((v2p>=1.0d0).or.(v2m>=1.0d0)) then
           IF_V_CHK: IF ( (v2p>=1.0d0) .OR. (v2m>=1.0d0) .OR. (v2p<0.0d0) .OR. (v2m<0.0d0) ) THEN

              ! v2  =          gxx(i,j,k)*velx(i,j,k)*velx(i,j,k) &
              !      & +       gyy(i,j,k)*vely(i,j,k)*vely(i,j,k) &
              !      & +       gzz(i,j,k)*velz(i,j,k)*velz(i,j,k) &
              !      & + 2.0d0*gxy(i,j,k)*velx(i,j,k)*vely(i,j,k) &
              !      & + 2.0d0*gxz(i,j,k)*velx(i,j,k)*velz(i,j,k) &
              !      & + 2.0d0*gyz(i,j,k)*vely(i,j,k)*velz(i,j,k)

              !$OMP CRITICAL
              ! WRITE( warnline, '( "Detect: ( v^2_plus, v^2_minus, v^2, x, y, z )", 6G16.7 )' ) v2p, v2m, v2, x(i,j,k), y(i,j,k), z(i,j,k)
              WRITE( warnline, '( "Detect: ( v^2_plus, v^2_minus, x, y, z )", 5G16.7 )' ) v2p, v2m, x(i,j,k), y(i,j,k), z(i,j,k)
              CALL CCTK_WARN( CCTK_WARN_COMPLAIN, warnline )
              ! call CCTK_WARN(2,'Setting the slope equal to zero in that cell FOR ALL THE QUANTITIES (metric included)')
              !$OMP END CRITICAL

              gxxr=gxx(i,j,k)
              gxyr=gxy(i,j,k)
              gxzr=gxz(i,j,k)
              gyyr=gyy(i,j,k)
              gyzr=gyz(i,j,k)
              gzzr=gzz(i,j,k)
              gxxl=gxx(i,j,k)
              gxyl=gxy(i,j,k)
              gxzl=gxz(i,j,k)
              gyyl=gyy(i,j,k)
              gyzl=gyz(i,j,k)
              gzzl=gzz(i,j,k)

              if (.not.(conformal_state .eq. 0)) then
                 psi4r = psi(i,j,k)
                 psi4l = psi(i,j,k)
                 gxxr = gxxr * psi4r
                 gxyr = gxyr * psi4r
                 gxzr = gxzr * psi4r
                 gyyr = gyyr * psi4r
                 gyzr = gyzr * psi4r
                 gzzr = gzzr * psi4r
                 gxxl = gxxl * psi4l
                 gxyl = gxyl * psi4l
                 gxzl = gxzl * psi4l
                 gyyl = gyyl * psi4l
                 gyzl = gyzl * psi4l
                 gzzl = gzzl * psi4l
              end if

              alpr=alp(i,j,k)
              alpl=alp(i,j,k)

              betaxr=betax(i,j,k)
              betayr=betay(i,j,k)
              betazr=betaz(i,j,k)
              betaxl=betax(i,j,k)
              betayl=betay(i,j,k)
              betazl=betaz(i,j,k)

              rhoplus(i,j,k) = rho(i,j,k)
              rhominus(i,j,k) = rho(i,j,k)
              velxplus(i,j,k) = velx(i,j,k)
              velxminus(i,j,k) = velx(i,j,k)
              velyplus(i,j,k) = vely(i,j,k)
              velyminus(i,j,k) = vely(i,j,k)
              velzplus(i,j,k) = velz(i,j,k)
              velzminus(i,j,k) = velz(i,j,k)
              epsplus(i,j,k) = eps(i,j,k)
              epsminus(i,j,k) = eps(i,j,k)
              Bxplus(i,j,k) = Bx(i,j,k)
              Bxminus(i,j,k) = Bx(i,j,k)
              Byplus(i,j,k) = By(i,j,k)
              Byminus(i,j,k) = By(i,j,k)
              Bzplus(i,j,k) = Bz(i,j,k)
              Bzminus(i,j,k) = Bz(i,j,k)


           END IF IF_V_CHK

        end if !end of the check on v^2 (used only if zero_slope_fix is equal to "yes")

        call SpatialDeterminant(gxxl,gxyl,gxzl,gyyl, gyzl,gzzl,avg_detl)
        call SpatialDeterminant(gxxr,gxyr,gxzr,gyyr, gyzr,gzzr,avg_detr)

        if(Spritz_eos_handle.ne.4) then !if EoS is not a hot one of nuclear tabulated type
          call prim2con_mhd(Spritz_eos_handle, gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, &
               avg_detl,densminus(i,j,k),sxminus(i,j,k),&
               syminus(i,j,k),szminus(i,j,k),tauminus(i,j,k),rhominus(i,j,k), &
               velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k),&
               epsminus(i,j,k),pressminus(i,j,k),w_lorentzminus(i, j, k), &
               Bxminus(i,j,k),Byminus(i,j,k),Bzminus(i,j,k),&
               dBxminus(i,j,k),dByminus(i,j,k),dBzminus(i,j,k),&
               alpl,betaxl,betayl,betazl,&
               Spritz_coldtab_rhomin,Spritz_coldtab_press_lo,&
               Spritz_coldtab_eps_lo,Spritz_coldtab_gamma_lo,&
               Spritz_coldtab_gammath,Spritz_coldtab_thfac)
          call prim2con_mhd(Spritz_eos_handle, gxxr,gxyr,gxzr,gyyr,gyzr,gzzr, &
               avg_detr, densplus(i,j,k),sxplus(i,j,k),&
               syplus(i,j,k),szplus(i,j ,k),tauplus(i,j,k),&
               rhoplus(i,j,k),velxplus(i,j,k),velyplus(i,j,k),&
               velzplus(i,j,k),epsplus(i,j,k),pressplus(i,j,k),&
               w_lorentzplus(i,j,k),Bxplus(i,j,k),Byplus(i,j,k),&
               Bzplus(i,j,k),dBxplus(i,j,k),dByplus(i,j,k),&
               dBzplus(i,j,k),alpr,betaxr,betayr,betazr,&
               Spritz_coldtab_rhomin,Spritz_coldtab_press_lo,&
               Spritz_coldtab_eps_lo,Spritz_coldtab_gamma_lo,&
               Spritz_coldtab_gammath,Spritz_coldtab_thfac)
        else !if using hot eos
          xtemp(1) = 0.5d0*(temperature(i,j,k) + &
                     temperature(i-xoffset,j-yoffset,k-zoffset))
          if (y_e_minus(i,j,k) .le. 0.0d0 .or. y_e_plus(i,j,k) .le. 0.0d0) then
            !$OMP CRITICAL
            write(warnline,'(a100,7g15.6)') '(y_e_minus,y_e_plus,x,y,z,rho)', y_e(i,j,k), y_e_minus(i,j,k), y_e_plus(i,j,k), x(i,j,k),y(i,j,k),z(i,j,k),rho(i,j,k)
            call CCTK_WARN(1, warnline)
            !$OMP END CRITICAL
          endif
          call prim2con_mhd_hot(Spritz_eos_handle, Constrain_Temp_now,&
               Spritz_reflevel, x(i,j,k),y(i,j,k),z(i,j,k), &
               gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, &
               avg_detl,densminus(i,j,k),sxminus(i,j,k),&
               syminus(i,j,k),szminus(i,j,k),tauminus(i,j,k),rhominus(i,j,k), &
               velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k),&
               epsminus(i,j,k),pressminus(i,j,k),w_lorentzminus(i, j, k), &
               Bxminus(i,j,k),Byminus(i,j,k),Bzminus(i,j,k),&
               dBxminus(i,j,k),dByminus(i,j,k),dBzminus(i,j,k),&
               alpl,betaxl,betayl,betazl, xtemp, y_e_minus(i,j,k))
          xtemp(1) = 0.5d0*(temperature(i,j,k) + &
                     temperature(i+xoffset,j+yoffset,k+zoffset))
          call prim2con_mhd_hot(Spritz_eos_handle, Constrain_Temp_now,&
               Spritz_reflevel, x(i,j,k),y(i,j,k),z(i,j,k), &
               gxxr,gxyr,gxzr,gyyr,gyzr,gzzr, &
               avg_detr, densplus(i,j,k),sxplus(i,j,k),&
               syplus(i,j,k),szplus(i,j ,k),tauplus(i,j,k),&
               rhoplus(i,j,k),velxplus(i,j,k),velyplus(i,j,k),&
               velzplus(i,j,k),epsplus(i,j,k),pressplus(i,j,k),&
               w_lorentzplus(i,j,k),Bxplus(i,j,k),Byplus(i,j,k),&
               Bzplus(i,j,k),dBxplus(i,j,k),dByplus(i,j,k),&
               dBzplus(i,j,k),alpr,betaxr,betayr,betazr, xtemp, y_e_plus(i,j,k))
          end if

      end do
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine primitive2conservative


 /*@@
   @routine    prim2con_mhd
   @date
   @author     Bruno Giacomazzo
   @desc
   Converts from primitive to conservative at a single point
   (GRMHD version)
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/

subroutine prim2con_mhd(handle, gxx, gxy, gxz, gyy, gyz, gzz, det, ddens, &
     dsx, dsy, dsz, dtau , drho, dvelx, dvely, dvelz, deps, dpress, w, &
     Bx, By, Bz, dBx, dBy, dBz, alp, betax, betay, betaz, &
     coldtab_rhomin, coldtab_press_lo, coldtab_eps_lo, coldtab_gamma_lo, &
     coldtab_gammath, coldtab_thfac)

  implicit none


  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL :: gxx, gxy, gxz, gyy, gyz, gzz, det
  CCTK_REAL :: ddens, dsx, dsy, dsz, dtau, drho, dvelx, dvely, dvelz,&
       deps, dpress, w, vlowx, vlowy, vlowz, dBx, dBy, dBz, &
       Blowx, Blowy, Blowz, bst, bslowx, bslowy, bslowz, bs2, &
       Bx, By, Bz, alp, betax, betay, betaz, betaxlow, betaylow, betazlow, &
       bsx, bsy, bsz
  CCTK_REAL :: coldtab_rhomin, coldtab_press_lo, coldtab_eps_lo, coldtab_gamma_lo
  CCTK_REAL :: coldtab_gammath, coldtab_thfac
  CCTK_INT :: handle, ierr

!#ifdef _EOS_BASE_INC_
!#undef _EOS_BASE_INC_
!#endif
!#include "EOS_Base.inc"

  !!!!!!
  ! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1)

  n = 1;keytemp = 0;anyerr = 0;keyerr = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
  ! end EOS Omni vars
  !!!!!

  w = 1.d0 / sqrt(1.d0 - (gxx*dvelx*dvelx + gyy*dvely*dvely + gzz &
       *dvelz*dvelz + 2*gxy*dvelx*dvely + 2*gxz*dvelx *dvelz + 2*gyz&
       *dvely*dvelz))


  !dpress = EOS_Pressure(handle, drho, deps)

  !!!!!!
  ! Routed through Spritz_coldtab_press (Spritz_coldtab_eos.c) instead of
  ! EOS_Omni_press directly: this is the hot interface-pressure calculation
  ! for every reconstructed cell face, every timestep. Below the
  ! cold_tabulated table's rho_min, uses the table's real local Gamma
  ! instead of EOS_Omni's own hardcoded Gamma=2; otherwise falls through to
  ! the unchanged real EOS_Omni_press call.
  xrho = drho
  xeps = deps
  call Spritz_coldtab_press(handle,keytemp,Eos_rf_prec,n,&
                      xrho,xeps,xtemp,xye,xpress,keyerr,anyerr,&
                      coldtab_rhomin,coldtab_press_lo,&
                      coldtab_eps_lo,coldtab_gamma_lo,&
                      coldtab_gammath,coldtab_thfac)
  dpress = xpress(1)
  !!!!!!

  !v_j
  vlowx = gxx*dvelx + gxy*dvely + gxz*dvelz
  vlowy = gxy*dvelx + gyy*dvely + gyz*dvelz
  vlowz = gxz*dvelx + gyz*dvely + gzz*dvelz

  !beta_j
  betaxlow = gxx*betax + gxy*betay + gxz*betaz
  betaylow = gxy*betax + gyy*betay + gyz*betaz
  betazlow = gxz*betax + gyz*betay + gzz*betaz


  !B_j
  Blowx = gxx*Bx + gxy*By + gxz*Bz
  Blowy = gxy*Bx + gyy*By + gyz*Bz
  Blowz = gxz*Bx + gyz*By + gzz*Bz

  !b^t
  bst = w*(Bx*vlowx + By*vlowy + Bz*vlowz) !this is b^0 * alp

  !b^j
  bsx = (Bx + bst*w*(dvelx-betax/alp))/w
  bsy = (By + bst*w*(dvely-betay/alp))/w
  bsz = (Bz + bst*w*(dvelz-betaz/alp))/w

  !b_j

  bslowx = betaxlow*bst/alp + gxx*bsx + gxy*bsy + gxz*bsz
  bslowy = betaylow*bst/alp + gxy*bsx + gyy*bsy + gyz*bsz
  bslowz = betazlow*bst/alp + gxz*bsx + gyz*bsy + gzz*bsz

  !b^mu b_mu
  bs2 = (Bx*Blowx + By*Blowy + Bz*Blowz + bst**2)/w**2

  ddens = sqrt(det) * drho * w
  dsx   = sqrt(det) * ((drho*(1+deps)+dpress+bs2)*w*w * vlowx - bst*bslowx)
  dsy   = sqrt(det) * ((drho*(1+deps)+dpress+bs2)*w*w * vlowy - bst*bslowy)
  dsz   = sqrt(det) * ((drho*(1+deps)+dpress+bs2)*w*w * vlowz - bst*bslowz)
  dtau  = sqrt(det) * ((drho*(1+deps)+dpress+bs2)*w*w - (dpress+0.5d0*bs2) - &
       bst**2) - ddens
  dBx = sqrt(det)*Bx
  dBy = sqrt(det)*By
  dBz = sqrt(det)*Bz

end subroutine prim2con_mhd


/*@@
  @routine    prim2con_mhd_hot
  @date
  @author     Lorenzo Sala, modifying code by Bruno Giacomazzo
  @desc
  Converts from primitive to conservative at a single point
  (GRMHD version with tabulated EOS)
  @enddesc
  @calls
  @calledby
  @history

  @endhistory

@@*/
subroutine prim2con_mhd_hot(handle, Constrain_Temp_now, Spritz_reflevel, x, y, z, &
     gxx, gxy, gxz, gyy, gyz, gzz, det, ddens, &
     dsx, dsy, dsz, dtau , drho, dvelx, dvely, dvelz, deps, dpress, w, &
     Bx, By, Bz, dBx, dBy, dBz, alp, betax, betay, betaz, temp, ye)

  implicit none


  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL :: gxx, gxy, gxz, gyy, gyz, gzz, det
  CCTK_REAL :: ddens, dsx, dsy, dsz, dtau, &
       dBvcx, dBvcy, dBvcz, &
       vlowx, vlowy, vlowz, dBx, dBy, dBz, &
       Blowx, Blowy, Blowz, bst, bslowx, bslowy, bslowz, bs2, &
       Bx, By, Bz, betaxlow, betaylow, betazlow, &
       bsx, bsy, bsz, temp, ye, drho, deps, dpress
  CCTK_REAL :: w, alp, betax, betay, betaz, &
               dvelx, dvely, dvelz, x, y, z
  CCTK_REAL, dimension(1) :: xtemp, xye, yein, xrho, xeps, xpress
  CCTK_INT :: handle, Spritz_reflevel, Constrain_Temp_now
  character(len=256) NaN_WarnLine
  character(len=512) warnline


! begin EOS Omni vars
  CCTK_INT  :: n, keytemp, anyerr, keyerr(1)
  CCTK_REAL :: temp0(1)
  n = 1; keytemp = 0; anyerr = 0; keyerr(1) = 0
! end EOS Omni vars

  temp0 = temp
  w = 1.d0 / sqrt(1.d0 - (gxx*dvelx*dvelx + gyy*dvely*dvely + gzz &
       *dvelz*dvelz + 2*gxy*dvelx*dvely + 2*gxz*dvelx *dvelz + 2*gyz&
       *dvely*dvelz))


   !!$ BEGIN: Check for NaN value
     if (w .ne. w) then
       !$OMP CRITICAL
       write(NaN_WarnLine,'(a100,6g15.6)') 'NaN produced in sqrt(): (gxx,gxy,gxz,gyy,gyz,gzz)', gxx, gxy, gxz, gyy, gyz, gzz
       call CCTK_WARN(1, NaN_WarnLine)
       write(NaN_WarnLine,'(a100,3g15.6)') 'NaN produced in sqrt(): (dvelx,dvely,dvelz)', dvelx, dvely, dvelz
       call CCTK_WARN(1, NaN_WarnLine)
       write(NaN_WarnLine,'(a100,3g15.6)') 'NaN produced in sqrt(): (x,y,z)', x, y, z
       call CCTK_WARN(1, NaN_WarnLine)
       !write(NaN_WarnLine,'(a100,g15.6)') 'NaN produced in sqrt(): w', w
       !call CCTK_WARN(Spritz_NaN_verbose, NaN_WarnLine)
       !$OMP END CRITICAL
     endif
   !!$ END: Check for NaN value

     xye = max(min(ye,Spritz_Y_e_max),Spritz_Y_e_min)
     yein = xye
     xrho = drho
     xeps = deps
     xtemp = max(min(temp,Spritz_max_temp),Spritz_hot_atmo_temp)

     keytemp = 0 !we are looking for temperature and pressure at given internal energy.
     if(Constrain_Temp_now.eq.1) keytemp = 1 !but, if we have to constrain temperature, we have to modify eps at given temperature
     call EOS_Omni_press(handle,keytemp,Eos_rf_prec,n,&
          xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)

     ! error handling
     if(anyerr.ne.0) then
        if(keyerr(1).eq.668) then
          !$OMP CRITICAL
          write(warnline, "(a33,3g15.6)") "Temperature too low at (x,y,z) = ", x, y, z
          call CCTK_WARN(2, warnline)
          write(warnline, "(a7,i5,i5)") "Code: ", keyerr(1)
          call CCTK_WARN(2, warnline)
          write(warnline, "(a48)") "Setting xtemp and temp to atmosphere temperature"
          call CCTK_WARN(2, warnline)
          !$OMP END CRITICAL
          xtemp = Spritz_hot_atmo_temp
          temp  = Spritz_hot_atmo_temp
        end if
        if(Spritz_reflevel.lt.Spritz_c2p_warn_from_reflevel) then
           ! in this case (coarse grid error that is hopefully restricted
           ! away), we use the average temperature between cells and call
           ! the EOS with keytemp=1
           keytemp=1 ! this means that we are not modifying temperature. We're simply setting up xeps and xpress.
           call EOS_Omni_press(handle,keytemp,Eos_rf_prec,n,&
                xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
           keytemp=0
           if(anyerr.ne.0) then
              !$OMP CRITICAL
              call CCTK_WARN(1,"EOS error in prim2con_hot: lev 2")
              write(warnline,"(1P10E15.6)") x,y,z
              call CCTK_WARN(1,warnline)
              write(warnline,"(1P10E15.6)") drho,deps,temp,ye,yein
              call CCTK_WARN(1,warnline)
              write(warnline,"(A7,i8)") "code: ",keyerr(1)
              call CCTK_WARN(1,warnline)
              write(warnline,"(A10,i5)") "reflevel: ", Spritz_reflevel
              call CCTK_WARN(1,warnline)
              !$OMP END CRITICAL
           endif
        else
           ! This is a way of recovering even on finer refinement levels:
           ! Use the average temperature at the interface instead of the
           ! reconstructed specific internal energy.
           !$OMP CRITICAL
           call CCTK_WARN(1,"EOS error in prim2con_hot: NOW using averaged temp!")
           write(warnline,"(1P10E15.6)") x,y,z
           call CCTK_WARN(1,warnline)
           write(warnline,"(1P10E15.6)") drho,deps,temp0,ye
           call CCTK_WARN(1,warnline)
           write(warnline,"(A7,i8)") "code: ",keyerr(1)
           call CCTK_WARN(1,warnline)
           write(warnline,"(A10,i5)") "reflevel: ", Spritz_reflevel
           call CCTK_WARN(1,warnline)
           !$OMP END CRITICAL
           keytemp=1 !not modifying temperature.
           xtemp = temp0
           call EOS_Omni_press(handle,keytemp,Eos_rf_prec,n,&
                xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
           xpress(1) = dpress !don't want to change pressure
           keytemp=0
           if(anyerr.ne.0) then
              !$OMP CRITICAL
              call CCTK_WARN(1,"EOS error in prim2con_hot")
              write(warnline,"(1P10E15.6)") x,y,z
              call CCTK_WARN(1,warnline)
              write(warnline,"(1P10E15.6)") drho,deps,temp,ye
              call CCTK_WARN(1,warnline)
              write(warnline,"(A7,i8)") "code: ",keyerr(1)
              call CCTK_WARN(1,warnline)
              write(warnline,"(A10,i5)") "reflevel: ", Spritz_reflevel
              call CCTK_WARN(1,warnline)
              call CCTK_ERROR("Aborting!!!")
              STOP
              !$OMP END CRITICAL
           endif
        endif
     endif

  dpress = xpress(1)
  deps   = xeps(1)   !this is unchanged if not constraining temperature
  temp   = xtemp(1)  !this is unchanged if     constraining temperature

  !v_j
  vlowx = gxx*dvelx + gxy*dvely + gxz*dvelz
  vlowy = gxy*dvelx + gyy*dvely + gyz*dvelz
  vlowz = gxz*dvelx + gyz*dvely + gzz*dvelz

  !beta_j
  betaxlow = gxx*betax + gxy*betay + gxz*betaz
  betaylow = gxy*betax + gyy*betay + gyz*betaz
  betazlow = gxz*betax + gyz*betay + gzz*betaz


  !B_j
  Blowx = gxx*Bx + gxy*By + gxz*Bz
  Blowy = gxy*Bx + gyy*By + gyz*Bz
  Blowz = gxz*Bx + gyz*By + gzz*Bz

  !b^t
  bst = w*(Bx*vlowx + By*vlowy + Bz*vlowz) !this is b^0 * alp

  !b^j
  bsx = (Bx + bst*w*(dvelx-betax/alp))/w
  bsy = (By + bst*w*(dvely-betay/alp))/w
  bsz = (Bz + bst*w*(dvelz-betaz/alp))/w

  !b_j

  bslowx = betaxlow*bst/alp + gxx*bsx + gxy*bsy + gxz*bsz
  bslowy = betaylow*bst/alp + gxy*bsx + gyy*bsy + gyz*bsz
  bslowz = betazlow*bst/alp + gxz*bsx + gyz*bsy + gzz*bsz

  !b^mu b_mu
  bs2 = (Bx*Blowx + By*Blowy + Bz*Blowz + bst**2)/w**2

  ddens = sqrt(det) * drho * w
  dsx   = sqrt(det) * ((drho*(1+deps)+dpress+bs2)*w*w * vlowx - bst*bslowx)
  dsy   = sqrt(det) * ((drho*(1+deps)+dpress+bs2)*w*w * vlowy - bst*bslowy)
  dsz   = sqrt(det) * ((drho*(1+deps)+dpress+bs2)*w*w * vlowz - bst*bslowz)
  dtau  = sqrt(det) * ((drho*(1+deps)+dpress+bs2)*w*w - (dpress+0.5d0*bs2) - &
       bst**2) - ddens
  dBx = sqrt(det)*Bx
  dBy = sqrt(det)*By
  dBz = sqrt(det)*Bz

end subroutine prim2con_mhd_hot

 /*@@
   @routine    Primitive2ConservativeCells
   @date       Sun Mar 10 21:16:20 2002
   @author
   @desc
   Wrapper function that converts primitive to conservative at the
     cell centres.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/


subroutine Primitive2ConservativeCells(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: i, j, k, mhd_stc
  CCTK_REAL :: det,psi4pt

  !we need a larger stencil, actually because of Constrained transport and high-order methods
  if (high_order .ne. 0) then
  		if (approximation_order == 4) mhd_stc = 1
  		if (approximation_order == 6) mhd_stc = 2
  else
  		mhd_stc = 1
  end if

  do k = Spritz_stencil-mhd_stc,cctk_lsh(3)-Spritz_stencil+1+mhd_stc
    do j = Spritz_stencil-mhd_stc,cctk_lsh(2)-Spritz_stencil+1+mhd_stc
      do i = Spritz_stencil-mhd_stc,cctk_lsh(1)-Spritz_stencil+1+mhd_stc

        if (conformal_state .eq. 0) then
          psi4pt = 1d0
          call SpatialDeterminant(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),&
               gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det)
        else
          psi4pt = psi(i,j,k)**4
          call SpatialDeterminant(psi4pt*gxx(i,j,k),psi4pt*gxy(i,j,k),&
               psi4pt*gxz(i,j,k),psi4pt*gyy(i,j,k),psi4pt*gyz(i,j,k),&
               psi4pt*gzz(i,j,k),det)
        end if

        call prim2con_mhd(Spritz_eos_handle,psi4pt*gxx(i,j,k),&
             psi4pt*gxy(i,j,k),psi4pt*gxz(i,j,k),&
             psi4pt*gyy(i,j,k),psi4pt*gyz(i,j,k),psi4pt*gzz(i,j,k),&
             det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
             tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
             eps(i,j,k),press(i,j,k),w_lorentz(i,j,k),&
             Bx(i,j,k),By(i,j,k),Bz(i,j,k),&
             dBx(i,j,k),dBy(i,j,k),dBz(i,j,k),&
             alp(i,j,k),betax(i,j,k),betay(i,j,k),betaz(i,j,k),&
             Spritz_coldtab_rhomin,Spritz_coldtab_press_lo,&
             Spritz_coldtab_eps_lo,Spritz_coldtab_gamma_lo,&
             Spritz_coldtab_gammath,Spritz_coldtab_thfac)

      end do
    end do
  end do

  if(evolve_Y_e.ne.0) then
     !$OMP PARALLEL DO PRIVATE(i, j, k)
     do k = 1,cctk_lsh(3)
        do j = 1,cctk_lsh(2)
           do i = 1,cctk_lsh(1)
              Y_e_con(i,j,k) = Y_e(i,j,k) * dens(i,j,k)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif

end subroutine Primitive2ConservativeCells


 /*@@
   @routine    Prim2ConservativePolytype
   @date       Tue Mar 19 22:52:21 2002
   @author     Ian Hawke
   @desc
   Same as first routine, only for polytropes.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/


subroutine Prim2ConservativePolytype(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  integer :: i, j, k, mhd_stc
  CCTK_REAL :: gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,avg_detl,&
       gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,avg_detr,psi4r,psi4l,&
       alpl,alpr,betaxl,betaxr,betayl,betayr,betazl,betazr,v2p,v2m, v2

  character(len=200) warnline


  !we need a larger stencil, actually because of Constrained transport and high-order methods
  if (high_order .ne. 0) then
  		if (approximation_order == 4) mhd_stc = 1
  		if (approximation_order == 6) mhd_stc = 2
  else
  		mhd_stc = 1
  end if


  !$OMP PARALLEL DO PRIVATE(gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,avg_detl,&
  !$OMP gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,avg_detr,psi4r,psi4l,&
  !$OMP alpl,alpr,betaxl,betaxr,betayl,betayr,betazl,betazr,v2p,v2m, v2, &
  !$OMP k, i, j, warnline)
  do k = Spritz_stencil-mhd_stc,cctk_lsh(3)-Spritz_stencil+1+mhd_stc
    do j = Spritz_stencil-mhd_stc,cctk_lsh(2)-Spritz_stencil+1+mhd_stc
      do i = Spritz_stencil-mhd_stc,cctk_lsh(1)-Spritz_stencil+1+mhd_stc

        gxxl = 0.5d0 * (gxx(i,j,k) + gxx(i-xoffset,j-yoffset,k-zoffset))
        gxyl = 0.5d0 * (gxy(i,j,k) + gxy(i-xoffset,j-yoffset,k-zoffset))
        gxzl = 0.5d0 * (gxz(i,j,k) + gxz(i-xoffset,j-yoffset,k-zoffset))
        gyyl = 0.5d0 * (gyy(i,j,k) + gyy(i-xoffset,j-yoffset,k-zoffset))
        gyzl = 0.5d0 * (gyz(i,j,k) + gyz(i-xoffset,j-yoffset,k-zoffset))
        gzzl = 0.5d0 * (gzz(i,j,k) + gzz(i-xoffset,j-yoffset,k-zoffset))
        gxxr = 0.5d0 * (gxx(i,j,k) + gxx(i+xoffset,j+yoffset,k+zoffset))
        gxyr = 0.5d0 * (gxy(i,j,k) + gxy(i+xoffset,j+yoffset,k+zoffset))
        gxzr = 0.5d0 * (gxz(i,j,k) + gxz(i+xoffset,j+yoffset,k+zoffset))
        gyyr = 0.5d0 * (gyy(i,j,k) + gyy(i+xoffset,j+yoffset,k+zoffset))
        gyzr = 0.5d0 * (gyz(i,j,k) + gyz(i+xoffset,j+yoffset,k+zoffset))
        gzzr = 0.5d0 * (gzz(i,j,k) + gzz(i+xoffset,j+yoffset,k+zoffset))

        alpl = 0.5d0 * (alp(i,j,k) + alp(i-xoffset,j-yoffset,k-zoffset))
        alpr = 0.5d0 * (alp(i,j,k) + alp(i+xoffset,j+yoffset,k+zoffset))
        betaxl = 0.5d0 * (betax(i,j,k) + betax(i-xoffset,j-yoffset,k-zoffset))
        betaxr = 0.5d0 * (betax(i,j,k) + betax(i+xoffset,j+yoffset,k+zoffset))
        betayl = 0.5d0 * (betay(i,j,k) + betay(i-xoffset,j-yoffset,k-zoffset))
        betayr = 0.5d0 * (betay(i,j,k) + betay(i+xoffset,j+yoffset,k+zoffset))
        betazl = 0.5d0 * (betaz(i,j,k) + betaz(i-xoffset,j-yoffset,k-zoffset))
        betazr = 0.5d0 * (betaz(i,j,k) + betaz(i+xoffset,j+yoffset,k+zoffset))

        if (.not.(conformal_state .eq. 0)) then
          psi4r = (0.5d0*(psi(i,j,k)+psi(i+xoffset,j+yoffset,k+zoffset)))**4
          psi4l = (0.5d0*(psi(i,j,k)+psi(i-xoffset,j-yoffset,k-zoffset)))**4
          gxxr = gxxr * psi4r
          gxyr = gxyr * psi4r
          gxzr = gxzr * psi4r
          gyyr = gyyr * psi4r
          gyzr = gyzr * psi4r
          gzzr = gzzr * psi4r
          gxxl = gxxl * psi4l
          gxyl = gxyl * psi4l
          gxzl = gxzl * psi4l
          gyyl = gyyl * psi4l
          gyzl = gyzl * psi4l
          gzzl = gzzl * psi4l
        end if

        if (zero_slope_fix /= 0) then !start of the check on v^2 (used only if zero_slope_fix is equal to "yes")

           !Check the values of v^2_plus and v^2_minus to see if they exceed the velocity of light

           v2p =          gxxr*velxplus(i,j,k)*velxplus(i,j,k) &
                & +       gyyr*velyplus(i,j,k)*velyplus(i,j,k) &
                & +       gzzr*velzplus(i,j,k)*velzplus(i,j,k) &
                & + 2.0d0*gxyr*velxplus(i,j,k)*velyplus(i,j,k) &
                & + 2.0d0*gxzr*velxplus(i,j,k)*velzplus(i,j,k) &
                & + 2.0d0*gyzr*velyplus(i,j,k)*velzplus(i,j,k)

           v2m =          gxxl*velxminus(i,j,k)*velxminus(i,j,k) &
                & +       gyyl*velyminus(i,j,k)*velyminus(i,j,k) &
                & +       gzzl*velzminus(i,j,k)*velzminus(i,j,k) &
                & + 2.0d0*gxyl*velxminus(i,j,k)*velyminus(i,j,k) &
                & + 2.0d0*gxzl*velxminus(i,j,k)*velzminus(i,j,k) &
                & + 2.0d0*gyzl*velyminus(i,j,k)*velzminus(i,j,k)

           ! if ((v2p>=1.0d0).or.(v2m>=1.0d0)) then
           IF_V_CHK: IF ( (v2p>=1.0d0) .OR. (v2m>=1.0d0) .OR. (v2p<0.0d0) .OR. (v2m<0.0d0) ) THEN

              ! v2 =           gxx(i,j,k)*velx(i,j,k)*velx(i,j,k) &
              !      & +       gyy(i,j,k)*vely(i,j,k)*vely(i,j,k) &
              !      & +       gzz(i,j,k)*velz(i,j,k)*velz(i,j,k) &
              !      & + 2.0d0*gxy(i,j,k)*velx(i,j,k)*vely(i,j,k) &
              !      & + 2.0d0*gxz(i,j,k)*velx(i,j,k)*velz(i,j,k) &
              !      & + 2.0d0*gyz(i,j,k)*vely(i,j,k)*velz(i,j,k)

              !$OMP CRITICAL
              ! WRITE( warnline, '( "Detect: ( v^2_plus, v^2_minus, v^2, x, y, z )", 6G16.7 )' ) v2p, v2m, v2, x(i,j,k), y(i,j,k), z(i,j,k)
              WRITE( warnline, '( "Detect: ( v^2_plus, v^2_minus, x, y, z )", 5G16.7 )' ) v2p, v2m, x(i,j,k), y(i,j,k), z(i,j,k)
              CALL CCTK_WARN( CCTK_WARN_COMPLAIN, warnline )
              ! call CCTK_WARN(2,'Setting the slope equal to zero in that cell FOR ALL THE QUANTITIES (metric included)')
              !$OMP END CRITICAL

              gxxr=gxx(i,j,k)
              gxyr=gxy(i,j,k)
              gxzr=gxz(i,j,k)
              gyyr=gyy(i,j,k)
              gyzr=gyz(i,j,k)
              gzzr=gzz(i,j,k)
              gxxl=gxx(i,j,k)
              gxyl=gxy(i,j,k)
              gxzl=gxz(i,j,k)
              gyyl=gyy(i,j,k)
              gyzl=gyz(i,j,k)
              gzzl=gzz(i,j,k)

              if (.not.(conformal_state .eq. 0)) then
                 psi4r = psi(i,j,k)
                 psi4l = psi(i,j,k)
                 gxxr = gxxr * psi4r
                 gxyr = gxyr * psi4r
                 gxzr = gxzr * psi4r
                 gyyr = gyyr * psi4r
                 gyzr = gyzr * psi4r
                 gzzr = gzzr * psi4r
                 gxxl = gxxl * psi4l
                 gxyl = gxyl * psi4l
                 gxzl = gxzl * psi4l
                 gyyl = gyyl * psi4l
                 gyzl = gyzl * psi4l
                 gzzl = gzzl * psi4l
              end if

              alpr=alp(i,j,k)
              alpl=alp(i,j,k)

              betaxr=betax(i,j,k)
              betayr=betay(i,j,k)
              betazr=betaz(i,j,k)
              betaxl=betax(i,j,k)
              betayl=betay(i,j,k)
              betazl=betaz(i,j,k)

              rhoplus(i,j,k) = rho(i,j,k)
              rhominus(i,j,k) = rho(i,j,k)
              velxplus(i,j,k) = velx(i,j,k)
              velxminus(i,j,k) = velx(i,j,k)
              velyplus(i,j,k) = vely(i,j,k)
              velyminus(i,j,k) = vely(i,j,k)
              velzplus(i,j,k) = velz(i,j,k)
              velzminus(i,j,k) = velz(i,j,k)
              epsplus(i,j,k) = eps(i,j,k)
              epsminus(i,j,k) = eps(i,j,k)
              Bxplus(i,j,k) = Bx(i,j,k)
              Bxminus(i,j,k) = Bx(i,j,k)
              Byplus(i,j,k) = By(i,j,k)
              Byminus(i,j,k) = By(i,j,k)
              Bzplus(i,j,k) = Bz(i,j,k)
              Bzminus(i,j,k) = Bz(i,j,k)

           END IF IF_V_CHK

        end if !end of the check on v^2 (used only if zero_slope_fix is equal to "yes")

        call SpatialDeterminant(gxxl,gxyl,gxzl,gyyl, gyzl,gzzl,avg_detl)
        call SpatialDeterminant(gxxr,gxyr,gxzr,gyyr, gyzr,gzzr,avg_detr)

        call prim2conpolytype_mhd(Spritz_eos_handle, gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, &
             avg_detl,densminus(i,j,k),sxminus(i,j,k),&
             syminus(i,j,k),szminus(i,j,k),tauminus(i,j,k),rhominus(i,j,k), &
             velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k),&
             epsminus(i,j,k),pressminus(i,j,k),w_lorentzminus(i, j, k), &
             Bxminus(i,j,k),Byminus(i,j,k),Bzminus(i,j,k),&
             dBxminus(i,j,k),dByminus(i,j,k),dBzminus(i,j,k), &
             alpl,betaxl,betayl,betazl)

        call prim2conpolytype_mhd(Spritz_eos_handle, gxxr,gxyr,gxzr,gyyr,gyzr,gzzr, &
             avg_detr, densplus(i,j,k),sxplus(i,j,k),&
             syplus(i,j,k),szplus(i,j ,k),tauplus(i,j,k),&
             rhoplus(i,j,k),velxplus(i,j,k),velyplus(i,j,k),&
             velzplus(i,j,k),epsplus(i,j,k),pressplus(i,j,k),&
             w_lorentzplus(i,j,k),Bxplus(i,j,k),Byplus(i,j,k),&
             Bzplus(i,j,k),dBxplus(i,j,k),dByplus(i,j,k),&
             dBzplus(i,j,k),alpr,betaxr,betayr,betazr)

      end do
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine Prim2ConservativePolytype


 /*@@
   @routine    prim2conpolytype_mhd
   @date
   @author     Bruno Giacomazzo
   @desc
   Converts from primitive to conservative at a single point
   (GRMHD version)
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/


subroutine prim2conpolytype_mhd(handle, gxx, gxy, gxz, gyy, gyz, gzz, det, ddens, &
     dsx, dsy, dsz, dtau , drho, dvelx, dvely, dvelz, deps, dpress, w, &
     Bx, By, Bz, dBx, dBy, dBz, alp, betax, betay, betaz)

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL :: gxx, gxy, gxz, gyy, gyz, gzz, det
  CCTK_REAL :: ddens, dsx, dsy, dsz, dtau, drho, dvelx, dvely, dvelz,&
       deps, dpress, w, vlowx, vlowy, vlowz, dBx, dBy, dBz,&
       Blowx, Blowy, Blowz, bst, bslowx, bslowy, bslowz, bs2, &
       Bx, By, Bz, alp, betax, betay, betaz, betaxlow, betaylow, betazlow, &
       bsx, bsy, bsz
  CCTK_INT :: handle, ierr

!#ifdef _EOS_BASE_INC_
!#undef _EOS_BASE_INC_
!#endif
!#include "EOS_Base.inc"

  !!!!!!
  ! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1)

  n = 1;keytemp = 0;anyerr = 0;keyerr = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
  ! end EOS Omni vars
  !!!!!

  w = 1.d0 / sqrt(1.d0 - (gxx*dvelx*dvelx + gyy*dvely*dvely + gzz &
       *dvelz*dvelz + 2*gxy*dvelx*dvely + 2*gxz*dvelx*dvelz + 2*gyz &
       *dvely*dvelz))

  !dpress = EOS_Pressure(handle, drho, deps)
  !deps = EOS_SpecificIntEnergy(handle, drho, dpress)

  !!!!!!
  ! NOTE: this used to also call EOS_Omni_EpsFromPress(handle,...,xpress,xeps,...)
  ! right after, overwriting deps with the "eps" implied by the press we just
  ! computed FROM that same deps. For any self-consistent EOS that round trip
  ! is a no-op: dpress was itself computed as a function of deps in the call
  ! below, so inverting it must return exactly deps again. EOS_Omni_EpsFromPress
  ! has no case(5) (cold_tabulated) and hard-aborts the run if this EOS is in
  ! use, so we drop the redundant inverse call entirely and leave deps as-is --
  ! mathematically identical result, and it now also works for cold_tabulated.
  xrho = drho
  xeps = deps
  call EOS_Omni_press(handle,keytemp,Eos_rf_prec,n,&
                      xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
  dpress = xpress(1)
  !!!!!!

  !v_j
  vlowx = gxx*dvelx + gxy*dvely + gxz*dvelz
  vlowy = gxy*dvelx + gyy*dvely + gyz*dvelz
  vlowz = gxz*dvelx + gyz*dvely + gzz*dvelz

  !beta_j
  betaxlow = gxx*betax + gxy*betay + gxz*betaz
  betaylow = gxy*betax + gyy*betay + gyz*betaz
  betazlow = gxz*betax + gyz*betay + gzz*betaz


  !B_j
  Blowx = gxx*Bx + gxy*By + gxz*Bz
  Blowy = gxy*Bx + gyy*By + gyz*Bz
  Blowz = gxz*Bx + gyz*By + gzz*Bz

  !b^t
  bst = w*(Bx*vlowx + By*vlowy + Bz*vlowz) !this is b^0 * alp

  !b^j
  bsx = (Bx + bst*w*(dvelx-betax/alp))/w
  bsy = (By + bst*w*(dvely-betay/alp))/w
  bsz = (Bz + bst*w*(dvelz-betaz/alp))/w

  !b_j
  bslowx = betaxlow*bst/alp + gxx*bsx + gxy*bsy + gxz*bsz
  bslowy = betaylow*bst/alp + gxy*bsx + gyy*bsy + gyz*bsz
  bslowz = betazlow*bst/alp + gxz*bsx + gyz*bsy + gzz*bsz

  !b^mu b_mu
  bs2 = (Bx*Blowx + By*Blowy + Bz*Blowz + bst**2)/w**2

  ddens = sqrt(det) * drho * w
  dsx   = sqrt(det) * ((drho*(1+deps)+dpress+bs2)*w*w * vlowx - bst*bslowx)
  dsy   = sqrt(det) * ((drho*(1+deps)+dpress+bs2)*w*w * vlowy - bst*bslowy)
  dsz   = sqrt(det) * ((drho*(1+deps)+dpress+bs2)*w*w * vlowz - bst*bslowz)
  dtau  = sqrt(det) * ((drho*(1+deps)+dpress+bs2)*w*w - (dpress+0.5d0*bs2) - &
       bst**2) - ddens
  dBx = sqrt(det)*Bx
  dBy = sqrt(det)*By
  dBz = sqrt(det)*Bz


end subroutine prim2conpolytype_mhd



 /*@@
   @routine    Primitive2ConservativePolyCells
   @date       Sun Mar 10 21:16:20 2002
   @author
   @desc
   Wrapper function that converts primitive to conservative at the
     cell centres.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/


subroutine Primitive2ConservativePolyCells(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: i, j, k, mhd_stc
  CCTK_REAL :: det

  !we need a larger stencil, actually because of Constrained transport and high-order methods
  if (high_order .ne. 0) then
  		if (approximation_order == 4) mhd_stc = 1
  		if (approximation_order == 6) mhd_stc = 2
  else
  		mhd_stc = 1
  end if


  do k = Spritz_stencil-mhd_stc,cctk_lsh(3)-Spritz_stencil+1+mhd_stc
    do j = Spritz_stencil-mhd_stc,cctk_lsh(2)-Spritz_stencil+1+mhd_stc
      do i = Spritz_stencil-mhd_stc,cctk_lsh(1)-Spritz_stencil+1+mhd_stc

        call SpatialDeterminant(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),&
             gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det)

        call prim2conpolytype_mhd(Spritz_eos_handle,gxx(i,j,k),gxy(i,j,k),&
             gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
             det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
             tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
             eps(i,j,k),press(i,j,k),w_lorentz(i,j,k),&
             Bx(i,j,k),By(i,j,k),Bz(i,j,k),&
             dBx(i,j,k),dBy(i,j,k),dBz(i,j,k),alp(i,j,k),&
             betax(i,j,k),betay(i,j,k),betaz(i,j,k))
      end do
    end do
  end do

  if(evolve_Y_e.ne.0) then
     !$OMP PARALLEL DO PRIVATE(i, j, k)
     do k = 1,cctk_lsh(3)
        do j = 1,cctk_lsh(2)
           do i = 1,cctk_lsh(1)
              Y_e_con(i,j,k) = Y_e(i,j,k) * dens(i,j,k)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif

end subroutine Primitive2ConservativePolyCells
