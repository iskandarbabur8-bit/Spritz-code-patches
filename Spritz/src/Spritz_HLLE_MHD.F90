  /*@@
   @file      Spritz_HLLE_MHD.F90
   @date
   @author    Bruno Giacomazzo
   @desc
   The HLLE solver for GRMHD. Called from the wrapper function, so works in
   all directions.
   @enddesc
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#define velx(i,j,k) vel(i,j,k,1)
#define vely(i,j,k) vel(i,j,k,2)
#define velz(i,j,k) vel(i,j,k,3)

! second order centered finite differencing
#define DIFF_X(VAR,i,j,k)                                                     \
  0.5d0/CCTK_DELTA_SPACE(1)*(VAR(i+1,j,k)-VAR(i-1,j,k))
#define DIFF_Y(VAR,i,j,k)                                                     \
  0.5d0/CCTK_DELTA_SPACE(2)*(VAR(i,j+1,k)-VAR(i,j-1,k))
#define DIFF_Z(VAR,i,j,k)                                                     \
  0.5d0/CCTK_DELTA_SPACE(3)*(VAR(i,j,k+1)-VAR(i,j,k-1))


 /*@@
   @routine    Spritz_HLLE_MHD
   @date
   @author     Bruno Giacomazzo, Lorenzo Sala, Beatrice Giudici
   @desc
   The HLLE solver for GRMHD with a CT scheme.
   @enddesc
   @calls
   @calledby
@@*/

subroutine Spritz_HLLE_MHD(CCTK_ARGUMENTS)
  use Interfaces, only:num_x_flux_MHD, eigenvalues_MHD
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: i, j, k, m
  CCTK_REAL, dimension(8) :: consp,consm_i,fplus,fminus
  CCTK_REAL, dimension(8) :: f1
  CCTK_REAL, dimension(4) :: lamplus,lamminus
  CCTK_REAL, dimension(8) :: qdiff
  CCTK_REAL ::  charmin, charmax, charpm,avg_alp,avg_det
  CCTK_REAL :: gxxh, gxyh, gxzh, gyyh, gyzh, gzzh, uxxh, uxyh, &
       uxzh, uyyh, uyzh, uzzh, avg_betax,avg_betay,avg_betaz, usendh, psi4h
  CCTK_REAL :: Bvminus, Bvplus, b0minus, b0plus, bs2minus, bs2plus, bsxlowplus, bsxlowminus, &
       bsylowminus, bsylowplus, bszlowminus, bszlowplus, bsxminus, bsxplus, &
       bsyminus, bsyplus, bszminus, bszplus, B2minus, B2plus, dpdrho, dpdeps, cs2minus, cs2plus
  CCTK_REAL :: presstmp, hminus, hplus, ptotminus, ptotplus, epstmp
  CCTK_REAL :: betax_low, betay_low, betaz_low
  CCTK_INT:: ierr
  CCTK_INT :: HO_stc = 0

!#include "EOS_Base.inc"

  !!!!!!
  ! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1),xcs2(1)
  CCTK_REAL :: xdpdeps(1),xdpdrho(1)

  n = 1;keytemp = 0;anyerr = 0;keyerr = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
  ! end EOS Omni vars
  !!!!!

!!$  write(*,*) 'Spritz_HLLE_MHD'
!!$  write(*,*) 'flux=',flux_direction

!! If high order approximation is enabled, we need a larger stencil
	if (high_order .ne. 0) then
!!$		write(*,*) "approx_order: ", approximation_order
		if (approximation_order == 4) HO_stc = 1
		if (approximation_order == 6) HO_stc = 2
	end if


  !$OMP PARALLEL DO PRIVATE(i, j, k, m, f1, lamminus, lamplus, &
  !$OMP consp, consm_i, fplus, fminus, qdiff, charmin, charmax, &
  !$OMP charpm, gxxh, gxyh, gxzh, gyyh, gyzh, gzzh, uxxh, uxyh, &
  !$OMP uxzh, uyyh, uyzh, uzzh, avg_betax,avg_betay,avg_betaz, &
  !$OMP usendh, psi4h, Bvminus, Bvplus, b0minus, b0plus, bs2minus, &
  !$OMP bs2plus, bsxlowplus, bsxlowminus, bsylowminus, bsylowplus, &
  !$OMP bszlowminus, bszlowplus, bsxminus, bsxplus, bsyminus, bsyplus, &
  !$OMP bszminus, bszplus, B2minus, B2plus, dpdrho, dpdeps, cs2minus, &
  !$OMP cs2plus, presstmp, hminus, hplus, ptotminus, ptotplus, &
  !$OMP betax_low, betay_low, betaz_low, avg_alp, avg_det, ierr,&
  !$OMP xrho, xeps, xpress, xtemp, xye, xdpdeps, xdpdrho, anyerr, &
  !$OMP keyerr)
  do k = Spritz_stencil - HO_stc, cctk_lsh(3) - Spritz_stencil + 1 + HO_stc
     do j = Spritz_stencil - HO_stc, cctk_lsh(2) - Spritz_stencil + 1 + HO_stc
        do i = Spritz_stencil - HO_stc, cctk_lsh(1) - Spritz_stencil + 1 + HO_stc

           f1 = 0.d0
           lamminus = 0.d0
           lamplus = 0.d0
           consp = 0.d0
           consm_i = 0.d0
           fplus = 0.d0
           fminus = 0.d0
           qdiff = 0.d0

!!$        Set the left (p for plus) and right (m_i for minus, i+1) states

           consp(1)   = densplus(i,j,k)
           consp(2)   = sxplus(i,j,k)
           consp(3)   = syplus(i,j,k)
           consp(4)   = szplus(i,j,k)
           consp(5)   = tauplus(i,j,k)
           consp(6)   = Bxplus(i,j,k)
           consp(7)   = Byplus(i,j,k)
           consp(8)   = Bzplus(i,j,k)

           consm_i(1) = densminus(i+xoffset,j+yoffset,k+zoffset)
           consm_i(2) = sxminus(i+xoffset,j+yoffset,k+zoffset)
           consm_i(3) = syminus(i+xoffset,j+yoffset,k+zoffset)
           consm_i(4) = szminus(i+xoffset,j+yoffset,k+zoffset)
           consm_i(5) = tauminus(i+xoffset,j+yoffset,k+zoffset)
           consm_i(6) = Bxminus(i+xoffset,j+yoffset,k+zoffset)
           consm_i(7) = Byminus(i+xoffset,j+yoffset,k+zoffset)
           consm_i(8) = Bzminus(i+xoffset,j+yoffset,k+zoffset)


!!$        Calculate various metric terms here.
!!$        Note also need the average of the lapse at the
!!$        left and right points.



              if (shift_state .ne. 0) then
                 avg_betax = 0.5d0 * (betax(i+xoffset,j+yoffset,k+zoffset) + &
                      betax(i,j,k))
                 avg_betay = 0.5d0 * (betay(i+xoffset,j+yoffset,k+zoffset) + &
                      betay(i,j,k))
                 avg_betaz = 0.5d0 * (betaz(i+xoffset,j+yoffset,k+zoffset) + &
                      betaz(i,j,k))
              else
                 avg_betax = 0.d0
                 avg_betay = 0.d0
                 avg_betaz = 0.d0
              end if

              avg_alp = 0.5 * (alp(i,j,k) + alp(i+xoffset,j+yoffset,k+zoffset))

              gxxh = 0.5d0 * (gxx(i+xoffset,j+yoffset,k+zoffset) + gxx(i,j,k))
              gxyh = 0.5d0 * (gxy(i+xoffset,j+yoffset,k+zoffset) + gxy(i,j,k))
              gxzh = 0.5d0 * (gxz(i+xoffset,j+yoffset,k+zoffset) + gxz(i,j,k))
              gyyh = 0.5d0 * (gyy(i+xoffset,j+yoffset,k+zoffset) + gyy(i,j,k))
              gyzh = 0.5d0 * (gyz(i+xoffset,j+yoffset,k+zoffset) + gyz(i,j,k))
              gzzh = 0.5d0 * (gzz(i+xoffset,j+yoffset,k+zoffset) + gzz(i,j,k))

              if (conformal_state .eq. 0) then
                 psi4h = 1.d0
              else
                 psi4h = (0.5d0*(psi(i,j,k)+psi(i+xoffset,j+yoffset,k+zoffset)))**4
              end if



           call SpatialDeterminant(psi4h*gxxh,psi4h*gxyh,psi4h*gxzh,&
                psi4h*gyyh,psi4h*gyzh,psi4h*gzzh,avg_det)

           !B^i v_i
           Bvminus = gxxh*Bxminus(i+xoffset,j+yoffset,k+zoffset)*velxminus(i+xoffset,j+yoffset,k+zoffset) + &
                gyyh*Byminus(i+xoffset,j+yoffset,k+zoffset)*velyminus(i+xoffset,j+yoffset,k+zoffset) + &
                gzzh*Bzminus(i+xoffset,j+yoffset,k+zoffset)*velzminus(i+xoffset,j+yoffset,k+zoffset) + &
                gxyh*(Bxminus(i+xoffset,j+yoffset,k+zoffset)*velyminus(i+xoffset,j+yoffset,k+zoffset) + &
                Byminus(i+xoffset,j+yoffset,k+zoffset)*velxminus(i+xoffset,j+yoffset,k+zoffset)) + &
                gxzh*(Bxminus(i+xoffset,j+yoffset,k+zoffset)*velzminus(i+xoffset,j+yoffset,k+zoffset) + &
                Bzminus(i+xoffset,j+yoffset,k+zoffset)*velxminus(i+xoffset,j+yoffset,k+zoffset)) + &
                gyzh*(Byminus(i+xoffset,j+yoffset,k+zoffset)*velzminus(i+xoffset,j+yoffset,k+zoffset) + &
                Bzminus(i+xoffset,j+yoffset,k+zoffset)*velyminus(i+xoffset,j+yoffset,k+zoffset))


           Bvplus = gxxh*Bxplus(i,j,k)*velxplus(i,j,k) + gyyh*Byplus(i,j,k)*velyplus(i,j,k) + &
                gzzh*Bzplus(i,j,k)*velzplus(i,j,k) + gxyh*(Bxplus(i,j,k)*velyplus(i,j,k) + &
                Byplus(i,j,k)*velxplus(i,j,k)) + gxzh*(Bxplus(i,j,k)*velzplus(i,j,k) + &
                Bzplus(i,j,k)*velxplus(i,j,k)) +  gyzh*(Byplus(i,j,k)*velzplus(i,j,k) + &
                Bzplus(i,j,k)*velyplus(i,j,k))



           b0minus = w_lorentzminus(i+xoffset,j+yoffset,k+zoffset)*Bvminus/avg_alp
           b0plus  = w_lorentzplus(i,j,k)*Bvplus/avg_alp

           !b^x
           bsxminus = Bxminus(i+xoffset,j+yoffset,k+zoffset)/w_lorentzminus(i+xoffset,j+yoffset,k+zoffset) + &
                Bvminus*w_lorentzminus(i+xoffset,j+yoffset,k+zoffset)*&
                (velxminus(i+xoffset,j+yoffset,k+zoffset)-avg_betax/avg_alp)
           bsxplus = Bxplus(i,j,k)/w_lorentzplus(i,j,k) + &
                Bvplus*w_lorentzplus(i,j,k)*(velxplus(i,j,k)-avg_betax/avg_alp)

           !b^y
           bsyminus = Byminus(i+xoffset,j+yoffset,k+zoffset)/w_lorentzminus(i+xoffset,j+yoffset,k+zoffset) + &
                Bvminus*w_lorentzminus(i+xoffset,j+yoffset,k+zoffset)*&
                (velyminus(i+xoffset,j+yoffset,k+zoffset)-avg_betay/avg_alp)
           bsyplus = Byplus(i,j,k)/w_lorentzplus(i,j,k) + &
                Bvplus*w_lorentzplus(i,j,k)*(velyplus(i,j,k)-avg_betay/avg_alp)

           !b^z
           bszminus = Bzminus(i+xoffset,j+yoffset,k+zoffset)/w_lorentzminus(i+xoffset,j+yoffset,k+zoffset) + &
                Bvminus*w_lorentzminus(i+xoffset,j+yoffset,k+zoffset)*&
                (velzminus(i+xoffset,j+yoffset,k+zoffset)-avg_betaz/avg_alp)
           bszplus = Bzplus(i,j,k)/w_lorentzplus(i,j,k) + &
                Bvplus*w_lorentzplus(i,j,k)*(velzplus(i,j,k)-avg_betaz/avg_alp)

           B2minus = gxxh*Bxminus(i+xoffset,j+yoffset,k+zoffset)**2 + &
                gyyh*Byminus(i+xoffset,j+yoffset,k+zoffset)**2 + &
                gzzh*Bzminus(i+xoffset,j+yoffset,k+zoffset)**2 + &
                2.0*(gxyh*Bxminus(i+xoffset,j+yoffset,k+zoffset)*Byminus(i+xoffset,j+yoffset,k+zoffset) + &
                gxzh*Bxminus(i+xoffset,j+yoffset,k+zoffset)*Bzminus(i+xoffset,j+yoffset,k+zoffset) + &
                gyzh*Byminus(i+xoffset,j+yoffset,k+zoffset)*Bzminus(i+xoffset,j+yoffset,k+zoffset))

           B2plus = gxxh*Bxplus(i,j,k)**2 + gyyh*Byplus(i,j,k)**2 + &
                gzzh*Bzplus(i,j,k)**2 + &
                2.0*(gxyh*Bxplus(i,j,k)*Byplus(i,j,k) + &
                gxzh*Bxplus(i,j,k)*Bzplus(i,j,k) + &
                gyzh*Byplus(i,j,k)*Bzplus(i,j,k))


           bs2minus = B2minus/w_lorentzminus(i+xoffset,j+yoffset,k+zoffset)**2 + &
                Bvminus**2
           bs2plus  = B2plus/w_lorentzplus(i,j,k)**2 + &
                Bvplus**2

           betax_low = gxxh*avg_betax + gxyh*avg_betay + gxzh*avg_betaz
           betay_low = gxyh*avg_betax + gyyh*avg_betay + gyzh*avg_betaz
           betaz_low = gxzh*avg_betax + gyzh*avg_betay + gzzh*avg_betaz


           !b_x
           bsxlowplus  = betax_low*b0plus  + gxxh*bsxplus  + gxyh*bsyplus  + gxzh*bszplus
           bsxlowminus = betax_low*b0minus + gxxh*bsxminus + gxyh*bsyminus + gxzh*bszminus
           !b_y
           bsylowplus  = betay_low*b0plus  + gxyh*bsxplus  + gyyh*bsyplus  + gyzh*bszplus
           bsylowminus = betay_low*b0minus + gxyh*bsxminus + gyyh*bsyminus + gyzh*bszminus
           !b_z
           bszlowplus  = betaz_low*b0plus  + gxzh*bsxplus  + gyzh*bsyplus  + gzzh*bszplus
           bszlowminus = betaz_low*b0minus + gxzh*bsxminus + gyzh*bsyminus + gzzh*bszminus

           !LSADD: this is the standard routine that finds h and cs2, without using the tabulated one.
             ! xtemp =temperature(i,j,k) !this is just a guess, or, if constraining, the actual value
             ! xrho = rhominus(i+xoffset,j+yoffset,k+zoffset)
             ! xeps = epsminus(i+xoffset,j+yoffset,k+zoffset)
             ! if(evolve_Y_e.eq.1) then
             !   xye  = y_e_minus(i+xoffset,j+yoffset,k+zoffset)
             ! end if
             !
             ! keytemp = 0
             ! if(Constrain_Temp_now.eq.1) keytemp = 1
             ! call EOS_Omni_press(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
             !       xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
             ! presstmp = xpress(1)
             ! if(CCTK_EQUALS(Spritz_eos_table,"nuc_eos")) keytemp = 1 !if nuc_eos, there's no need in finding T again
             ! call EOS_Omni_DPressByDRho(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
  			     !     xrho,xeps,xtemp,xye,xdpdrho,keyerr,anyerr)
             ! dpdrho = xdpdrho(1)
             ! call EOS_Omni_DPressByDEps(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
  			     !     xrho,xeps,xtemp,xye,xdpdeps,keyerr,anyerr)
             ! dpdeps = xdpdeps(1)
             ! !!!!!!
             ! hminus = 1.0d0 + epsminus(i+xoffset,j+yoffset,k+zoffset) + presstmp/rhominus(i+xoffset,j+yoffset,k+zoffset)
             ! cs2minus = (dpdrho + presstmp * dpdeps / (rhominus(i+xoffset,j+yoffset,k+zoffset)**2))/hminus
             ! ptotminus = presstmp + Spritz_mhd_coeff*0.5d0*bs2minus
             !
             ! xrho = rhoplus(i,j,k)
             ! xeps = epsplus(i,j,k)
             ! if(evolve_Y_e.eq.1) then
             !   xye  = y_e_plus(i,j,k)
             ! end if
             !
             ! keytemp = 0
             ! if(Constrain_Temp_now.eq.1) keytemp = 1
             ! call EOS_Omni_press(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
             !       xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
             ! presstmp = xpress(1)
             ! if(CCTK_EQUALS(Spritz_eos_table,"nuc_eos")) keytemp = 1 !if nuc_eos, there's no need in finding T again
             ! call EOS_Omni_DPressByDRho(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
  			     !     xrho,xeps,xtemp,xye,xdpdrho,keyerr,anyerr)
             ! dpdrho = xdpdrho(1)
             ! call EOS_Omni_DPressByDEps(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
  			     !     xrho,xeps,xtemp,xye,xdpdeps,keyerr,anyerr)
             ! dpdeps = xdpdeps(1)
             ! keytemp = 0
             ! !!!!!!
             ! hplus = 1.0d0 + epsplus(i,j,k) + presstmp/rhoplus(i,j,k)
             ! cs2plus = (dpdrho + presstmp * dpdeps / (rhoplus(i,j,k)**2))/hplus
             ! ptotplus = presstmp + Spritz_mhd_coeff*0.5d0*bs2plus

             if(CCTK_EQUALS(Spritz_eos_table,"nuc_eos")) then
               xtemp =temperature(i,j,k) !this is just a guess, or, if constraining, the actual value
               xrho = rhominus(i+xoffset,j+yoffset,k+zoffset)
               xeps = epsminus(i+xoffset,j+yoffset,k+zoffset)
               if(evolve_Y_e.eq.1) then
                 xye  = y_e_minus(i+xoffset,j+yoffset,k+zoffset)
               end if

               keytemp = 0
               if(Constrain_Temp_now.eq.1) keytemp = 1
               call EOS_Omni_press(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                     xrho,xeps,xtemp,xye,xpress,keyerr,anyerr) !find temperature, pressure and, eventually, energy
               presstmp = xpress(1)
               epstmp   = xeps(1) !we don't care about temperature, it is stored in xtemp
               keytemp  = 1 !now we perfectly know temperature, let's find cs2
               call EOS_Omni_cs2(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                  xrho,xeps,xtemp,xye,xcs2,keyerr,anyerr)
               keytemp = 0
               !!!!!!
               hminus = 1.0d0 + epstmp + presstmp/rhominus(i+xoffset,j+yoffset,k+zoffset)
               cs2minus = xcs2(1)
               ptotminus = presstmp + Spritz_mhd_coeff*0.5d0*bs2minus

               xrho = rhoplus(i,j,k)
               xeps = epsplus(i,j,k)
               if(evolve_Y_e.eq.1) then
                 xye  = y_e_plus(i,j,k)
               end if

               keytemp = 0
               if(Constrain_Temp_now.eq.1) keytemp = 1
               call EOS_Omni_press(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                     xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
               presstmp = xpress(1)
               epstmp   = xeps(1) !we don't care about temperature, it is stored in xtemp
               keytemp  = 1 !now we perfectly know temperature, let's find cs2
               call EOS_Omni_cs2(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                  xrho,xeps,xtemp,xye,xcs2,keyerr,anyerr)
               keytemp = 0
               !!!!!!
               hplus = 1.0d0 + epstmp + presstmp/rhoplus(i,j,k)
               cs2plus = xcs2(1)
               ptotplus = presstmp + Spritz_mhd_coeff*0.5d0*bs2plus
             else !non-nuclear EoS
               xrho = rhominus(i+xoffset,j+yoffset,k+zoffset)
               xeps = epsminus(i+xoffset,j+yoffset,k+zoffset)

               keytemp = 0
               call EOS_Omni_press(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                     xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
               presstmp = xpress(1)
               call EOS_Omni_DPressByDRho(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                  xrho,xeps,xtemp,xye,xdpdrho,keyerr,anyerr)
               dpdrho = xdpdrho(1)
               call EOS_Omni_DPressByDEps(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                  xrho,xeps,xtemp,xye,xdpdeps,keyerr,anyerr)
               dpdeps = xdpdeps(1)
               !!!!!!
               hminus = 1.0d0 + epsminus(i+xoffset,j+yoffset,k+zoffset) + presstmp/rhominus(i+xoffset,j+yoffset,k+zoffset)
               cs2minus = (dpdrho + presstmp * dpdeps / (rhominus(i+xoffset,j+yoffset,k+zoffset)**2))/hminus
               ptotminus = presstmp + Spritz_mhd_coeff*0.5d0*bs2minus

               xrho = rhoplus(i,j,k)
               xeps = epsplus(i,j,k)

               keytemp = 0
               call EOS_Omni_press(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                     xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
               presstmp = xpress(1)
               call EOS_Omni_DPressByDRho(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                  xrho,xeps,xtemp,xye,xdpdrho,keyerr,anyerr)
               dpdrho = xdpdrho(1)
               call EOS_Omni_DPressByDEps(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                  xrho,xeps,xtemp,xye,xdpdeps,keyerr,anyerr)
               dpdeps = xdpdeps(1)
               keytemp = 0
               !!!!!!
               hplus = 1.0d0 + epsplus(i,j,k) + presstmp/rhoplus(i,j,k)
               cs2plus = (dpdrho + presstmp * dpdeps / (rhoplus(i,j,k)**2))/hplus
               ptotplus = presstmp + Spritz_mhd_coeff*0.5d0*bs2plus
             end if

              call UpperMetric(uxxh, uxyh, uxzh, uyyh, uyzh, uzzh, &
                   avg_det,psi4h*gxxh, psi4h*gxyh, psi4h*gxzh, &
                   psi4h*gyyh, psi4h*gyzh, psi4h*gzzh)

!!$        Eigenvalues and fluxes either side of the cell interface


!!$              write(*,*) 'i=',i,'j=',j,'k=',k, 'flux_direction=',flux_direction

              if (flux_direction == 1) then              
                 call eigenvalues_MHD(avg_alp,bs2minus,b0minus,cs2minus,&
                      velxminus(i+xoffset,j+yoffset,k+zoffset),&
                      Bxminus(i+xoffset,j+yoffset,k+zoffset),bsxminus,avg_betax,&
                      uxxh,rhominus(i+xoffset,j+yoffset,k+zoffset),hminus,&
                      w_lorentzminus(i+xoffset,j+yoffset,k+zoffset),lamminus,x(i,j,k),y(i,j,k),z(i,j,k))
                 call eigenvalues_MHD(avg_alp,bs2plus,b0plus,cs2plus,&
                      velxplus(i,j,k),Bxplus(i,j,k),bsxplus,avg_betax,&
                      uxxh,rhoplus(i,j,k),hplus,&
                      w_lorentzplus(i,j,k),lamplus,x(i,j,k),y(i,j,k),z(i,j,k))

                 call num_x_flux_MHD(consp(1),consp(2),consp(3),consp(4),consp(5),&
                      consp(6),consp(7),consp(8),fplus(1),fplus(2),fplus(3),fplus(4),fplus(5),&
                      fplus(6),fplus(7),fplus(8),velxplus(i,j,k),velyplus(i,j,k),velzplus(i,j,k),&
                      avg_betax,avg_betay,avg_betaz,avg_alp,avg_det,w_lorentzplus(i,j,k),ptotplus,&
                      b0plus,bsxlowplus,bsylowplus,bszlowplus)                     

                 call num_x_flux_MHD(consm_i(1),consm_i(2),consm_i(3),consm_i(4),consm_i(5),&
                      consm_i(6),consm_i(7),consm_i(8),fminus(1),fminus(2),fminus(3),fminus(4),fminus(5),&
                      fminus(6),fminus(7),fminus(8),velxminus(i+xoffset,j+yoffset,k+zoffset),&
                      velyminus(i+xoffset,j+yoffset,k+zoffset),velzminus(i+xoffset,j+yoffset,k+zoffset),&
                      avg_betax,avg_betay,avg_betaz,avg_alp,avg_det,&
                      w_lorentzminus(i+xoffset,j+yoffset,k+zoffset),ptotminus,&
                      b0minus,bsxlowminus,bsylowminus,bszlowminus)                      
              	  
              else if (flux_direction == 2) then
                 call eigenvalues_MHD(avg_alp,bs2minus,b0minus,cs2minus,&
                      velyminus(i+xoffset,j+yoffset,k+zoffset),&
                      Byminus(i+xoffset,j+yoffset,k+zoffset),bsyminus,avg_betay,&
                      uyyh,rhominus(i+xoffset,j+yoffset,k+zoffset),hminus,&
                      w_lorentzminus(i+xoffset,j+yoffset,k+zoffset),lamminus,x(i,j,k),y(i,j,k),z(i,j,k))
                 call eigenvalues_MHD(avg_alp,bs2plus,b0plus,cs2plus,&
                      velyplus(i,j,k),Byplus(i,j,k),bsyplus,avg_betay,&
                      uyyh,rhoplus(i,j,k),hplus,&
                      w_lorentzplus(i,j,k),lamplus,x(i,j,k),y(i,j,k),z(i,j,k))

                 call num_x_flux_MHD(consp(1),consp(3),consp(2),consp(4),consp(5),&
                      consp(7),consp(6),consp(8),fplus(1),fplus(3),fplus(2),fplus(4),fplus(5),&
                      fplus(7),fplus(6),fplus(8),velyplus(i,j,k),velxplus(i,j,k),velzplus(i,j,k),&
                      avg_betay,avg_betax,avg_betaz,avg_alp,avg_det,w_lorentzplus(i,j,k),ptotplus,&
                      b0plus,bsylowplus,bsxlowplus,bszlowplus)

                 call num_x_flux_MHD(consm_i(1),consm_i(3),consm_i(2),consm_i(4),consm_i(5),&
                      consm_i(7),consm_i(6),consm_i(8),fminus(1),fminus(3),fminus(2),fminus(4),fminus(5),&
                      fminus(7),fminus(6),fminus(8),velyminus(i+xoffset,j+yoffset,k+zoffset),&
                      velxminus(i+xoffset,j+yoffset,k+zoffset),velzminus(i+xoffset,j+yoffset,k+zoffset),&
                      avg_betay,avg_betax,avg_betaz,avg_alp,avg_det, &
                      w_lorentzminus(i+xoffset,j+yoffset,k+zoffset),ptotminus,&
                      b0minus,bsylowminus,bsxlowminus,bszlowminus)
              
              else if (flux_direction == 3) then

                 call eigenvalues_MHD(avg_alp,bs2minus,b0minus,cs2minus,&
                      velzminus(i+xoffset,j+yoffset,k+zoffset),&
                      Bzminus(i+xoffset,j+yoffset,k+zoffset),bszminus,avg_betaz,&
                      uzzh,rhominus(i+xoffset,j+yoffset,k+zoffset),hminus,&
                      w_lorentzminus(i+xoffset,j+yoffset,k+zoffset),lamminus,x(i,j,k),y(i,j,k),z(i,j,k))


                 call eigenvalues_MHD(avg_alp,bs2plus,b0plus,cs2plus,&
                      velzplus(i,j,k),Bzplus(i,j,k),bszplus,avg_betaz,&
                      uzzh,rhoplus(i,j,k),hplus,&
                      w_lorentzplus(i,j,k),lamplus,x(i,j,k),y(i,j,k),z(i,j,k))

                 call num_x_flux_MHD(consp(1),consp(4),consp(3),consp(2),consp(5),&
                      consp(8),consp(7),consp(6),fplus(1),fplus(4),fplus(3),fplus(2),fplus(5),&
                      fplus(8),fplus(7),fplus(6),velzplus(i,j,k),velyplus(i,j,k),velxplus(i,j,k),&
                      avg_betaz,avg_betay,avg_betax,avg_alp,avg_det,w_lorentzplus(i,j,k),ptotplus,&
                      b0plus,bszlowplus,bsylowplus,bsxlowplus)

                 call num_x_flux_MHD(consm_i(1),consm_i(4),consm_i(3),consm_i(2),consm_i(5),&
                      consm_i(8),consm_i(7),consm_i(6),fminus(1),fminus(4),fminus(3),fminus(2),fminus(5),&
                      fminus(8),fminus(7),fminus(6),velzminus(i+xoffset,j+yoffset,k+zoffset),&
                      velyminus(i+xoffset,j+yoffset,k+zoffset),velxminus(i+xoffset,j+yoffset,k+zoffset),&
                      avg_betaz,avg_betay,avg_betax,avg_alp,avg_det,&
                      w_lorentzminus(i+xoffset,j+yoffset,k+zoffset),ptotminus,&
                      b0minus,bszlowminus,bsylowminus,bsxlowminus)
              
              else
                 call CCTK_Warn(0, "Flux direction not x,y,z")
              end if

              if (CCTK_EQUALS(flux_type,"HLLE")) then
!!$              Find minimum and maximum wavespeeds
                 charmin = min(0.d0, lamplus(1), lamplus(2), lamplus(3), &
                      lamplus(4),  lamminus(1),lamminus(2),lamminus(3),&
                      lamminus(4))

                 charmax = max(0.d0, lamplus(1), lamplus(2), lamplus(3), &
                      lamplus(4),  lamminus(1),lamminus(2),lamminus(3),&
                      lamminus(4))

                 charpm = charmax - charmin

!!$              Calculate flux by standard formula

                 do m = 1,8

                    qdiff(m) = consm_i(m) - consp(m)

                    !consm_i and consp(6:8) are Bx,By,Bz (not dBx,dBy,dBz)
                    if(m>5) qdiff(m)=sqrt(avg_det)*qdiff(m)

                    f1(m) = (charmax * fplus(m) - charmin * fminus(m) + &
                         charmax * charmin * qdiff(m)) / charpm

                 end do

              else if (CCTK_EQUALS(flux_type,"LxF")) then

!!$              Find minimum and maximum wavespeeds for Lax-Friedrichs formula
                 charmax = max( abs(0.d0), abs(lamplus(1)), abs(lamplus(2)), abs(lamplus(3)), &
                                abs(lamplus(4)),  abs(lamminus(1)),abs(lamminus(2)),abs(lamminus(3)),&
                                abs(lamminus(4)) )
                 charmin = -1.0d0 * charmax

                 charpm = 2.0d0

!!$              Calculate flux by Lax-Friedrichs formula, avoiding problems if charmax = 0.0d0

                 do m = 1,8

                    qdiff(m) = consm_i(m) - consp(m)

                    !consm_i and consp(6:8) are Bx,By,Bz (not dBx,dBy,dBz)
                    if(m>5) qdiff(m)=sqrt(avg_det)*qdiff(m)

                    f1(m) = (fplus(m) + fminus(m) - &
                         charmax * qdiff(m)) / charpm

                 end do

              else
                 call CCTK_Warn(0, "flux_type not recognized")
              end if

		     densflux(i,j,k) = f1(1)
		     sxflux(i,j,k)   = f1(2)
		     syflux(i,j,k)   = f1(3)
		     szflux(i,j,k)   = f1(4)
		     tauflux(i,j,k)  = f1(5)			  
           Bxflux(i,j,k)   = f1(6)*(1-xoffset)!force the x flux of Bx to be zero
           Byflux(i,j,k)   = f1(7)*(1-yoffset)!force the y flux of By to be zero
           Bzflux(i,j,k)   = f1(8)*(1-zoffset)!force the z flux of Bz to be zero

!!$           densflux(i,j,k) = 0.0d0! f1(1)
!!$           sxflux(i,j,k)   = 0.0d0!f1(2)
!!$           syflux(i,j,k)   = 0.0d0!f1(3)
!!$           szflux(i,j,k)   = 0.0d0!f1(4)
!!$           tauflux(i,j,k)  = 0.0d0!f1(5)
!!$           Bxflux(i,j,k)   = 0.0d0!f1(6)*(1-xoffset)!force the x flux of Bx to be zero
!!$           Byflux(i,j,k)   = 0.0d0!f1(7)*(1-yoffset)!force the y flux of By to be zero
!!$           Bzflux(i,j,k)   = 0.0d0!f1(8)*(1-zoffset)!force the z flux of Bz to be zero
!!$           if ((x(i,j,k)==0.0d0).and.(y(i,j,k)==0.0d0).and.(z(i,j,k)==0.0d0)) call CCTK_WARN(1,"Setting the fluxes to zero")

           Bxflux(i,j,k)   = avg_alp*Bxflux(i,j,k)
           Byflux(i,j,k)   = avg_alp*Byflux(i,j,k)
           Bzflux(i,j,k)   = avg_alp*Bzflux(i,j,k)
	
			  if(evolve_Y_e.eq.1) then
				 if (densflux(i, j, k) > 0.d0) then
					 Y_e_con_flux(i, j, k) = &
							Y_e_plus(i, j, k) * &
							densflux(i, j, k)
				 else
					 Y_e_con_flux(i, j, k) = &
						  Y_e_minus(i + xoffset, j + yoffset, k + zoffset) * &
							densflux(i, j, k)
				 endif
!!$ Computing the electron fraction using directly the Riemann solver
!				 if (flux_direction == 1) then
!				 	veltilde_p = (velxplus(i,j,k) - avg_betax/avg_alp)
!				 	veltilde_m = (velxminus(i+xoffset, j+yoffset, k+zoffset) - avg_betax/avg_alp)
!				 else if (flux_direction == 2) then
!				 	veltilde_p = (velyplus(i,j,k) - avg_betay/avg_alp)
!				 	veltilde_m = (velyminus(i+xoffset, j+yoffset, k+zoffset) - avg_betay/avg_alp)			 
!				 else if (flux_direction == 3) then
!				 	veltilde_p = (velzplus(i,j,k) - avg_betaz/avg_alp)
!				 	veltilde_m = (velzminus(i+xoffset, j+yoffset, k+zoffset) - avg_betaz/avg_alp)
!				 else
!				 	call CCTK_Warn(0, "Flux direction not x,y,z")
!				 end if

!				 FYe0plus = densplus(i,j,k) * Y_e_plus(i,j,k)
!				 FYe0minus = densminus(i+xoffset, j+yoffset, k+zoffset) * &
!				           Y_e_minus(i+xoffset, j+yoffset, k+zoffset)
	 
!				 FYei_plus = veltilde_p * F0plus
!				 FYei_minus = veltilde_m * F0minus
		 
!				 if (CCTK_EQUALS(flux_type,"HLLE")) then
!				 	Y_e_con_flux(i,j,k) = (charmax * Fi_plus - charmin * Fi_minus + &
!				 									charmax * charmin * (F0minus - F0plus)) / charpm			 
!				 else if (CCTK_EQUALS(flux_type,"LxF") then
!				 	Y_e_con_flux(i,j,k) = (Fi_plus + Fi_minus - &
!				 									charmax*(F0minus - F0plus)) / charpm
!				 else
!				 	call CCTK_Warn(0, "flux_type not recognized")
!				 end if
				 	
  			  endif
	        
!!$ If high order approximation is enabled, we must keep track of the fluxes before the modification
!!$ so we use different variables for a preliminar storing
			  if (high_order .ne. 0) then
				  dens_HO(i,j,k) = f1(1)*avg_alp
				  sx_HO(i,j,k)   = f1(2)*avg_alp
				  sy_HO(i,j,k)   = f1(3)*avg_alp
				  sz_HO(i,j,k)   = f1(4)*avg_alp
				  tau_HO(i,j,k)  = f1(5)*avg_alp
				  
				  if (evolve_Y_e.eq.1) then
				  	Y_e_HO(i,j,k) = Y_e_con_flux(i,j,k)*avg_alp
				  end if
			  end if
        
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  !$OMP BARRIER


  !Constrained Transport Scheme to keep the divergence of the magnetic field constant
  ! (zero if it was zero at t=0)
  ! Computing the Electric field


  !$OMP PARALLEL DO PRIVATE(i, j)
  do k = Spritz_stencil, cctk_lsh(3) - Spritz_stencil
     do j = Spritz_stencil, cctk_lsh(2) - Spritz_stencil
        do i = Spritz_stencil, cctk_lsh(1) - Spritz_stencil
         
           if(Spritz_Avec_Staggered.ne.0) then
              if (flux_direction == 1) then
                 Ey(i,j,k) = Ey(i,j,k) + 0.25d0*(Bzflux(i,j,k)+Bzflux(i,j,k+1))
                 Ez(i,j,k) = Ez(i,j,k) + 0.25d0*(-Byflux(i,j,k)-Byflux(i,j+1,k))
              else if (flux_direction == 2) then
                 Ex(i,j,k) = Ex(i,j,k) + 0.25d0*(-Bzflux(i,j,k)-Bzflux(i,j,k+1))
                 Ez(i,j,k) = Ez(i,j,k) + 0.25d0*(Bxflux(i,j,k)+Bxflux(i+1,j,k))
              else if (flux_direction == 3) then
                 Ex(i,j,k) = 0.25d0*(Byflux(i,j,k)+Byflux(i,j+1,k))
                 Ey(i,j,k) = 0.25d0*(-Bxflux(i,j,k)-Bxflux(i+1,j,k))
                 Ez(i,j,k) = 0.0d0
              else
                 call CCTK_Warn(0, "Flux direction not x,y,z")
              end if
           else
              if (flux_direction == 1) then
                 Ey(i,j,k) = Ey(i,j,k) + 0.25d0*(Bzflux(i,j,k)+Bzflux(i-1,j,k))
                 Ez(i,j,k) = Ez(i,j,k) + 0.25d0*(-Byflux(i,j,k)-Byflux(i-1,j,k))
              else if (flux_direction == 2) then
                 Ex(i,j,k) = Ex(i,j,k) + 0.25d0*(-Bzflux(i,j,k)-Bzflux(i,j-1,k))
                 Ez(i,j,k) = Ez(i,j,k) + 0.25d0*(Bxflux(i,j,k)+Bxflux(i,j-1,k))
              else if (flux_direction == 3) then
                 Ex(i,j,k) = 0.25d0*(Byflux(i,j,k)+Byflux(i,j,k-1))
                 Ey(i,j,k) = 0.25d0*(-Bxflux(i,j,k)-Bxflux(i,j,k-1))
                 Ez(i,j,k) = 0.0d0
              else
                 call CCTK_Warn(0, "Flux direction not x,y,z")
              end if
           end if

        end do
     end do
  end do
  !$OMP END PARALLEL DO

!! If high_order is enabled, then the fluxes must be corrected in order to achieve
!! high-order convergence 
if (high_order .ne. 0) then

	if (flux_direction == 1) then
		!$OMP PARALLEL DO PRIVATE(i,j,k)
		do k = Spritz_stencil, cctk_lsh(3) - Spritz_stencil + 1
			do j = Spritz_stencil, cctk_lsh(2) - Spritz_stencil + 1
				call HO_correction(densflux(:,j,k), dens_HO(:,j,k), cctk_lsh(1))
				call HO_correction(sxflux(:,j,k), sx_HO(:,j,k), cctk_lsh(1))
				call HO_correction(syflux(:,j,k), sy_HO(:,j,k), cctk_lsh(1))
				call HO_correction(szflux(:,j,k), sz_HO(:,j,k), cctk_lsh(1))
				call HO_correction(tauflux(:,j,k), tau_HO(:,j,k), cctk_lsh(1))
				
				if(evolve_Y_e.eq.1) then
					call HO_correction(Y_e_con_flux(:,j,k), Y_e_HO(:,j,k), cctk_lsh(1))
				end if
			end do
		end do
		!$OMP END PARALLEL DO
	
	else if (flux_direction == 2) then
		!$OMP PARALLEL DO PRIVATE(i,j,k)
		do k = Spritz_stencil, cctk_lsh(3) - Spritz_stencil + 1
			do j = Spritz_stencil, cctk_lsh(1) - Spritz_stencil + 1
				call HO_correction(densflux(j,:,k), dens_HO(j,:,k), cctk_lsh(2))
				call HO_correction(sxflux(j,:,k), sx_HO(j,:,k), cctk_lsh(2))
				call HO_correction(syflux(j,:,k), sy_HO(j,:,k), cctk_lsh(2))
				call HO_correction(szflux(j,:,k), sz_HO(j,:,k), cctk_lsh(2))
				call HO_correction(tauflux(j,:,k), tau_HO(j,:,k), cctk_lsh(2))
				
				if(evolve_Y_e.eq.1) then
					call HO_correction(Y_e_con_flux(j,:,k), Y_e_HO(j,:,k), cctk_lsh(2))
				end if
			end do
		end do
		!$OMP END PARALLEL DO
	
	else if (flux_direction == 3) then
		!$OMP PARALLEL DO PRIVATE(i,j,k)
		do k = Spritz_stencil, cctk_lsh(2) - Spritz_stencil + 1
			do j = Spritz_stencil, cctk_lsh(1) - Spritz_stencil + 1
				call HO_correction(densflux(j,k,:), dens_HO(j,k,:), cctk_lsh(3))
				call HO_correction(sxflux(j,k,:), sx_HO(j,k,:), cctk_lsh(3))
				call HO_correction(syflux(j,k,:), sy_HO(j,k,:), cctk_lsh(3))
				call HO_correction(szflux(j,k,:), sz_HO(j,k,:), cctk_lsh(3))
				call HO_correction(tauflux(j,k,:), tau_HO(j,k,:), cctk_lsh(3))
				
				if(evolve_Y_e.eq.1) then
					call HO_correction(Y_e_con_flux(j,k,:), Y_e_HO(j,k,:), cctk_lsh(3))
				end if
			end do
		end do
		!$OMP END PARALLEL DO
	
	else
		!$OMP CRITICAL
		call CCTK_Warn(0, "Flux direction not x,y,z")
		!$OMP END CRITICAL
	end if
  
end if


end subroutine Spritz_HLLE_MHD

 /*@@
   @routine    HO_correction
   @date
   @author     Beatrice Giudici
   @desc
   Compute the high order correction to the fluxes, which can be 4th or 6th order accurate
   (Del Zanna et al. 2007, Equation A.21)
   @enddesc
   @calls
   @calledby
@@*/

subroutine HO_correction(f_hat, f, N)
	implicit none
	
	DECLARE_CCTK_PARAMETERS
	
	CCTK_INT :: N, i
	CCTK_REAL, dimension(N) :: f_hat, f
	
	do i = Spritz_stencil, N - Spritz_stencil + 1
		if (approximation_order == 4) then
			f_hat(i) = (13.0d0/12.0d0)*f(i) - (1.0d0/24.0d0)*(f(i-1) + f(i+1))
		
		else if (approximation_order == 6) then
			f_hat(i) = (1067.0d0/960.0d0)*f(i) - (29.0d0/480.0d0)*(f(i-1) + f(i+1)) &
						+ (3.0d0/640.0d0)*(f(i-2) + f(i+2))
		
		else
			call CCTK_Warn(0, "Approximation order must be 4 or 6")
		end if
	end do
	
	return
	
end subroutine HO_correction


 /*@@
   @routine    eigenvalues_MHD
   @date
   @author     Bruno Giacomazzo
   @desc
   Compute the fast and slow magnetosonic eigenvalues
   @enddesc
   @calls
   @calledby
@@*/

subroutine eigenvalues_MHD(alp,b2,b0,cs2,vx,Bx,bsx,betax,uxx,rho,h,wlorentz,lam,x,y,z)
  use Interfaces, only:quartic, zrhqr, cubic
  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL :: alp,b2,b0,cs2,vx,Bx,bsx,betax,uxx,rho,h,wlorentz
  CCTK_REAL,dimension(4) :: lam,dump
  CCTK_REAL,dimension(0:4) :: a
  CCTK_REAL,intent(IN) :: x,y,z

  CCTK_REAL :: btilde0,btildex,btilde2,eps2,rightalfven,leftalfven,det

  CCTK_REAL,dimension(0:2) :: c
  CCTK_REAL,dimension(3) :: W

  character(len=200) warnline

  btilde0 = b0/sqrt(rho*h+b2)
  btildex = bsx/sqrt(rho*h+b2)
  btilde2 = b2/(rho*h+b2)
  eps2=cs2+btilde2-cs2*btilde2


!!$  !special relativistic case
!!$  a(0)=btildex**2*cs2-eps2*wlorentz**2*vx**2-(eps2-1.0D0)*wlorentz**4*vx**4
!!$  a(1)=-2.0D0*btilde0*btildex*cs2+2.0D0*eps2*wlorentz**2*vx+4.0D0*(eps2-1.0D0)*wlorentz**4*vx**3
!!$  a(2)=(btilde0**2-btildex**2)*cs2+wlorentz**2*(6.0D0*wlorentz**2*vx**2+ &
!!$       eps2*(-1.0D0+(1.0D0-6.0D0*wlorentz**2)*vx**2))
!!$  a(3)=2.0D0*btilde0*btildex*cs2+2.0D0*wlorentz**2*(-eps2+2.0D0*(-1.0D0+eps2)*wlorentz**2)*vx
!!$  a(4)=-btilde0**2*cs2+wlorentz**4+eps2*(wlorentz**2-wlorentz**4)


  !approximate values for the fast waves
  !These are obtained imposing Bx=0 and b0=0 (i.e. B^j v_j=0)
  !Analogous to what done by Leismann et al. (2005) A&A 436, 503


  a(0) = (b2 + cs2*h*rho)*(betax**2 - alp**2*uxx) - &
       (-1 + cs2)*h*rho*(betax - alp*vx)**2*wlorentz**2

  a(1) = 2*betax*(b2 + cs2*h*rho) - &
       2*(-1 + cs2)*h*rho*(betax - alp*vx)*wlorentz**2

  a(2) = b2 + h*rho*(cs2 + wlorentz**2 - cs2*wlorentz**2)

  a(3) = 0.0d0

  a(4) = 0.0d0

  det = a(1)**2-4.0d0*a(2)*a(0)

  if(det<0.0d0) then
     det=0.0d0
     !lam(1)=0.0d0
     !lam(2)=0.0d0
     !lam(3)=0.0d0
     !lam(4)=0.0d0
     !return
  end if


  lam(1)=(-a(1)+sqrt(det))/(2.0d0*a(2))
  lam(2)=lam(1)
  lam(3)=(-a(1)-sqrt(det))/(2.0d0*a(2))
  lam(4)=lam(3)

  lam = lam/alp



  if(maxval(abs(lam))>1.0d0 .and. &
      (cs2.le.1.0d0.and.CCTK_EQUALS(Spritz_eos_table,"nuc_eos") .or. &
       .not.CCTK_EQUALS(Spritz_eos_table,"nuc_eos"))) then !LSMOD, if cs2 is faster than light in nuclear tables, don't complain
     !$OMP CRITICAL
     write(warnline,'(a50)') 'Case a4*x^4+a3*x^3+a2*x^2+a1*x+a0==0:'
     call CCTK_WARN(Spritz_NaN_verbose,warnline)
     write(warnline,'(a20,3g16.7)') 'in x,y,z=',x,y,z
     call CCTK_WARN(Spritz_NaN_verbose,warnline)
     write(warnline,'(a20,5g16.7)') 'a(0:4)=',a
     call CCTK_WARN(Spritz_NaN_verbose,warnline)
     write(warnline,'(a20,4g16.7)') 'lam(1:4)=',lam
     call CCTK_WARN(Spritz_NaN_verbose,warnline)
     call CCTK_WARN(Spritz_NaN_verbose,"Eigenvalues exceed the velocity of light!")
     !call CCTK_WARN(Spritz_NaN_verbose,"We set them to zero and move on")
     !$OMP END CRITICAL
     !!lam(1)=0.0d0
     !!lam(2)=0.0d0
     !!lam(3)=0.0d0
     !!lam(4)=0.0d0
     !!return
  end if



end subroutine eigenvalues_MHD



 /*@@
   @routine    num_x_flux_MHD
   @date
   @author     Bruno Giacomazzo
   @desc
   Fluxes along the x direction
   (but can also be used along other directions changing the order of input and output variables)
   @enddesc
   @calls
   @calledby
@@*/


subroutine num_x_flux_MHD(dens,sx,sy,sz,tau,Bx,By,Bz,densf,sxf,syf,szf,tauf,Bxf,Byf,Bzf,&
     velx,vely,velz,avg_betax,avg_betay,avg_betaz,avg_alp,avg_det,w_lorentz,ptot,b0,bxlow,bylow,bzlow)
  implicit none

  DECLARE_CCTK_PARAMETERS

  CCTK_REAL :: dens,sx,sy,sz,tau,Bx,By,Bz
  CCTK_REAL :: densf,sxf,syf,szf,tauf,Bxf,Byf,Bzf
  CCTK_REAL :: avg_betax,avg_betay,avg_betaz,avg_alp,avg_det
  CCTK_REAL :: velx,vely,velz,ptot,b0,bxlow,bylow,bzlow
  CCTK_REAL :: velxtilde,velytilde,velztilde,w_lorentz

  velxtilde = (velx - avg_betax/avg_alp)
  velytilde = (vely - avg_betay/avg_alp)
  velztilde = (velz - avg_betaz/avg_alp)

  !dens, sx, sy, sz, tau are already multiplied by sqrt(avg_det)

  densf = dens*velxtilde
  sxf   = sx*velxtilde + sqrt(avg_det)*(ptot - Spritz_mhd_coeff*bxlow*Bx/w_lorentz)
  syf   = sy*velxtilde - Spritz_mhd_coeff*sqrt(avg_det)*(bylow*Bx/w_lorentz)
  szf   = sz*velxtilde - Spritz_mhd_coeff*sqrt(avg_det)*(bzlow*Bx/w_lorentz)
  tauf  = tau*velxtilde + sqrt(avg_det)*(ptot*velx - Spritz_mhd_coeff*avg_alp*b0*Bx/w_lorentz)
  Bxf   = 0.0d0
  Byf   = sqrt(avg_det)*(velxtilde*By - velytilde*Bx)
  Bzf   = sqrt(avg_det)*(velxtilde*Bz - velztilde*Bx)

end subroutine num_x_flux_MHD

 /*@@
   @routine    SpritzMHDCT
   @date
   @author     Bruno Giacomazzo
   @desc
   Recompute the magnetic field fluxes using a Constrained Transport scheme
   @enddesc
   @calls
   @calledby
@@*/

subroutine SpritzMHDCT(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

!#ifdef _EOS_BASE_INC_
!#undef _EOS_BASE_INC_
!#endif
!#include "EOS_Base.inc"

  integer :: i, j, k, m
  CCTK_REAL :: idx, beta_avec_plus, beta_avec_minus, det_plus, det_minus, Avec_plus, Avec_minus
  ! character(len=200) warnline
  CCTK_REAL :: uxx_plus,uxy_plus,uxz_plus,uyy_plus,uyz_plus,uzz_plus
  CCTK_REAL :: uxx_minus,uxy_minus,uxz_minus,uyy_minus,uyz_minus,uzz_minus

  CCTK_REAL :: curlv_x, curlv_y, curlv_z, curlv, Blowx, Blowy, Blowz, Bnorm_local

  CCTK_REAL :: rho_avg, press_avg, rho_vol_avg, press_vol_avg, dens_avg, tau_avg, sx_avg, sy_avg, sz_avg
  CCTK_REAL :: dBx_avg, dBy_avg, dBz_avg, det, uxx, uxy, uxz, uyy, uyz, uzz
  CCTK_REAL :: psi4pt, Bx_avg, By_avg, Bz_avg, Bxlow_avg, Bylow_avg, Bzlow_avg, B2_avg
  CCTK_REAL :: velx_avg, vely_avg, velz_avg, eps_avg, w_lorentz_avg, delta_w, pmin
  CCTK_REAL :: vlowx, vlowy, vlowz, bst, b2small_local, eos_gamma

  ! Temporary variables to store interpolated values that take into account staggering    !
  ! Basically, we wakt all the following variables to be stored at cell edges, in order   !
  ! to compute the derivatives of scalar potential on cell vertices...                    !
  !                                                                                       !
  CCTK_REAL :: F_p, G_p                                                                   !
  CCTK_REAL :: F_m, G_m                                                                   !
  CCTK_REAL :: Fbeta, Fbeta_1, Fbeta_2                                                    !
  CCTK_REAL :: Axup, Ayup, Azup, beta_avec                                                !
  CCTK_REAL :: alp_vert                                                                   !
  !                                                                                       !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical :: epsnegative, c2perror, outside_ah

  integer :: ll, mm, nn, number_of_zones
  character(len=100) warnline

  !!!!!!
  ! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1),one(1)

  n = 1;keytemp = 0;anyerr = 0;keyerr = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
  one = 1.0d0
  ! end EOS Omni vars
  !!!!!

  !pmin=EOS_Pressure(Spritz_polytrope_handle, Spritz_rho_min, 1.0d0)

  !!!!!!
  xrho = Spritz_rho_min
  call EOS_Omni_press(Spritz_polytrope_handle,keytemp,Eos_rf_prec,n,&
                 xrho,one,xtemp,xye,xpress,keyerr,anyerr)
  pmin = xpress(1)
  !!!!!!

  !trick to obtain the value of eos_ideal_fluid_gamma from thorn EOS_Ideal
  !eos_gamma = EOS_Pressure(Spritz_eos_handle,1.0d0,1.0d0)+1.0d0

  !!!!!!
  call EOS_Omni_press(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                 one,one,xtemp,xye,xpress,keyerr,anyerr)
  eos_gamma = xpress(1) + 1.0d0
  !!!!!!

  if (Spritz_use_zrake_macfadyen_amplification.ne.0) then
     write(warnline,'(a33,g16.7)') 'DEBUG MESSAGE: eos_gamma=', eos_gamma
     call CCTK_WARN(1,warnline)
  end if

!!$  write(*,*) 'SpritzMHDCT'

  if (Spritz_use_zrake_macfadyen_amplification.ne.0) eps_KT=0.0d0

  !we evolve the vector potential and then in c2p we compute B from it
  !$OMP PARALLEL DO PRIVATE(k, j, i, xoffset, yoffset, zoffset, idx, &
  !$OMP beta_avec_plus, beta_avec_minus, det_plus, &
  !$OMP uxx_plus, uxy_plus, uxz_plus, uyy_plus, uyz_plus, uzz_plus, &
  !$OMP det_minus, uxx_minus, uxy_minus, uxz_minus, uyy_minus, uyz_minus, uzz_minus, &
  !$OMP Avec_plus, Avec_minus, &
  !$OMP Blowx, Blowy, Blowz, Bnorm_local, curlv, curlv_x, curlv_y, curlv_z, &
  !$OMP rho_avg, press_avg, rho_vol_avg, press_vol_avg, dens_avg, tau_avg, sx_avg, sy_avg, sz_avg, &
  !$OMP dBx_avg, dBy_avg, dBz_avg, ll, mm, nn, det, uxx, uxy, uxz, uyy, uyz, uzz, &
  !$OMP psi4pt, Bx_avg, By_avg, Bz_avg, Bxlow_avg, Bylow_avg, Bzlow_avg, B2_avg, &
  !$OMP velx_avg, vely_avg, velz_avg, eps_avg, w_lorentz_avg, &
  !$OMP epsnegative, c2perror, outside_ah, delta_w, warnline, &
  !$OMP vlowx, vlowy, vlowz, bst, b2small_local, number_of_zones)
  do k = Spritz_stencil+1, cctk_lsh(3) - Spritz_stencil
     do j = Spritz_stencil+1, cctk_lsh(2) - Spritz_stencil
        do i = Spritz_stencil+1, cctk_lsh(1) - Spritz_stencil

           Avec_x_rhs(i,j,k) = -1.0d0*Ex(i,j,k) !A_x
           Avec_y_rhs(i,j,k) = -1.0d0*Ey(i,j,k) !A_y
           Avec_z_rhs(i,j,k) = -1.0d0*Ez(i,j,k) !A_z

           !Amplify the vector potential using a subgrid model that
           ! takes into account Zrake & MacFadyen results
           ! TO BE CHECKED BECAUSE OF STAGGERING...
           if (Spritz_use_zrake_macfadyen_amplification.ne.0) then
              !v_j
              vlowx = gxx(i,j,k)*velx(i,j,k) + gxy(i,j,k)*vely(i,j,k) + gxz(i,j,k)*velz(i,j,k)
              vlowy = gxy(i,j,k)*velx(i,j,k) + gyy(i,j,k)*vely(i,j,k) + gyz(i,j,k)*velz(i,j,k)
              vlowz = gxz(i,j,k)*velx(i,j,k) + gyz(i,j,k)*vely(i,j,k) + gzz(i,j,k)*velz(i,j,k)
              !compute |B| (used in the saturation factor)
!!!B_j
              Blowx = gxx(i,j,k)*Bx(i,j,k) + gxy(i,j,k)*By(i,j,k) + gxz(i,j,k)*Bz(i,j,k)
              Blowy = gxy(i,j,k)*Bx(i,j,k) + gyy(i,j,k)*By(i,j,k) + gyz(i,j,k)*Bz(i,j,k)
              Blowz = gxz(i,j,k)*Bx(i,j,k) + gyz(i,j,k)*By(i,j,k) + gzz(i,j,k)*Bz(i,j,k)
!!!B^j B_j
              Bnorm_local = sqrt(Bx(i,j,k)*Blowx + By(i,j,k)*Blowy + Bz(i,j,k)*Blowz)
              !b^t
              bst = w_lorentz(i,j,k)*(Bx(i,j,k)*vlowx + By(i,j,k)*vlowy + Bz(i,j,k)*vlowz) !this is b^0 * alp
              !b^mu b_mu
              b2small_local = (Bnorm_local**2 + bst**2)/w_lorentz(i,j,k)**2


              !compute curl(v)
              curlv_x = DIFF_Y(velz,i,j,k)-DIFF_Z(vely,i,j,k)
              curlv_y = DIFF_Z(velx,i,j,k)-DIFF_X(velz,i,j,k)
              curlv_z = DIFF_X(vely,i,j,k)-DIFF_Y(velx,i,j,k)
              curlv = sqrt(curlv_x**2 + curlv_y**2 + curlv_z**2)

              if (curlv>0.5) then

!!$                 if(rho(i,j,k) > 100.0d0*Spritz_rho_min) then
!!$                    curlv = sqrt(curlv_x**2 + curlv_y**2 + curlv_z**2)
!!$                 else
!!$                    curlv = 0.0d0
!!$                 end if

                 !compute kinetic energy of the turbulent flow (see appendix of http://arxiv.org/abs/1302.7306)
                 rho_avg       = 0.0d0 !initialize to zero
                 press_avg     = 0.0d0
                 w_lorentz_avg = 0.0d0
                 velx_avg      = 0.0d0
                 vely_avg      = 0.0d0
                 velz_avg      = 0.0d0
                 eps_avg       = 0.0d0
                 rho_vol_avg   = 0.0d0
                 press_vol_avg = 0.0d0
                 dens_avg      = 0.0d0
                 tau_avg       = 0.0d0
                 sx_avg        = 0.0d0
                 sy_avg        = 0.0d0
                 sz_avg        = 0.0d0
                 dBx_avg       = 0.0d0
                 dBy_avg       = 0.0d0
                 dBz_avg       = 0.0d0
                 !compute the average quantities
                 do ll=-zrake_stencil,zrake_stencil,1
                    do mm=-zrake_stencil,zrake_stencil,1
                       do nn=-zrake_stencil,zrake_stencil,1
                          !primitive averages to be used as initial guesses in the C2P routine
                          rho_avg       = rho_avg   + rho(i+ll,j+mm,k+nn)
                          press_avg     = press_avg + press(i+ll,j+mm,k+nn)
                          w_lorentz_avg = w_lorentz_avg + w_lorentz(i+ll,j+mm,k+nn)
                          velx_avg      = velx_avg + velx(i+ll,j+mm,k+nn)
                          vely_avg      = vely_avg + vely(i+ll,j+mm,k+nn)
                          velz_avg      = velz_avg + velz(i+ll,j+mm,k+nn)
                          eps_avg       = eps_avg  + eps(i+ll,j+mm,k+nn)
                          !we need to keep these in memory to compute delta_w
                          rho_vol_avg   = rho_vol_avg   + rho(i+ll,j+mm,k+nn)
                          press_vol_avg = press_vol_avg + press(i+ll,j+mm,k+nn)
                          !Volume average of conserved variables
                          dens_avg      = dens_avg  + dens(i+ll,j+mm,k+nn)
                          tau_avg       = tau_avg   + tau(i+ll,j+mm,k+nn)
                          sx_avg        = sx_avg    + scon(i+ll,j+mm,k+nn,1)
                          sy_avg        = sy_avg    + scon(i+ll,j+mm,k+nn,2)
                          sz_avg        = sz_avg    + scon(i+ll,j+mm,k+nn,3)
                          dBx_avg       = dBx_avg   + dBx(i+ll,j+mm,k+nn)
                          dBy_avg       = dBy_avg   + dBy(i+ll,j+mm,k+nn)
                          dBz_avg       = dBz_avg   + dBz(i+ll,j+mm,k+nn)
                       end do
                    end do
                 end do
                 number_of_zones = (2*zrake_stencil+1)**3
                 rho_avg       = rho_avg/number_of_zones
                 press_avg     = press_avg/number_of_zones
                 w_lorentz_avg = w_lorentz_avg/number_of_zones
                 velx_avg      = velx_avg/number_of_zones
                 vely_avg      = vely_avg/number_of_zones
                 velz_avg      = velz_avg/number_of_zones
                 eps_avg       = eps_avg/number_of_zones
                 rho_vol_avg   = rho_vol_avg/number_of_zones
                 press_vol_avg = press_vol_avg/number_of_zones
                 dens_avg      = dens_avg/number_of_zones
                 tau_avg       = tau_avg/number_of_zones
                 sx_avg        = sx_avg/number_of_zones
                 sy_avg        = sy_avg/number_of_zones
                 sz_avg        = sz_avg/number_of_zones
                 dBx_avg       = dBx_avg/number_of_zones
                 dBy_avg       = dBy_avg/number_of_zones
                 dBz_avg       = dBz_avg/number_of_zones

                 !Get primitive variables from averaged conservatives
                 if (conformal_state .eq. 0) then
                    call SpatialDeterminant(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),&
                         gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det)
                    call UpperMetric(uxx,uxy,uxz,uyy,uyz,uzz,det,&
                         gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),&
                         gyz(i,j,k),gzz(i,j,k))
                 else
                    psi4pt = psi(i,j,k)**4
                    call SpatialDeterminant(psi4pt*gxx(i,j,k),psi4pt*gxy(i,j,k),&
                         psi4pt*gxz(i,j,k),psi4pt*gyy(i,j,k),&
                         psi4pt*gyz(i,j,k),psi4pt*gzz(i,j,k),det)
                    call UpperMetric(uxx,uxy,uxz,uyy,uyz,uzz,det,&
                         psi4pt*gxx(i,j,k),psi4pt*gxy(i,j,k),psi4pt*gxz(i,j,k),&
                         psi4pt*gyy(i,j,k),psi4pt*gyz(i,j,k),psi4pt*gzz(i,j,k))
                 end if

                 !B^j
                 Bx_avg = dBx_avg/sqrt(det)
                 By_avg = dBy_avg/sqrt(det)
                 Bz_avg = dBz_avg/sqrt(det)
                 !B_j
                 Bxlow_avg = Bx_avg*gxx(i,j,k)+By_avg*gxy(i,j,k)+Bz_avg*gxz(i,j,k)
                 Bylow_avg = Bx_avg*gxy(i,j,k)+By_avg*gyy(i,j,k)+Bz_avg*gyz(i,j,k)
                 Bzlow_avg = Bx_avg*gxz(i,j,k)+By_avg*gyz(i,j,k)+Bz_avg*gzz(i,j,k)

                 !B^j B_j
                 B2_avg = gxx(i,j,k) * Bx_avg**2 + &
                      gyy(i,j,k) * By_avg**2 + &
                      gzz(i,j,k) * Bz_avg**2 + &
                      2.d0 * gxy(i,j,k) * Bx_avg * By_avg + &
                      2.d0 * gxz(i,j,k) * Bx_avg * Bz_avg + &
                      2.d0 * gyz(i,j,k) * By_avg * Bz_avg

                 epsnegative = .false.
                 c2perror = .false.
                 outside_ah = .false.

                 call Con2Prim_pt_mhd(Spritz_eos_handle, dens_avg,sx_avg,sy_avg, &
                      sz_avg,dBx_avg,dBy_avg,dBz_avg,&
                      tau_avg,rho_avg,velx_avg,vely_avg, &
                      velz_avg,eps_avg,press_avg,w_lorentz_avg, &
                      uxx,uxy,uxz,uyy,uyz,uzz,det,x(i,j,k),y(i,j,k), &
                      z(i,j,k),r(i,j,k),epsnegative,Spritz_rho_min, &
                      Bx_avg,By_avg,Bz_avg,B2_avg,Bxlow_avg,Bylow_avg,Bzlow_avg,Spritz_polytrope_handle,pmin, &
                      Spritz_reflevel, c2perror, outside_ah, mhd_c2p_mask(i,j,k))

                 if (epsnegative.or.c2perror) then
                    !$OMP CRITICAL
                    if (epsnegative) call CCTK_WARN(1,"epsnegative")
                    if (c2perror) call CCTK_WARN(1,"c2perror")
                    call CCTK_WARN(1,'Spritz_HLLE_MHD: stopping the code.')
                    write(warnline,'(a28,i2)') 'on carpet reflevel: ',Spritz_reflevel
                    call CCTK_WARN(1,warnline)
                    write(warnline,'(a20,3g16.7)') 'xyz location: ',&
                         x(i,j,k),y(i,j,k),z(i,j,k)
                    call CCTK_WARN(1,warnline)
                    write(warnline,'(a20,g16.7)') 'radius: ',r(i,j,k)
                    call CCTK_WARN(1,warnline)
                    write(warnline,'(a20,3g16.7)') 'velocities: ',&
                         velx_avg,vely_avg,velz_avg
                    call CCTK_WARN(1,warnline)
                    write(warnline,'(a20,3g16.7)') 'magnetic field: ',&
                         Bx_avg,By_avg,Bz_avg
                    call CCTK_WARN(1,warnline)
                    write(warnline,'(a30,2g16.7)') 'dens_avg, rho_avg: ',&
                         dens_avg, rho_avg
                    call CCTK_WARN(1,warnline)
                    write(warnline,'(a30,3g16.7)') 'tau_avg, eps_avg, press_avg: ',&
                         tau_avg, eps_avg, press_avg
                    call CCTK_WARN(1,warnline)
                    write(warnline,'(a20,g16.7)') 'w_lorentz_avg: ',&
                         w_lorentz_avg
                    call CCTK_WARN(1,warnline)
                    call CCTK_WARN(1,"C2P Failed in HLLE_MHD when applying Zrake-MacFadyen source term")
                    !$OMP END CRITICAL
                    curlv = 0.0d0
                 end if

                 !W=rho+rho*eps=rho+press/(eos_gamma-1)
                 !for gamma=2, W=rho+press
                 delta_w = (rho_avg + press_avg/(eos_gamma-1.0d0)) - (rho_vol_avg + press_vol_avg/(eos_gamma-1.0d0))
                 eps_KT(i,j,k) = zrake_c2*delta_w !min(delta_w,1.0e-6) !b2_max=1.0e-6 which means B_max~1.0e-3~1e17G

                 if (delta_w<=0.0d0) then
                    curlv = 0.0d0
                 else
                    !multiply for saturation factor (hard coded to 1e16 G) and density (weighted vorticity)
                    !eps_KT = 1.0e-4
                    curlv = zrake_c1 * max(curlv-0.5,0.0d0) * max(1-1.0d4*Spritz_rho_min/rho(i,j,k),0.0d0) * max(1.0d0-b2small_local/eps_KT(i,j,k),0.0d0)
                 end if

                 !update source terms for the vector potential
                 Avec_x_rhs(i,j,k) = Avec_x_rhs(i,j,k) + curlv*Avec_x(i,j,k) !A_x
                 Avec_y_rhs(i,j,k) = Avec_y_rhs(i,j,k) + curlv*Avec_y(i,j,k) !A_y
                 Avec_z_rhs(i,j,k) = Avec_z_rhs(i,j,k) + curlv*Avec_z(i,j,k) !A_z
              end if
           end if

        end do
     end do
  end do
  !$OMP END PARALLEL DO
  !$OMP BARRIER

  !Initialize fnctions for computing update terms to 0
  F_p = 0.0d0
  G_p = 0.0d0
  F_m = 0.0d0
  G_m = 0.0d0
  Fbeta = 0.0d0
  Fbeta_1 = 0.0d0
  Fbeta_2 = 0.0d0
  Fx = 0.0d0
  Fy = 0.0d0
  Fz = 0.0d0
  Fbetax = 0.0d0
  Fbetay = 0.0d0
  Fbetaz = 0.0d0
  G = 0.0d0
  Ax_center = 0.0d0
  Ay_center = 0.0d0
  Az_center = 0.0d0
  psi_mhd_center = 0.0d0
  Axup = 0.0d0
  Ayup = 0.0d0
  Azup = 0.0d0
  beta_avec = 0.0d0
  alp_vert = 0.0d0

  if((Spritz_mhd_fluxCD_Lorenz.ne.0)) then

    !call CCTK_INFO("Using Lorenz Gauge in HLLE_MHD...entering...")

    ! Nested do loop for evaluation of update terms for Avec and psi_mhd in the
    ! generalized lorenz gauge case.
    !$OMP PARALLEL DO PRIVATE(k, j, i, det, uxx, uxy, uxz, uyy, uyz, uzz, &
    !$OMP Axup, Ayup, Azup, beta_avec)
    do k = Spritz_stencil-1, cctk_lsh(3) - Spritz_stencil + 3
      do j = Spritz_stencil-1, cctk_lsh(2) - Spritz_stencil + 3
         do i = Spritz_stencil-1, cctk_lsh(1) - Spritz_stencil + 3

                 ! Interpolating functions at cells center...
                 ! NOTE: Because of staggered prolongation we have:
                 !       A_x(i,j,k) = A_x(i,j+1/2,k+1/2)
                 !       A_y(i,j,k) = A_y(i+1/2,j,k+1/2)
                 !       A_z(i,j,k) = A_z(i+1/2,j+1/2,k)
                 !       Read the comments on "InterpEdges2Center" function...
                 Ax_center(i,j,k) = InterpEdges2Center(Avec_x(i,j-1,k-1),Avec_x(i,j,k-1),&
                                                       Avec_x(i,j-1,k),Avec_x(i,j,k))
                 Ay_center(i,j,k) = InterpEdges2Center(Avec_y(i-1,j,k-1),Avec_y(i,j,k-1),&
                                                       Avec_y(i-1,j,k),Avec_y(i,j,k))
                 Az_center(i,j,k) = InterpEdges2Center(Avec_z(i-1,j-1,k),Avec_z(i,j-1,k),&
                                                       Avec_z(i-1,j,k),Avec_z(i,j,k))
                 ! NOTE: Because of staggered prolongation we have:
                 !       psi_mhd(i,j,k) = psi_mhd(i+1/2,j+1/2,k+1/2)
                 !       Read the comments on "InterpVertices2Centers" function...
                 psi_mhd_center(i,j,k) = InterpVertices2Centers(psi_mhd(i-1,j-1,k-1),psi_mhd(i,j-1,k-1),&
                                                                psi_mhd(i-1,j,k-1),psi_mhd(i,j,k-1),&
                                                                psi_mhd(i-1,j-1,k),psi_mhd(i,j-1,k),&
                                                                psi_mhd(i-1,j,k),psi_mhd(i,j,k))

                 ! Now let's compute the values of function that must be derived in order to
                 ! obtain the update terms for vector and scalar potential. In particular,
                 ! Fi and Fbetai are the function used for the scalar potential's update term, while G
                 ! is the function for the components of the vector potential's update terms.
                 ! We want all of them computed at cells' centers.
                 ! Namely:
                 ! Fx = (alpha * sqrt(det) * A^x )
                 ! Fy = (alpha * sqrt(det) * A^y )
                 ! Fz = (alpha * sqrt(det) * A^z )
                 ! Fbetax = (- beta^x * psi_mhd)
                 ! Fbetay = (- beta^y * psi_mhd)
                 ! Fbetaz = (- beta^z * psi_mhd)
                 ! G = (alpha * psi_mhd / sqrt(det) - beta^j * A_j )

                 call SpatialDeterminant(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),&
                                         gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det)

                 call UpperMetric(uxx,uxy,uxz,uyy,uyz,uzz,det,&
                                  gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),&
                                  gyz(i,j,k),gzz(i,j,k))

                 Axup = uxx*Ax_center(i,j,k) + uxy*Ay_center(i,j,k) + uxz*Az_center(i,j,k)
                 Ayup = uxy*Ax_center(i,j,k) + uyy*Ay_center(i,j,k) + uyz*Az_center(i,j,k)
                 Azup = uxz*Ax_center(i,j,k) + uyz*Ay_center(i,j,k) + uzz*Az_center(i,j,k)

                 beta_avec = betax(i,j,k)*Ax_center(i,j,k) + &
                             betay(i,j,k)*Ay_center(i,j,k) + &
                             betaz(i,j,k)*Az_center(i,j,k)

                 Fx(i,j,k) = alp(i,j,k)*sqrt(det)*Axup
                 Fy(i,j,k) = alp(i,j,k)*sqrt(det)*Ayup
                 Fz(i,j,k) = alp(i,j,k)*sqrt(det)*Azup
                 Fbetax(i,j,k) = -betax(i,j,k)*psi_mhd_center(i,j,k)
                 Fbetay(i,j,k) = -betay(i,j,k)*psi_mhd_center(i,j,k)
                 Fbetaz(i,j,k) = -betaz(i,j,k)*psi_mhd_center(i,j,k)
                 G(i,j,k) = alp(i,j,k)*psi_mhd_center(i,j,k)/sqrt(det) - beta_avec

          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !Make omp process synchronize, because Cactus doens't do it
    ! authomatically within a single subroutine. This is done,
    ! because the arrays defined in the previous loop are used
    ! in the next loop to compute update terms.
    !$OMP BARRIER

    ! Initialize rhs of scalar potential to 0
    psi_mhdrhs = 0.0d0

    ! From the previous nested do loop, we have populated arrays that must be used for obtaining
    ! the updated terms in case of Lorenz Gauge. The values within these arrays are computed
    ! at cells' centers. Now, we will interpolate these values at cells's edges (for F)
    ! and vertices (for G) and use them for numerical differentiation, finally obtaining the
    ! update terms computed at vertices (psi_mhd) or edges (A).
    !$OMP PARALLEL DO PRIVATE(k, j, i, m, idx, F_p, G_p, F_m, G_m,&
    !$OMP Fbeta, Fbeta_1, Fbeta_2, alp_vert)
    do k = Spritz_stencil+1, cctk_lsh(3) - Spritz_stencil
       do j = Spritz_stencil+1, cctk_lsh(2) - Spritz_stencil
          do i = Spritz_stencil+1, cctk_lsh(1) - Spritz_stencil
             !use the Lorenz gauge to evolve the vector potential
             do m = 3, 1, -1

                ! Here we apply centered finite difference for numerical differentiation
                idx = 1.0d0/(CCTK_DELTA_SPACE(m))

                ! G_p should be in principle the same in every direction, but for safety reason we
                ! just recompute that in every direction. That should not be too expensive...
                ! G_m, F_m , Fbeta_2, Fbeta_1 and Fbeta depend on the differentiation direction.
                if (m==1) then

                  ! Interpolate function G on vertex (i+1/2,j+1/2,k+1/2) to obtain G_p
                  G_p = InterpCenters2Vertex(G(i,j,k),G(i,j+1,k),&
                                             G(i,j,k+1),G(i,j+1,k+1),&
                                             G(i+1,j,k),G(i+1,j+1,k),&
                                             G(i+1,j,k+1),G(i+1,j+1,k+1))   ! => G(i+1/2,j+1/2,k+1/2)

                  ! Interpolate functions Fx from cells' centers to cells' edges
                  ! Keep i index fixed for the x direction:
                  F_m = InterpCenters2Edge(Fx(i,j,k),Fx(i,j+1,k),&
                                           Fx(i,j,k+1),Fx(i,j+1,k+1))     ! ==> Fx(i,j+1/2,k+1/2)
                  F_p = InterpCenters2Edge(Fx(i+1,j,k),Fx(i+1,j+1,k),&
                                           Fx(i+1,j,k+1),Fx(i+1,j+1,k+1)) ! ==> Fx(i+1,j+1/2,k+1/2)

                  ! Interpolate functions Fbetax from cells' centers to cells' vertices
                  ! i.e. compute the terms for the 2nd-order upwinded derivative in the x direction:
                  if (betax(i,j,k) .lt. 0.0d0) then
                    Fbeta = InterpCenters2Vertex(Fbetax(i,j,k),Fbetax(i,j+1,k),&
                                                 Fbetax(i,j,k+1),Fbetax(i,j+1,k+1),&
                                                 Fbetax(i+1,j,k),Fbetax(i+1,j+1,k),&
                                                 Fbetax(i+1,j,k+1),Fbetax(i+1,j+1,k+1))          ! ==> Fbetax(i+1/2,j+1/2,k+1/2)
                    Fbeta_1 = InterpCenters2Vertex(Fbetax(i-1,j,k),Fbetax(i-1,j+1,k),&
                                                   Fbetax(i-1,j,k+1),Fbetax(i-1,j+1,k+1),&
                                                   Fbetax(i,j,k),Fbetax(i,j+1,k),&
                                                   Fbetax(i,j,k+1),Fbetax(i,j+1,k+1))            ! ==> Fbetax(i-1/2,j+1/2,k+1/2)
                    Fbeta_2 = InterpCenters2Vertex(Fbetax(i-2,j,k),Fbetax(i-2,j+1,k),&
                                                   Fbetax(i-2,j,k+1),Fbetax(i-2,j+1,k+1),&
                                                   Fbetax(i-1,j,k),Fbetax(i-1,j+1,k),&
                                                   Fbetax(i-1,j,k+1),Fbetax(i-1,j+1,k+1))        ! ==> Fbetax(i-3/2,j+1/2,k+1/2)
                  else
                    Fbeta = -1.0d0*InterpCenters2Vertex(Fbetax(i,j,k),Fbetax(i,j+1,k),&
                                                        Fbetax(i,j,k+1),Fbetax(i,j+1,k+1),&
                                                        Fbetax(i+1,j,k),Fbetax(i+1,j+1,k),&
                                                        Fbetax(i+1,j,k+1),Fbetax(i+1,j+1,k+1))   ! ==> - Fbetax(i+1/2,j+1/2,k+1/2)
                    Fbeta_1 = -1.0d0*InterpCenters2Vertex(Fbetax(i+1,j,k),Fbetax(i+1,j+1,k),&
                                                          Fbetax(i+1,j,k+1),Fbetax(i+1,j+1,k+1),&
                                                          Fbetax(i+2,j,k),Fbetax(i+2,j+1,k),&
                                                          Fbetax(i+2,j,k+1),Fbetax(i+2,j+1,k+1)) ! ==> - Fbetax(i+3/2,j+1/2,k+1/2)
                    Fbeta_2 = -1.0d0*InterpCenters2Vertex(Fbetax(i+2,j,k),Fbetax(i+2,j+1,k),&
                                                          Fbetax(i+2,j,k+1),Fbetax(i+2,j+1,k+1),&
                                                          Fbetax(i+3,j,k),Fbetax(i+3,j+1,k),&
                                                          Fbetax(i+3,j,k+1),Fbetax(i+3,j+1,k+1)) ! ==> - Fbetax(i+5/2,j+1/2,k+1/2)
                  end if

                  ! Interpolate function G on vertex (i-1/2,j+1/2,k+1/2) to obtain G_m on x direction
                  G_m = InterpCenters2Vertex(G(i-1,j,k),G(i-1,j+1,k),&
                                             G(i-1,j,k+1),G(i-1,j+1,k+1),&
                                             G(i,j,k),G(i,j+1,k),&
                                             G(i,j,k+1),G(i,j+1,k+1))     ! => G(i-1/2,j+1/2,k+1/2)

                  ! Numerical differentition for Avec_x_rhs
                  Avec_x_rhs(i,j,k) = Avec_x_rhs(i,j,k) - idx*( G_p - G_m )

                  ! Numerical differentition for psi_mhd_rhs (F_p, F_m, Fbeta, Fbeta_1 and Fbeta_2 here, depend on the direction)
                  psi_mhdrhs(i,j,k) = psi_mhdrhs(i,j,k) - idx*( F_p - F_m ) - idx*0.5d0*( Fbeta_2 - 4.0d0*Fbeta_1 + 3.0d0*Fbeta )

                  ! NOTE: m=1 is the last flux_direction, so only here we should add the lorenz_damp_factor term
                  !       to the computed value of psi_mhdrhs(i,j,k).
                  ! Computing the Generalized Lorenz Gauge (the one with the damping factor) to the scalar
                  ! potential update term. This term can be written as:
                  !
                  !     - lorenz_damp_fac * alp(i,j,k) * psi_mhd(i,j,k)
                  !
                  ! where alp(i,j,k) is computed at cells' centers. Thus we should interpolate
                  ! these last values at cells' vertices, before adding this term
                  alp_vert = InterpCenters2Vertex(alp(i,j,k),alp(i,j+1,k),&
                                                  alp(i,j,k+1),alp(i,j+1,k+1),&
                                                  alp(i+1,j,k),alp(i+1,j+1,k),&
                                                  alp(i+1,j,k+1),alp(i+1,j+1,k+1)) ! ==> alp(i+1/2,j+1/2,k+1/2)

                  psi_mhdrhs(i,j,k) = psi_mhdrhs(i,j,k) - (lorenz_damp_fac*alp_vert*psi_mhd(i,j,k))

                else if (m==2) then

                  ! Interpolate function G on vertex (i+1/2,j+1/2,k+1/2) to obtain G_p
                  G_p = InterpCenters2Vertex(G(i,j,k),G(i,j+1,k),&
                                             G(i,j,k+1),G(i,j+1,k+1),&
                                             G(i+1,j,k),G(i+1,j+1,k),&
                                             G(i+1,j,k+1),G(i+1,j+1,k+1))   ! => G(i+1/2,j+1/2,k+1/2)

                  ! Interpolate functions Fy from cells' centers to cells' edges
                  ! Keep j index fixed for the y direction:
                  F_m = InterpCenters2Edge(Fy(i,j,k),Fy(i,j,k+1),&
                                           Fy(i+1,j,k),Fy(i+1,j,k+1))     ! ==> Fy(i+1/2,j,k+1/2)
                  F_p = InterpCenters2Edge(Fy(i,j+1,k),Fy(i,j+1,k+1),&
                                           Fy(i+1,j+1,k),Fy(i+1,j+1,k+1)) ! ==> Fy(i+1/2,j+1,k+1/2)

                  ! Interpolate functions Fbetay from cells' centers to cells' vertices
                  ! i.e. compute the terms for the 2nd-order upwinded derivative in the y direction:
                  if (betay(i,j,k) .lt. 0.0d0) then
                    Fbeta = InterpCenters2Vertex(Fbetay(i,j,k),Fbetay(i+1,j,k),&
                                                 Fbetay(i,j,k+1),Fbetay(i+1,j,k+1),&
                                                 Fbetay(i,j+1,k),Fbetay(i+1,j+1,k),&
                                                 Fbetay(i,j+1,k+1),Fbetay(i+1,j+1,k+1))          ! ==> Fbetay(i+1/2,j+1/2,k+1/2)
                    Fbeta_1 = InterpCenters2Vertex(Fbetay(i,j-1,k),Fbetay(i+1,j-1,k),&
                                                   Fbetay(i,j-1,k+1),Fbetay(i+1,j-1,k+1),&
                                                   Fbetay(i,j,k),Fbetay(i+1,j,k),&
                                                   Fbetay(i,j,k+1),Fbetay(i+1,j,k+1))            ! ==> Fbetay(i+1/2,j-1/2,k+1/2)
                    Fbeta_2 = InterpCenters2Vertex(Fbetay(i,j-2,k),Fbetay(i+1,j-2,k),&
                                                   Fbetay(i,j-2,k+1),Fbetay(i+1,j-2,k+1),&
                                                   Fbetay(i,j-1,k),Fbetay(i+1,j-1,k),&
                                                   Fbetay(i,j-1,k+1),Fbetay(i+1,j-1,k+1))        ! ==> Fbetay(i+1/2,j-3/2,k+1/2)
                  else
                    Fbeta = -1.0d0*InterpCenters2Vertex(Fbetay(i,j,k),Fbetay(i,j,k+1),&
                                                        Fbetay(i+1,j,k),Fbetay(i+1,j,k+1),&
                                                        Fbetay(i,j+1,k),Fbetay(i,j+1,k+1),&
                                                        Fbetay(i+1,j+1,k),Fbetay(i+1,j+1,k+1))   ! ==> - Fbetay(i+1/2,j+1/2,k+1/2)
                    Fbeta_1 = -1.0d0*InterpCenters2Vertex(Fbetay(i,j+1,k),Fbetay(i+1,j+1,k),&
                                                          Fbetay(i,j+1,k+1),Fbetay(i+1,j+1,k+1),&
                                                          Fbetay(i,j+2,k),Fbetay(i+1,j+2,k),&
                                                          Fbetay(i,j+2,k+1),Fbetay(i+1,j+2,k+1)) ! ==> - Fbetay(i+1/2,j+3/2,k+1/2)
                    Fbeta_2 = -1.0d0*InterpCenters2Vertex(Fbetay(i,j+2,k),Fbetay(i,j+2,k+1),&
                                                          Fbetay(i+1,j+2,k),Fbetay(i+1,j+2,k+1),&
                                                          Fbetay(i,j+3,k),Fbetay(i,j+3,k+1),&
                                                          Fbetay(i+1,j+3,k),Fbetay(i+1,j+3,k+1)) ! ==> - Fbetay(i+1/2,j+5/2,k+1/2)
                  end if

                  ! Interpolate function G on vertex (i+1/2,j-1/2,k+1/2) to obtain G_m on y direction
                  G_m = InterpCenters2Vertex(G(i,j-1,k),G(i+1,j-1,k),&
                                             G(i,j-1,k+1),G(i+1,j-1,k+1),&
                                             G(i,j,k),G(i,j,k+1),&
                                             G(i+1,j,k),G(i+1,j,k+1))     ! => G(i+1/2,j-1/2,k+1/2)

                  ! Numerical differentition for Avec_y_rhs
                  Avec_y_rhs(i,j,k) = Avec_y_rhs(i,j,k) - idx*( G_p - G_m )

                  ! Numerical differentition for psi_mhd_rhs (F_p, F_m, Fbeta, Fbeta_1 and Fbeta_2 here, depend on the direction)
                  psi_mhdrhs(i,j,k) = psi_mhdrhs(i,j,k) - idx*( F_p - F_m ) - idx*0.5d0*( Fbeta_2 - 4.0d0*Fbeta_1 + 3.0d0*Fbeta )

                else if (m==3) then

                  ! Interpolate function G on vertex (i+1/2,j+1/2,k+1/2) to obtain G_p
                  G_p = InterpCenters2Vertex(G(i,j,k),G(i,j+1,k),&
                                             G(i,j,k+1),G(i,j+1,k+1),&
                                             G(i+1,j,k),G(i+1,j+1,k),&
                                             G(i+1,j,k+1),G(i+1,j+1,k+1))   ! => G(i+1/2,j+1/2,k+1/2)

                  ! Interpolate functions Fz from cells' centers to cells' edges
                  ! Keep k index fixed for the z direction:
                  F_m = InterpCenters2Edge(Fz(i,j,k),Fz(i,j+1,k),&
                                           Fz(i+1,j,k),Fz(i+1,j+1,k))     ! ==> Fz(i+1/2,j+1/2,k)
                  F_p = InterpCenters2Edge(Fz(i,j,k+1),Fz(i,j+1,k+1),&
                                           Fz(i+1,j,k+1),Fz(i+1,j+1,k+1)) ! ==> Fz(i+1/2,j+1/2,k+1)

                  ! Interpolate functions Fbetaz from cells' centers to cells' vertices
                  ! i.e. compute the terms for the 2nd order upwinded derivative in the z direction:
                  if (betaz(i,j,k) .lt. 0.0d0) then
                    Fbeta = InterpCenters2Vertex(Fbetaz(i,j,k),Fbetaz(i,j+1,k),&
                                                 Fbetaz(i+1,j,k),Fbetaz(i+1,j+1,k),&
                                                 Fbetaz(i,j,k+1),Fbetaz(i,j+1,k+1),&
                                                 Fbetaz(i+1,j,k+1),Fbetaz(i+1,j+1,k+1))          ! ==> Fbetaz(i+1/2,j+1/2,k+1/2)
                    Fbeta_1 = InterpCenters2Vertex(Fbetaz(i,j,k-1),Fbetaz(i,j+1,k-1),&
                                                   Fbetaz(i+1,j,k-1),Fbetaz(i+1,j+1,k-1),&
                                                   Fbetaz(i,j,k),Fbetaz(i,j+1,k),&
                                                   Fbetaz(i+1,j,k),Fbetaz(i+1,j+1,k))            ! ==> Fbetaz(i+1/2,j+1/2,k-1/2)
                    Fbeta_2 = InterpCenters2Vertex(Fbetaz(i,j,k-2),Fbetaz(i,j+1,k-2),&
                                                   Fbetaz(i+1,j,k-2),Fbetaz(i+1,j+1,k-2),&
                                                   Fbetaz(i,j,k-1),Fbetaz(i,j+1,k-1),&
                                                   Fbetaz(i+1,j,k-1),Fbetaz(i+1,j+1,k-1))        ! ==> Fbetaz(i+1/2,j+1/2,k-3/2)
                  else
                    Fbeta = -1.0d0*InterpCenters2Vertex(Fbetaz(i,j,k),Fbetaz(i,j+1,k),&
                                                        Fbetaz(i+1,j,k),Fbetaz(i+1,j+1,k),&
                                                        Fbetaz(i,j,k+1),Fbetaz(i,j+1,k+1),&
                                                        Fbetaz(i+1,j,k+1),Fbetaz(i+1,j+1,k+1))   ! ==> - Fbetaz(i+1/2,j+1/2,k+1/2)
                    Fbeta_1 = -1.0d0*InterpCenters2Vertex(Fbetaz(i,j,k+1),Fbetaz(i,j+1,k+1),&
                                                          Fbetaz(i+1,j,k+1),Fbetaz(i+1,j+1,k+1),&
                                                          Fbetaz(i,j,k+2),Fbetaz(i,j+1,k+2),&
                                                          Fbetaz(i+1,j,k+2),Fbetaz(i+1,j+1,k+2)) ! ==> - Fbetaz(i+1/2,j+1/2,k+3/2)
                    Fbeta_2 = -1.0d0*InterpCenters2Vertex(Fbetaz(i,j,k+2),Fbetaz(i,j+1,k+2),&
                                                          Fbetaz(i+1,j,k+2),Fbetaz(i+1,j+1,k+2),&
                                                          Fbetaz(i,j,k+3),Fbetaz(i,j+1,k+3),&
                                                          Fbetaz(i+1,j,k+3),Fbetaz(i+1,j+1,k+3)) ! ==> - Fbetaz(i+1/2,j+1/2,k+5/2)
                  end if

                  ! Interpolate function G on vertex (i+1/2,j+1/2,k-1/2) to obtain G_m on z direction
                  G_m = InterpCenters2Vertex(G(i,j,k-1),G(i+1,j,k-1),&
                                             G(i,j+1,k-1),G(i+1,j+1,k-1),&
                                             G(i,j,k),G(i,j+1,k),&
                                             G(i+1,j,k),G(i+1,j+1,k))     ! => G(i+1/2,j+1/2,k-1/2)

                  ! Numerical differentition for Avec_z_rhs
                  Avec_z_rhs(i,j,k) = Avec_z_rhs(i,j,k) - idx*( G_p - G_m )

                  ! Numerical differentition for psi_mhd_rhs (F_p, F_m, Fbeta, Fbeta_1 and Fbeta_2 here, depend on the direction)
                  ! NOTE: m=3 is the first in the flux_direction loop, so here we simply assing to psi_mhdrhs(i,j,k) its
                  !       first value.
                  psi_mhdrhs(i,j,k) = - idx*( F_p - F_m ) - idx*0.5d0*( Fbeta_2 - 4.0d0*Fbeta_1 + 3.0d0*Fbeta )
                else
                  !$OMP CRITICAL
                  call CCTK_WARN(0,"Flux direction not x,y,z")
                  !$OMP END CRITICAL
                end if

             end do !loop on m end
          end do !loop on i end
       end do !loop on j end
    end do !loop on k end
    !$OMP END PARALLEL DO
  !else

  !  call CCTK_INFO("Adopting Algebraic Gauge in HLLE_MHD...exiting")

  end if

  Ex = 0.0d0
  Ey = 0.0d0
  Ez = 0.0d0
  Bxrhs = 0.0d0
  Byrhs = 0.0d0
  Bzrhs = 0.0d0
  return

  contains
    function InterpCenters2Edge(C1,C2,C3,C4) result(Centers2Edge)
      implicit none
      CCTK_REAL, intent(IN) :: C1, C2, C3, C4
      CCTK_REAL :: Centers2Edge

      ! Given values at points at 4 neighbouring cells's centers, this function computes the linear
      ! interpolation of these values at the edge between them. Namely, let F(i,j,k) be computed
      ! at the cell center (i,j,k). Then the following holds:
      !
      ! F(i,j+1/2,k) = 1/2 * ( F(i,j,k) + F(i,j+1,k) )
      ! F(i,j+1/2,k+1/2) = 1/2 * ( F(i,j+1/2,k) + F(i,j+1/2,k+1) )
      !
      ! Straightforward substitutions lead to:
      !
      ! F(i,j+1/2,k+1/2) = 1/4 * ( F(i,j,k) + F(i,j+1,k) + F(i,j,k+1) + F(i,j+1,k+1) )

      Centers2Edge = 0.25d0 * ( C1 + C2 + C3 + C4 )

    end function InterpCenters2Edge

    function InterpCenters2Vertex(C1,C2,C3,C4,C5,C6,C7,C8) result(Centers2Vertex)
      implicit none
      CCTK_REAL, intent(IN) :: C1, C2, C3, C4, C5, C6, C7, C8
      CCTK_REAL :: Centers2Vertex

      ! Given values computed at 8 neighbouring cells' centers, this function computes the linear
      ! interpolation of these values at the vertex between them. Namely, let F(i+1/2,j+1/2,k+1/2)
      ! be computed at a cell vertex. Then the following hold:
      !
      ! F(i+1/2,j+1/2,k+1/2) = 1/8 * ( F(i,j,k) + F(i,j+1,k) + F(i,j,k+1) + F(i,j+1,k+1)
      !                                F(i+1,j,k) + F(i+1,j+1,k) + F(i+1,j,k+1) + F(i+1,j+1,k+1) )

      Centers2Vertex = 0.125d0 * ( C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 )

    end function InterpCenters2Vertex

    function InterpVertices2Centers(V1,V2,V3,V4,V5,V6,V7,V8) result(Vertices2Center)
      implicit none
      CCTK_REAL, intent(IN) :: V1, V2, V3, V4, V5, V6, V7, V8
      CCTK_REAL :: Vertices2Center

      ! Given values computed at vertices of a cell, this function computes the linear
      ! interpolation of these values at the center of the cell. Namely, let F(i-1/2,j-1/2,k-1/2),
      ! F(i+1/2,j-1/2,k-1/2), F(i+1/2,j+1/2,k-1/2), F(i-1/2,j+1/2,k-1/2), F(i-1/2,j-1/2,k+1/2),
      ! F(i+1/2,j-1/2,k+1/2), F(i+1/2,j+1/2,k+1/2), F(i-1/2,j+1/2,k+1/2) define the vertices of
      ! cell (i,j,k). Then the following hold:
      !
      ! Lower Face :
      ! F(i,j,k-1/2) = 1/4 * ( F(i-1/2,j-1/2,k-1/2) + F(i+1/2,j-1/2,k-1/2) +
      !                        F(i+1/2,j+1/2,k-1/2) + F(i-1/2,j+1/2,k-1/2) )
      !
      ! Upper Face:
      ! F(i,j,k+1/2) = 1/4 * ( F(i-1/2,j-1/2,k+1/2) + F(i+1/2,j-1/2,k+1/2) +
      !                        F(i+1/2,j+1/2,k+1/2) + F(i-1/2,j+1/2,k+1/2) )
      !
      ! At cell center:
      ! F(i,j,k) = 1/2 * ( F(i,j,k-1/2) + F(i,j,k+1/2) ) =
      !          = 1/8 * ( F(i-1/2,j-1/2,k-1/2) + F(i+1/2,j-1/2,k-1/2) +
      !                    F(i+1/2,j+1/2,k-1/2) + F(i-1/2,j+1/2,k-1/2) +
      !                    F(i-1/2,j-1/2,k+1/2) + F(i+1/2,j-1/2,k+1/2) +
      !                    F(i+1/2,j+1/2,k+1/2) + F(i-1/2,j+1/2,k+1/2) )

      Vertices2Center = 0.125d0 * ( V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 )

    end function InterpVertices2Centers

    function InterpEdges2Center(E1,E2,E3,E4) result(Edges2Center)
      implicit none
      CCTK_REAL, intent(IN) :: E1, E2, E3, E4
      CCTK_REAL :: Edges2Center

      ! Given values computed at vertices of a cell, this function computes the linear
      ! interpolation of these values at the center of the cell. E. g. let F(i,j-1/2,k-1/2),
      ! F(i,j+1/2,k-1/2), F(i,j-1/2,k+1/2), F(i,j+1/2,k+1/2) define the edges of
      ! cell (i,j,k). Then the following hold:
      !
      ! Lower Face :
      ! F(i,j,k-1/2) = 1/2 * ( F(i,j-1/2,k-1/2) + F(i,j+1/2,k-1/2) )
      !
      ! Upper Face:
      ! F(i,j,k+1/2) = 1/2 * ( F(i,j-1/2,k+1/2) + F(i,j+1/2,k+1/2) )
      !
      ! At cell center:
      ! F(i,j,k) = 1/2 * ( F(i,j,k-1/2) + F(i,j,k+1/2) ) =
      !          = 1/4 * ( F(i,j-1/2,k-1/2) + F(i,j+1/2,k-1/2) +
      !                    F(i,j-1/2,k+1/2) + F(i,j+1/2,k+1/2) )

      Edges2Center = 0.25d0 * ( E1 + E2 + E3 + E4 )

    end function InterpEdges2Center

end subroutine SpritzMHDCT

