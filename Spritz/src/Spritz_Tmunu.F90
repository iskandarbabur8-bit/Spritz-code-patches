 /*@@
   @file      Spritz_Tmunu.F90
   @date      Thu Apr 16 19:38:40 2009
   @author    Ian Hawke
   @histpry
    Apr. 2009: Luca Baiotti copied and adapted for the Tmunu-thorn mechanism the original include file
   @desc 
     The calculation of the stress energy tensor.
     The version used here was worked out by Miguel Alcubierre. I
     think it was an extension of the routine from GR3D, written
     by Mark Miller.
     C version added by Ian Hawke.

     Lower components of the stress-energy tensor obtained from
     the hydro variables.  The components are given by:

     T      =  rho h  u   u    +  P  g
      mu nu            mu  nu         mu nu

     where rho is the energy density of the fluid, h the enthalpy
     and P the pressure.  The enthalpy is given in terms of the
     basic variables as:

     h  =  1  +  e  +  P/rho

     with e the internal energy (eps here).

     In the expresion for T_{mu,nu} we also have the four-velocity
     of the fluid given by (v_i is the 3-velocity field):

                                 i
     u  =  W ( - alpha +  v  beta  )
      0                    i

     u  =  W v
      i       i
                                                i  -1/2
     with W the Lorentz factor:   W = ( 1 -  v v  )
                                              i

     and where alpha and beta are the lapse and shift vector.

     Finally, the 4 metric is given by

                   2             i
     g   =  - alpha  + beta  beta
      00                   i

     g   =  beta
      0i        i


     g   =  gamma      (the spatial metric)
      ij        ij


   @enddesc 
 @@*/
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define velx(i,j,k) vel(i,j,k,1)
#define vely(i,j,k) vel(i,j,k,2)
#define velz(i,j,k) vel(i,j,k,3)

 subroutine Spritz_Tmunu(CCTK_ARGUMENTS)

   implicit none
   
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS

   CCTK_REAL velxlow, velylow, velzlow
   CCTK_REAL betaxlow, betaylow, betazlow, beta2
   CCTK_REAL utlow, uxlow, uylow, uzlow
   CCTK_REAL rhoenthalpy
   CCTK_REAL ux,uy,uz,bst,bsx,bsy,bsz,bs2
   CCTK_REAL bstlow,bsxlow,bsylow,bszlow
   CCTK_REAL Spritz_mhd_coeff_value
   CCTK_INT i,j,k
!   CCTK_INT :: type_bits, atmosphere

   !if Tmunu is not stored in memory do nothing

   if (stress_energy_state.eq.0) then
      !$OMP CRITICAL
      call CCTK_WARN(1,"You are not using Tmunu!")
      !$OMP END CRITICAL
      return
   end if

   !Do not compute Tmunu if you don't want to use that in the Einstein equations
   if (Spritz_use_Tmunu.eq.0) then
      eTtt = 0.0d0
      eTtx = 0.0d0
      eTty = 0.0d0
      eTtz = 0.0d0
      eTxx = 0.0d0
      eTxy = 0.0d0
      eTxz = 0.0d0
      eTyy = 0.0d0
      eTyz = 0.0d0
      eTzz = 0.0d0
      return
   end if

   Spritz_mhd_coeff_value = 1.0d0
   
   !$OMP PARALLEL DO PRIVATE(k, j, i, velxlow, velylow, velzlow, betaxlow, &
   !$OMP betaylow, betazlow, beta2, rhoenthalpy, utlow, uxlow, uylow, uzlow, &
   !$OMP ux, uy, uz, bst, bsx, bsy, bsz, bstlow, bsxlow, bsylow, bszlow, bs2)
   do k = 1, cctk_lsh(3)
      do j = 1, cctk_lsh(2)
         do i = 1, cctk_lsh(1)
            
            if (conformal_state .eq. 0) then
               velxlow = gxx(i,j,k)*velx(i,j,k) + gxy(i,j,k)*vely(i,j,k) + gxz(i,j,k)*velz(i,j,k)
               velylow = gxy(i,j,k)*velx(i,j,k) + gyy(i,j,k)*vely(i,j,k) + gyz(i,j,k)*velz(i,j,k)
               velzlow = gxz(i,j,k)*velx(i,j,k) + gyz(i,j,k)*vely(i,j,k) + gzz(i,j,k)*velz(i,j,k)
            else
               velxlow = psi(i,j,k)**4 * (gxx(i,j,k)*velx(i,j,k) + gxy(i,j,k)*vely(i,j,k) + gxz(i,j,k)*velz(i,j,k) )
               velylow = psi(i,j,k)**4 * (gxy(i,j,k)*velx(i,j,k) + gyy(i,j,k)*vely(i,j,k) + gyz(i,j,k)*velz(i,j,k) )
               velzlow = psi(i,j,k)**4 * (gxz(i,j,k)*velx(i,j,k) + gyz(i,j,k)*vely(i,j,k) + gzz(i,j,k)*velz(i,j,k) )
            end if

            !        Calculate lower components and square of shift vector.
            
            if (shift_state .ne. 0) then
               
               if (conformal_state .eq. 0) then
                  betaxlow = gxx(i,j,k)*betax(i,j,k) + gxy(i,j,k)*betay(i,j,k) + gxz(i,j,k)*betaz(i,j,k)
                  betaylow = gxy(i,j,k)*betax(i,j,k) + gyy(i,j,k)*betay(i,j,k) + gyz(i,j,k)*betaz(i,j,k)
                  betazlow = gxz(i,j,k)*betax(i,j,k) + gyz(i,j,k)*betay(i,j,k) + gzz(i,j,k)*betaz(i,j,k)
               else
                  betaxlow = psi(i,j,k)**4 * (gxx(i,j,k)*betax(i,j,k) + gxy(i,j,k)*betay(i,j,k) + gxz(i,j,k)*betaz(i,j,k) )
                  betaylow = psi(i,j,k)**4 * (gxy(i,j,k)*betax(i,j,k) + gyy(i,j,k)*betay(i,j,k) + gyz(i,j,k)*betaz(i,j,k) )
                  betazlow = psi(i,j,k)**4 * (gxz(i,j,k)*betax(i,j,k) + gyz(i,j,k)*betay(i,j,k) + gzz(i,j,k)*betaz(i,j,k) )
               end if

               beta2 = betax(i,j,k)*betaxlow + betay(i,j,k)*betaylow + betaz(i,j,k)*betazlow 
               
            else
               
               betaxlow = 0.0D0
               betaylow = 0.0D0
               betazlow = 0.0D0
               
               beta2 = 0.0D0
               
            end if

            !        Calculate the specific relativistic enthalpy times rho times the
            !        square of the lorentz factor.

            rhoenthalpy = w_lorentz(i,j,k)**2*(rho(i,j,k)*(1.0D0 + eps(i,j,k)) + press(i,j,k))

            !        Calculate lower components of 4-velocity (without the Lorent factor).
            
            utlow = (-alp(i,j,k) + velx(i,j,k)*betaxlow + vely(i,j,k)*betaylow + velz(i,j,k)*betazlow)
            
            uxlow = velxlow
            uylow = velylow
            uzlow = velzlow
            
            
            !     Calculate the upper components of 4-velocity (without the Lorentz factor)
!!$         ut = 1.0d0/alp(i,j,k) !UNUSED
            ux = velx(i,j,k)-betax(i,j,k)/alp(i,j,k)
            uy = vely(i,j,k)-betay(i,j,k)/alp(i,j,k)
            uz = velz(i,j,k)-betaz(i,j,k)/alp(i,j,k)
            
            
            !     Calculate the 4 vector b of the magnetic field
            bst=w_lorentz(i,j,k)*(Bx(i,j,k)*velxlow + By(i,j,k)*velylow + Bz(i,j,k)*velzlow)/alp(i,j,k)
            bsx=(Bx(i,j,k)+alp(i,j,k)*bst*w_lorentz(i,j,k)*ux)/(w_lorentz(i,j,k))
            bsy=(By(i,j,k)+alp(i,j,k)*bst*w_lorentz(i,j,k)*uy)/(w_lorentz(i,j,k))
            bsz=(Bz(i,j,k)+alp(i,j,k)*bst*w_lorentz(i,j,k)*uz)/(w_lorentz(i,j,k))
            
            !     Calculate the lower components of b
            
            if (conformal_state .eq. 0) then
               bstlow = bst*(-alp(i,j,k)**2+beta2)+bsx*betaxlow + bsy*betaylow + bsz*betazlow
               bsxlow = bst*betaxlow+gxx(i,j,k)*bsx + gxy(i,j,k)*bsy + gxz(i,j,k)*bsz
               bsylow = bst*betaylow+gxy(i,j,k)*bsx + gyy(i,j,k)*bsy + gyz(i,j,k)*bsz
               bszlow = bst*betazlow+gxz(i,j,k)*bsx + gyz(i,j,k)*bsy + gzz(i,j,k)*bsz
            else
               bstlow = bst*(-alp(i,j,k)**2+beta2)+bsx*betaxlow + bsy*betaylow + bsz*betazlow
               
               bsxlow = bst*betaxlow+psi(i,j,k)**4*(gxx(i,j,k)*bsx + gxy(i,j,k)*bsy + gxz(i,j,k)*bsz)
               bsylow = bst*betaylow+psi(i,j,k)**4*(gxy(i,j,k)*bsx + gyy(i,j,k)*bsy + gyz(i,j,k)*bsz)
               bszlow = bst*betazlow+psi(i,j,k)**4*(gxz(i,j,k)*bsx + gyz(i,j,k)*bsy + gzz(i,j,k)*bsz)
            end if
          
            !     Calculate b^2 = b^mu b_mu
            bs2 = bst*bstlow+bsx*bsxlow+bsy*bsylow+bsz*bszlow
        

            !        Calculate Tmunu (the lower components!).

            eTtt(i,j,k) = eTtt(i,j,k) + rhoenthalpy*utlow**2 + press(i,j,k)*(beta2 - alp(i,j,k)**2)
            
            eTtx(i,j,k) = eTtx(i,j,k) + rhoenthalpy*utlow*uxlow + press(i,j,k)*betaxlow
            eTty(i,j,k) = eTty(i,j,k) + rhoenthalpy*utlow*uylow + press(i,j,k)*betaylow
            eTtz(i,j,k) = eTtz(i,j,k) + rhoenthalpy*utlow*uzlow + press(i,j,k)*betazlow

            !adding the MHD part
            eTtt(i,j,k) = eTtt(i,j,k) + Spritz_mhd_coeff_value*((w_lorentz(i,j,k)**2*utlow**2 + 0.5*(-alp(i,j,k)**2+beta2))*bs2-bstlow*bstlow)
            eTtx(i,j,k) = eTtx(i,j,k) + Spritz_mhd_coeff_value*((w_lorentz(i,j,k)**2*utlow*uxlow + 0.5*betaxlow)*bs2-bstlow*bsxlow)
            eTty(i,j,k) = eTty(i,j,k) + Spritz_mhd_coeff_value*((w_lorentz(i,j,k)**2*utlow*uylow + 0.5*betaylow)*bs2-bstlow*bsylow)
            eTtz(i,j,k) = eTtz(i,j,k) + Spritz_mhd_coeff_value*((w_lorentz(i,j,k)**2*utlow*uzlow + 0.5*betazlow)*bs2-bstlow*bszlow)
            
       
            if (conformal_state .eq. 0) then
               
               eTxx(i,j,k) = eTxx(i,j,k) + rhoenthalpy*uxlow**2 + press(i,j,k)*gxx(i,j,k)
               eTyy(i,j,k) = eTyy(i,j,k) + rhoenthalpy*uylow**2 + press(i,j,k)*gyy(i,j,k)
               eTzz(i,j,k) = eTzz(i,j,k) + rhoenthalpy*uzlow**2 + press(i,j,k)*gzz(i,j,k)
               
               eTxy(i,j,k) = eTxy(i,j,k) + rhoenthalpy*uxlow*uylow + press(i,j,k)*gxy(i,j,k)
               eTxz(i,j,k) = eTxz(i,j,k) + rhoenthalpy*uxlow*uzlow + press(i,j,k)*gxz(i,j,k)
               eTyz(i,j,k) = eTyz(i,j,k) + rhoenthalpy*uylow*uzlow + press(i,j,k)*gyz(i,j,k)
     
               !adding the MHD part
               eTxx(i,j,k) = eTxx(i,j,k) + Spritz_mhd_coeff_value*((w_lorentz(i,j,k)**2*uxlow**2 + 0.5*gxx(i,j,k))*bs2-bsxlow**2)
               eTyy(i,j,k) = eTyy(i,j,k) + Spritz_mhd_coeff_value*((w_lorentz(i,j,k)**2*uylow**2 + 0.5*gyy(i,j,k))*bs2-bsylow**2)
               eTzz(i,j,k) = eTzz(i,j,k) + Spritz_mhd_coeff_value*((w_lorentz(i,j,k)**2*uzlow**2 + 0.5*gzz(i,j,k))*bs2-bszlow**2)
               
               eTxy(i,j,k) = eTxy(i,j,k) + Spritz_mhd_coeff_value*((w_lorentz(i,j,k)**2*uxlow*uylow+0.5*gxy(i,j,k))*bs2- bsxlow*bsylow)
               eTxz(i,j,k) = eTxz(i,j,k) + Spritz_mhd_coeff_value*((w_lorentz(i,j,k)**2*uxlow*uzlow+0.5*gxz(i,j,k))*bs2- bsxlow*bszlow)
               eTyz(i,j,k) = eTyz(i,j,k) + Spritz_mhd_coeff_value*((w_lorentz(i,j,k)**2*uylow*uzlow+0.5*gyz(i,j,k))*bs2- bsylow*bszlow) 
               
            else
          
               eTxx(i,j,k) = eTxx(i,j,k) + rhoenthalpy*uxlow**2 + psi(i,j,k)**4*press(i,j,k)*gxx(i,j,k)
               eTyy(i,j,k) = eTyy(i,j,k) + rhoenthalpy*uylow**2 + psi(i,j,k)**4*press(i,j,k)*gyy(i,j,k)
               eTzz(i,j,k) = eTzz(i,j,k) + rhoenthalpy*uzlow**2 + psi(i,j,k)**4*press(i,j,k)*gzz(i,j,k)

               eTxy(i,j,k) = eTxy(i,j,k) + rhoenthalpy*uxlow*uylow + psi(i,j,k)**4*press(i,j,k)*gxy(i,j,k)
               eTxz(i,j,k) = eTxz(i,j,k) + rhoenthalpy*uxlow*uzlow + psi(i,j,k)**4*press(i,j,k)*gxz(i,j,k)
               eTyz(i,j,k) = eTyz(i,j,k) + rhoenthalpy*uylow*uzlow + psi(i,j,k)**4*press(i,j,k)*gyz(i,j,k)
               
               !add the MHD part
               eTxx(i,j,k) = eTxx(i,j,k) + Spritz_mhd_coeff_value*((w_lorentz(i,j,k)**2*uxlow**2 + 0.5*psi(i,j,k)**4*gxx(i,j,k))*bs2-bsxlow**2)
               eTyy(i,j,k) = eTyy(i,j,k) + Spritz_mhd_coeff_value*((w_lorentz(i,j,k)**2*uylow**2 + 0.5*psi(i,j,k)**4*gyy(i,j,k))*bs2-bsylow**2)
               eTzz(i,j,k) = eTzz(i,j,k) + Spritz_mhd_coeff_value*((w_lorentz(i,j,k)**2*uzlow**2 + 0.5*psi(i,j,k)**4*gzz(i,j,k))*bs2-bszlow**2)
               
               eTxy(i,j,k) = eTxy(i,j,k) + Spritz_mhd_coeff_value*((w_lorentz(i,j,k)**2*uxlow*uylow+0.5*psi(i,j,k)**4*gxy(i,j,k))*bs2 - bsxlow*bsylow)
               eTxz(i,j,k) = eTxz(i,j,k) + Spritz_mhd_coeff_value*((w_lorentz(i,j,k)**2*uxlow*uzlow+0.5*psi(i,j,k)**4*gxz(i,j,k))*bs2 - bsxlow*bszlow)
               eTyz(i,j,k) = eTyz(i,j,k) + Spritz_mhd_coeff_value*((w_lorentz(i,j,k)**2*uylow*uzlow+0.5*psi(i,j,k)**4*gyz(i,j,k))*bs2 - bsylow*bszlow)
               
            end if

         end do
      end do
   end do
   !$OMP END PARALLEL DO
   return
 
 end subroutine Spritz_Tmunu
