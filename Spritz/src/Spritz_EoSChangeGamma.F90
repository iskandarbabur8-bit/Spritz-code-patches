/*@@
@file      Spritz_EOSResetHydro.F90
@date      Sat Jan 26 01:36:57 2002
@author     Ian Hawke
@desc 
   This routine will reset the specific internal energy using the polytropic
   EOS that will be used at evolution. This is wanted if the EoS changes
   between setting up the initial data and evolving
  @enddesc 
  @@*/
  
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
 
#define velx(i,j,k) vel(i,j,k,1)
#define vely(i,j,k) vel(i,j,k,2)
#define velz(i,j,k) vel(i,j,k,3)
#define sx(i,j,k) scon(i,j,k,1)
#define sy(i,j,k) scon(i,j,k,2)
#define sz(i,j,k) scon(i,j,k,3)

 /*@@
   @routine    Spritz_EOSResetHydro
   @date       Sat Jan 26 01:38:12 2002
   @author     Ian Hawke
   @desc 
   see above
   @enddesc 
   @calls     
   @calledby   
   @history 
   
   @endhistory 

@@*/


subroutine Spritz_EoSChangeGamma(CCTK_ARGUMENTS)

  USE EOS_Polytrope_Scalars

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: i, j, k
  CCTK_REAL :: det
  
  CCTK_REAL :: local_gamma
  
!#include "CactusEOS/EOS_Base/src/EOS_Base.inc"

  !!!!!!
  ! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1),one(1)
  
  n = 1;keytemp = 0;anyerr = 0;keyerr = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
  one = 1.0d0
  ! end EOS Omni vars
  !!!!!

!!$  Set up the fluid constants

  !local_Gamma = 1.d0 + EOS_Pressure(Spritz_polytrope_handle, 1.d0, 1.d0) / &
  !     EOS_SpecificIntEnergy(Spritz_polytrope_handle, 1.d0, 1.d0)
  
  !!!!!!
  call EOS_Omni_press(Spritz_polytrope_handle,keytemp,Eos_rf_prec,n,&
                 one,one,xtemp,xye,xpress,keyerr,anyerr)
  call EOS_Omni_EpsFromPress(Spritz_polytrope_handle,keytemp,Eos_rf_prec,n,&
			one,xeps,xtemp,xye,one,xeps,keyerr,anyerr)
  local_gamma = 1.0d0 + xpress(1)/xeps(1)
  !!!!!!

!!$  Change the pressure and specific internal energy

  press = p_geom_factor * eos_k_cgs * &
       (rho * rho_geom_factor_inv)**local_Gamma 
  eps = press / (rho * (local_Gamma - 1.d0))

!!$  Get the conserved variables. Hardwired polytrope EoS!!!
!!$  Note that this call also sets pressure and eps
    
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
        
        call SpatialDeterminant(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),&
             gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det)

        if (Spritz_mhd_on .ne. 0) then
           call prim2conpolytype(Spritz_polytrope_handle,gxx(i,j,k),gxy(i,j,k),&
                gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
                det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
                tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
                eps(i,j,k),press(i,j,k),w_lorentz(i,j,k),Bx(i,j,k),By(i,j,k),&
                Bz(i,j,k))
        else
           call prim2conpolytype(Spritz_polytrope_handle,gxx(i,j,k),gxy(i,j,k),&
                gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
                det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
                tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
                eps(i,j,k),press(i,j,k),w_lorentz(i,j,k))
        end if

      end do
    end do
  end do


end subroutine Spritz_EoSChangeGamma

 /*@@
   @routine    Spritz_EoSChangeK
   @date       Mon Oct 20 12:56:14 2003
   @author     Ian Hawke
   @desc 
   Reset the hydro variables when K is changed.
   Unlike the routine above, this actually gives a solution to
   the constraints.

   Only two cases are given as the general case is transcendental.
   We find this by holding rho * enthalpy fixed and assuming a
   polytropic EOS.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Spritz_EoSChangeK(CCTK_ARGUMENTS)

  USE EOS_Polytrope_Scalars

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: i, j, k
  CCTK_REAL :: det
  
  CCTK_REAL :: local_gamma, local_k
  
  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: Q

!#include "CactusEOS/EOS_Base/src/EOS_Base.inc"

  !!!!!!
  ! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1),one(1)
  
  n = 1;keytemp = 0;anyerr = 0;keyerr = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
  one = 1.0d0
  ! end EOS Omni vars
  !!!!!

!!$  Set up the fluid constants

  !local_Gamma = 1.d0 + EOS_Pressure(Spritz_polytrope_handle, 1.d0, 1.d0) / &
  !     EOS_SpecificIntEnergy(Spritz_polytrope_handle, 1.d0, 1.d0)

  !local_K = EOS_Pressure(Spritz_polytrope_handle, 1.d0, 1.d0)

  !!!!!!
  call EOS_Omni_press(Spritz_polytrope_handle,keytemp,Eos_rf_prec,n,&
                 one,one,xtemp,xye,xpress,keyerr,anyerr)
  call EOS_Omni_EpsFromPress(Spritz_polytrope_handle,keytemp,Eos_rf_prec,n,&
			one,xeps,xtemp,xye,one,xeps,keyerr,anyerr)
  local_gamma = 1.0d0 + xpress(1)/xeps(1)
  local_k =xpress(1) 
  !!!!!!

  if (abs(local_Gamma - 2.d0) < 1.d-10) then

    rho = -0.5d0/local_k+sqrt(0.25d0/local_k**2+(rho+initial_k*rho**2)/local_k)

  else if (abs(local_Gamma - 3.d0) < 1.d-10) then

!!$ This code is probably just wrong. We have never used it anyway.
    Q = -9.d0 * local_k**2 * rho * (2.d0 + 3.d0 * initial_k * rho**2) + &
         sqrt(local_k**3 * (32.d0 + 81.d0 * local_k * rho**2 * &
         (2.d0 + 3.d0 * initial_k * rho**2)**2))
    
    Q = Q**(1.d0/3.d0)
    
    rho = (2**(7.d0/3.d0) * local_k - 2**(2.d0/3.d0) * Q**2) / &
         (6.d0 * local_k * Q)

  else
    call CCTK_WARN(0, "EoSChangeK only knows how to do Gamma=2 or 3!")
  end if
  
  press = local_k * rho**local_gamma
  eps = (local_gamma - 1.d0) * local_k * rho**local_gamma

  call Primitive2ConservativePolyCells(CCTK_ARGUMENTS)

end subroutine Spritz_EoSChangeK



 /*@@
   @routine    Spritz_EoSChangeGammaK_Shibata
   @date       Jan. 2005
   @author     Christian D. Ott
   @desc 
   Reset the hydro variables when K and Gamma are changed.

   This is according to Shibata in astro-ph/0412243 (PRD71 024014) in
   which he switches K and Gamma after initial data setup,
   but keeps the internal energy constant.

   Note: this works only with one refinement level

   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Spritz_EoSChangeGammaK_Shibata(CCTK_ARGUMENTS)

  USE EOS_Polytrope_Scalars

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: i, j, k
  CCTK_REAL :: det
  
  CCTK_REAL :: local_Gamma, local_k, eos_k_initial_cgs

  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: Q

  character(len=100) infoline

!#include "CactusEOS/EOS_Base/src/EOS_Base.inc"
  
  !!!!!!
  ! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1),one(1)
  
  n = 1;keytemp = 0;anyerr = 0;keyerr = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
  one = 1.0d0
  ! end EOS Omni vars
  !!!!!

!!$  Set up the fluid constants

  call CCTK_INFO("Pulling the rug...by changing K and Gamma")

  !local_Gamma = 1.d0 + EOS_Pressure(Spritz_polytrope_handle, 1.d0, 1.d0) / &
  !     EOS_SpecificIntEnergy(Spritz_polytrope_handle, 1.d0, 1.d0)

  !local_K = EOS_Pressure(Spritz_polytrope_handle, 1.d0, 1.d0)

  !!!!!!
  call EOS_Omni_press(Spritz_polytrope_handle,keytemp,Eos_rf_prec,n,&
                 one,one,xtemp,xye,xpress,keyerr,anyerr)
  call EOS_Omni_EpsFromPress(Spritz_polytrope_handle,keytemp,Eos_rf_prec,n,&
			one,xeps,xtemp,xye,one,xeps,keyerr,anyerr)
  local_gamma = 1.0d0 + xpress(1)/xeps(1)
  local_k =xpress(1) 
  !!!!!!

  eos_k_initial_cgs = initial_k * rho_geom_factor**initial_Gamma / p_geom_factor

  press = (local_Gamma - 1.d0) / (initial_Gamma - 1.0d0 ) * p_geom_factor * eos_k_initial_cgs * &
  	     (rho * rho_geom_factor_inv) ** initial_Gamma

  eps = p_geom_factor * eos_k_initial_cgs * & 
     (rho * rho_geom_factor_inv) ** initial_Gamma / &
     (rho * (initial_Gamma - 1.0d0))

  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
        
        call SpatialDeterminant(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),&
             gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det)

        IF ( Spritz_mhd_on /= 0 ) THEN
           CALL prim2con(                                                   &
                & Spritz_eos_handle,                                        &
                & gxx(i,j,k), gxy(i,j,k), gxz(i,j,k),                       &
                & gyy(i,j,k), gyz(i,j,k), gzz(i,j,k), det,                  &
                & dens(i,j,k), sx(i,j,k), sy(i,j,k), sz(i,j,k), tau(i,j,k), &
                & rho(i,j,k), velx(i,j,k), vely(i,j,k), velz(i,j,k),        &
                & eps(i,j,k), press(i,j,k), w_lorentz(i,j,k),               &
                & Bx(i,j,k), By(i,j,k), Bz(i,j,k) )
        ELSE
           CALL prim2con(                                                   &
                & Spritz_eos_handle,                                        &
                & gxx(i,j,k), gxy(i,j,k), gxz(i,j,k),                       &
                & gyy(i,j,k), gyz(i,j,k), gzz(i,j,k), det,                  &
                & dens(i,j,k), sx(i,j,k), sy(i,j,k), sz(i,j,k), tau(i,j,k), &
                & rho(i,j,k), velx(i,j,k), vely(i,j,k), velz(i,j,k),        &
                & eps(i,j,k), press(i,j,k), w_lorentz(i,j,k),               &
                & Spritz_reflevel )
        END IF

      end do
    end do
  end do


end subroutine Spritz_EoSChangeGammaK_Shibata
