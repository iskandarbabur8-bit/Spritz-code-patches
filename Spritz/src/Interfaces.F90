#include "cctk.h"

module Interfaces

  interface
     subroutine num_x_flux_MHD(dens,sx,sy,sz,tau,Bx,By,Bz,densf,sxf,syf,szf,tauf,Bxf,Byf,Bzf,&
          velx,vely,velz,avg_betax,avg_betay,avg_betaz,avg_alp,avg_det,w_lorentz,ptot,b0,bxlow,bylow,bzlow)
       implicit none
       CCTK_REAL :: dens,sx,sy,sz,tau,Bx,By,Bz
       CCTK_REAL :: densf,sxf,syf,szf,tauf,Bxf,Byf,Bzf
       CCTK_REAL :: avg_betax,avg_betay,avg_betaz,avg_alp,avg_det
       CCTK_REAL :: velx,vely,velz,ptot,b0,bxlow,bylow,bzlow
       CCTK_REAL :: w_lorentz
     end subroutine num_x_flux_MHD
  end interface

  interface
     subroutine cubic(c,W,check)
       implicit none
       CCTK_REAL,dimension(0:2),intent(IN)::c
       CCTK_REAL,dimension(1:3),intent(OUT)::W
       logical(KIND(.true.)),intent(OUT),optional::check
     end subroutine cubic
  end interface

  interface
     subroutine quartic(a,rtr,x,y,z)
       implicit none
       CCTK_REAL,dimension(0:4),intent(IN)::a
       CCTK_REAL,dimension(4),intent(OUT)::rtr
       CCTK_REAL,intent(IN)::x,y,z
     end subroutine quartic
  end interface

  interface
     subroutine quartic_new(a,rtr,x,y,z)
       implicit none
       CCTK_REAL,dimension(0:4),intent(INOUT)::a
       CCTK_REAL,dimension(4),intent(OUT)::rtr
       CCTK_REAL,intent(IN)::x,y,z
     end subroutine quartic_new
  end interface

  interface
     function cuberoot(num)
       implicit none
       CCTK_REAL::num,cuberoot
     end function cuberoot
  end interface

  interface
     subroutine eigenvalues_MHD(alp,b2,b0,cs2,vx,Bx,bsx,betax,uxx,rho,h,wlorentz,lam,x,y,z)
       implicit none
       CCTK_REAL :: alp,b2,b0,cs2,vx,Bx,bsx,betax,uxx,rho,h,wlorentz
       CCTK_REAL,dimension(4) :: lam
       CCTK_REAL,intent(IN) :: x,y,z
     end subroutine eigenvalues_MHD
  end interface

  interface
     subroutine jacobian_mhd(handle,udens,s2,utau,B2,BS,unk,fvec,fjac)
       implicit none
       CCTK_INT :: handle
       CCTK_REAL :: udens,s2,utau,B2,BS
       CCTK_REAL,dimension(2) :: unk,fvec
       CCTK_REAL,dimension(2,2) :: fjac
     end subroutine jacobian_mhd
  end interface

  interface
     subroutine Con2Prim_pt_hot_Palenzuela(handle, dens, sx, sy, sz, dBx, dBy, dBz, &
          tau, rho, velx, vely, velz, &
          eps, press, w_lorentz, ye, temp, entropy, uxx, uxy, uxz, uyy, &
          uyz, uzz, det, x, y, z, r, Spritz_rho_min, &
          Bx, By, Bz, B2, Bxlow, Bylow, Bzlow, &
          Spritz_reflevel)
       implicit none
       CCTK_REAL:: dens, sx, sy, sz, dBx, dBy, dBz, tau, rho, velx, vely, velz, eps, &
            press, uxx, uxy, uxz, uyy, uyz, uzz, det, w_lorentz, x, &
            y, z, r, Spritz_rho_min, Bx, By, Bz, B2, Bxlow, Bylow, Bzlow, &
            udens, usx, usy, usz, utau, s2, BS, ye, temp, entropy
       CCTK_INT:: i, handle, Spritz_reflevel, error_level
     end subroutine Con2Prim_pt_hot_Palenzuela
  end interface

  interface
    subroutine BrentMethod ( a, b, cons_q, cons_r, cons_s, cons_t, cons_udens, ye, temp, machep, t, x)
      implicit none
      CCTK_REAL, intent(in)  :: a, b, cons_q, cons_r, cons_s, cons_t, cons_udens, ye, machep, t
      CCTK_REAL, intent(out) :: x
      CCTK_REAL, intent(inout) :: temp
    end subroutine
  end interface

  interface
    subroutine f_Palenzuela_value ( x, q, r, s, t, udens, ye, temp, fval )
      implicit none
      CCTK_REAL, intent(in)    :: x, q, r, s, t, udens, ye
      CCTK_REAL, intent(inout) :: temp
      CCTK_REAL, intent(out)   :: fval
    end subroutine
  end interface

  interface
     subroutine Con2Prim_pt_cold_ConstrainT(handle, dens, sx, sy, sz, dBx, dBy, dBz, &
          tau, rho, velx, vely, velz, &
          eps, press, w_lorentz, ye, temp, entropy, uxx, uxy, uxz, uyy, &
          uyz, uzz, det, x, y, z, r, Spritz_rho_min, &
          Bx, By, Bz, B2, Bxlow, Bylow, Bzlow, &
          Spritz_reflevel)
       implicit none
       CCTK_REAL:: dens, sx, sy, sz, dBx, dBy, dBz, tau, rho, velx, vely, velz, eps, &
            press, uxx, uxy, uxz, uyy, uyz, uzz, det, w_lorentz, x, &
            y, z, r, Spritz_rho_min, Bx, By, Bz, B2, Bxlow, Bylow, Bzlow, &
            udens, usx, usy, usz, utau, s2, BS, ye, temp, entropy
       CCTK_INT:: i, handle, Spritz_reflevel, error_level
     end subroutine
  end interface

  interface
    subroutine BrentMethodCold ( a, b, udens, s2, B2, BS, ye, temp, machep, t, x)
      implicit none
      CCTK_REAL, intent(in)  :: a, b, udens, s2, B2, BS, ye, temp, machep, t
      CCTK_REAL, intent(out) :: x
    end subroutine
  end interface

  interface
    subroutine f_cold_eos_mhd_value ( w_lorentz, udens, s2, B2, BS, ye, temp, fval )
      implicit none
      CCTK_REAL, intent(inout)    :: w_lorentz
      CCTK_REAL, intent(in)       :: udens, s2, B2, BS, ye, temp
      CCTK_REAL, intent(out)      :: fval
    end subroutine
  end interface

  interface
     subroutine jacobian_mhd_1eq(handle,udens,s2,utau,B2,BS,unk,fvec,fjac)
       implicit none
       CCTK_INT :: handle
       CCTK_REAL :: udens,s2,utau,B2,BS,unk,fvec,fjac
     end subroutine jacobian_mhd_1eq
  end interface

  interface
     subroutine Con2Prim_pt_mhd(handle, dens, sx, sy, sz, dBx, dBy, dBz, &
          tau, rho, velx, vely, &
          velz, eps, press, w_lorentz, uxx, uxy, uxz, uyy, &
          uyz, uzz, det, x, y, z, r, epsnegative,Spritz_rho_min, &
          Bx, By, Bz, B2, Bxlow, Bylow, Bzlow,Spritz_polytrope_handle,pmin,Spritz_reflevel,c2perror,outside_ah, mhd_c2p_mask_local)
       implicit none
       CCTK_REAL:: dens, sx, sy, sz, dBx, dBy, dBz, tau, rho, velx, vely, velz, eps, &
            press, uxx, uxy, uxz, uyy, uyz, uzz, det, w_lorentz, x, &
            y, z, r, Spritz_rho_min, Bx, By, Bz, B2, Bxlow, Bylow, Bzlow, pmin
       CCTK_INT:: handle, Spritz_polytrope_handle, Spritz_reflevel
       CCTK_REAL, optional:: mhd_c2p_mask_local
       logical epsnegative, c2perror, outside_ah
     end subroutine Con2Prim_pt_mhd
  end interface

  interface
     subroutine Con2Prim_ptPolytype_mhd(handle, dens, sx, sy, sz, &
          dBx, dBy, dBz, tau, rho, &
          velx, vely, velz, eps, press, w_lorentz, uxx, uxy, uxz, uyy, &
          uyz, uzz, det, x, y, z, r, Spritz_rho_min, &
          Bx, By, Bz, B2, Bxlow, Bylow, Bzlow, Spritz_reflevel, c2perror, mhd_c2p_mask_local)
       implicit none
       CCTK_REAL:: dens, sx, sy, sz, dBx, dBy, dBz, tau, rho, velx, vely, velz, eps, &
            press, uxx, uxy, uxz, uyy, uyz, uzz, det, w_lorentz, x, &
            y, z, r, Spritz_rho_min, Bx, By, Bz, B2, Bxlow, Bylow, Bzlow
       CCTK_INT :: handle, Spritz_reflevel
       CCTK_REAL, optional:: mhd_c2p_mask_local
       logical, optional::c2perror
     end subroutine Con2Prim_ptPolytype_mhd
  end interface

  interface
     subroutine jacobian_mhd_polytype(handle,udens,s2,B2,BS,unk,fvec,fjac)
       implicit none
       CCTK_INT :: handle
       CCTK_REAL :: udens,s2,utau,B2,BS
       CCTK_REAL :: unk,fvec
       CCTK_REAL :: fjac
     end subroutine jacobian_mhd_polytype
  end interface

  interface
     SUBROUTINE zrhqr(a,rtr,rti)
       IMPLICIT NONE
       CCTK_REAL, DIMENSION(:), INTENT(IN) :: a
       CCTK_REAL, DIMENSION(:), INTENT(OUT) :: rtr,rti
     end SUBROUTINE zrhqr
  end interface

  interface
     subroutine con2prim_rtsafe(handle, dens, sx, sy, sz, tau, rho, velx, vely, &
          velz, eps, press, w_lorentz, uxx, uxy, uxz, uyy, &
          uyz, uzz, det, x, y, z, r, epsnegative,Spritz_rho_min, &
          Bx, By, Bz, B2, Bxlow, Bylow, Bzlow, Spritz_reflevel)
       implicit none
       CCTK_REAL:: dens, sx, sy, sz, tau, rho, velx, vely, velz, eps, &
            press, uxx, uxy, uxz, uyy, uyz, uzz, det, w_lorentz, x, &
            y, z, r, Spritz_rho_min, Bx, By, Bz, B2, Bxlow, Bylow, Bzlow, &
            udens, usx, usy, usz, utau, s2, BS, pold, rhoold, w_lorentzold, w, ww2
       CCTK_INT:: handle, Spritz_reflevel
       character(len=300) warnline
       logical epsnegative
     end subroutine con2prim_rtsafe
  end interface

  interface
     SUBROUTINE funcd(w_lorentz,gamma1,B2,udens,utau,BS,S2,fval,fderiv,err)
       IMPLICIT NONE
       CCTK_REAL, INTENT(IN) :: w_lorentz,gamma1,B2,udens,utau,BS,S2
       CCTK_REAL, INTENT(OUT) :: fval,fderiv
       logical,INTENT(OUT),OPTIONAL:: err
     end SUBROUTINE funcd
  end interface

  interface
     subroutine Pfunc(w_lorentz,P,Pderiv,gamma1,B2,udens,utau,BS)
       implicit none
       CCTK_REAL,intent(IN)::w_lorentz
       CCTK_REAL,intent(OUT)::P,Pderiv
       CCTK_REAL, intent(IN)::gamma1,B2,udens,utau,BS
     end subroutine Pfunc
  end interface

!!$  interface
!!$     subroutine prim2con_mhd(handle, gxx, gxy, gxz, gyy, gyz, gzz, det, ddens, &
!!$          dsx, dsy, dsz, dtau , drho, dvelx, dvely, dvelz, deps, dpress, w, &
!!$          Bx, By, Bz, dBx, dBy, dBz, alp, betax, betay, betaz)
!!$       implicit none
!!$       CCTK_REAL :: gxx, gxy, gxz, gyy, gyz, gzz, det
!!$       CCTK_REAL :: ddens, dsx, dsy, dsz, dtau, drho, dvelx, dvely, dvelz,&
!!$            deps, dpress, w, vlowx, vlowy, vlowz, dBx, dBy, dBz, &
!!$            Blowx, Blowy, Blowz, bst, bslowx, bslowy, bslowz, bs2, &
!!$            Bx, By, Bz, alp, betax, betay, betaz
!!$       CCTK_INT :: handle
!!$     end subroutine prim2con_mhd
!!$  end interface

!!$  interface
!!$     subroutine prim2conpolytype_mhd(handle, gxx, gxy, gxz, gyy, gyz, gzz, det, ddens, &
!!$          dsx, dsy, dsz, dtau , drho, dvelx, dvely, dvelz, deps, dpress, w, &
!!$          Bx, By, Bz, dBx, dBy, dBz, alp, betax, betay, betaz)
!!$       implicit none
!!$       CCTK_REAL :: gxx, gxy, gxz, gyy, gyz, gzz, det
!!$       CCTK_REAL :: ddens, dsx, dsy, dsz, dtau, drho, dvelx, dvely, dvelz,&
!!$            deps, dpress, w, vlowx, vlowy, vlowz, dBx, dBy, dBz,&
!!$            Blowx, Blowy, Blowz, bst, bslowx, bslowy, bslowz, bs2, &
!!$            Bx, By, Bz, alp, betax, betay, betaz, betaxlow, betaylow, betazlow
!!$       CCTK_INT :: handle
!!$     end subroutine prim2conpolytype_mhd
!!$  end interface

  interface
     subroutine SimplePPM_1d_MHD(handle,poly,&
          nx,dx,rho,velx,vely,velz,eps,press,Bx,By,Bz,rhominus,&
          velxminus,velyminus,velzminus,epsminus,Bxminus,Byminus,Bzminus,&
          rhoplus,velxplus,velyplus,&
          velzplus,epsplus,Bxplus,Byplus,Bzplus,trivial_rp, &
          gxx, gxy, gxz, gyy, gyz, gzz, psi4, beta, alp, w_lorentz, &
          dir, ni, nj, nrx, nry, nrz, Spritz_Cs2_min, &
          psi_mhd, psi_mhd_minus, psi_mhd_plus )
       implicit none
       CCTK_INT :: handle,poly,nx
       CCTK_REAL :: dx
       CCTK_REAL, dimension(nx) :: rho,velx,vely,velz,eps,Bx,By,Bz
       CCTK_REAL, dimension(nx) :: rhominus,velxminus,velyminus,velzminus,epsminus,Bxminus,Byminus,Bzminus
       CCTK_REAL, dimension(nx) :: rhoplus,velxplus,velyplus,velzplus,epsplus,Bxplus,Byplus,Bzplus
       CCTK_REAL, dimension(nx) :: press
       logical, dimension(nx) :: trivial_rp
       CCTK_REAL, dimension(nx) :: gxx, gxy, gxz, gyy, gyz, gzz, &
            psi4, beta, alp, w_lorentz
       CCTK_INT :: dir, ni, nj, nrx, nry, nrz
       CCTK_REAL, dimension(nrx, nry, nrz) :: ev_l, ev_r, xw
       CCTK_REAL :: Spritz_Cs2_min
       CCTK_REAL, dimension( nx),OPTIONAL :: psi_mhd, psi_mhd_minus, psi_mhd_plus
     end subroutine SimplePPM_1d_MHD
  end interface
  
  interface
  		subroutine wenoz_reconstruct_1D(N, v, vminus, vplus, trivial_rp)
  			implicit none
	
			CCTK_INT :: N, i
			CCTK_REAL :: tau
			CCTK_REAL, dimension(N) :: v, vminus, vplus
			CCTK_REAL, dimension(2,3) :: alphaZ, omegaZ !! 1->plus, 2->minus
			CCTK_REAL, dimension(3) :: beta
			logical, dimension(N) :: trivial_rp
			CCTK_INT :: HO_stc
  		end subroutine wenoz_reconstruct_1D
  end interface



  interface
     subroutine Con2Prim_pt_mhd_1eq(handle, dens, sx, sy, sz, dBx, dBy, dBz, &
          tau, rho, velx, vely, &
          velz, eps, press, w_lorentz, uxx, uxy, uxz, uyy, &
          uyz, uzz, det, x, y, z, r, epsnegative,Spritz_rho_min, &
          Bx, By, Bz, B2, Bxlow, Bylow, Bzlow, Spritz_polytrope_handle, pmin, &
          Spritz_reflevel, c2perror, outside_ah, mhd_c2p_mask_local)
       implicit none
       CCTK_REAL :: dens, sx, sy, sz, dBx, dBy, dBz, tau, rho, velx, vely, velz, eps, &
       press, uxx, uxy, uxz, uyy, uyz, uzz, det, w_lorentz, x, &
       y, z, r, Spritz_rho_min, Bx, By, Bz, B2, Bxlow, Bylow, Bzlow, pmin
       CCTK_INT :: Spritz_polytrope_handle, Spritz_reflevel, handle
       CCTK_REAL, optional :: mhd_c2p_mask_local
       logical :: epsnegative, c2perror, outside_ah
     end subroutine Con2Prim_pt_mhd_1eq
  end interface

end module Interfaces
