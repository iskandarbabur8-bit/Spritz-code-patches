 /*@@
   @file      Spritz_MHD_Con2Prim.F90
   @date
   @author
   @desc
   The routines for converting conservative to primitive variables in GRMHD.
   @enddesc
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    Conservative2Primitive_mhd
   @date
   @author     Bruno Giacomazzo
   @desc
   Wrapper routine that converts from conserved to primitive variables
   at every grid cell centre.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/

subroutine Conservative2Primitive_mhd(CCTK_ARGUMENTS)

  use Interfaces, only: Con2Prim_pt_mhd, Con2Prim_ptPolytype_mhd, Con2Prim_pt_mhd_1eq

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

!#ifdef _EOS_BASE_INC_
!#undef _EOS_BASE_INC_
!#endif
!#include "EOS_Base.inc"

  integer :: i, j, k, nx, ny, nz
  CCTK_REAL :: uxx, uxy, uxz, uyy, uyz, uzz, det, psi4pt, pmin, B2
  CCTK_REAL :: Bxlow, Bylow, Bzlow, idx, idy, idz, divB2_temp
  logical :: epsnegative,err,c2perror,outside_ah
  character(len=100) warnline
  character(len=200) infoline


  CCTK_INT :: k_lowerbound, k_upperbound, j_lowerbound, j_upperbound, i_lowerbound, i_upperbound
  !the following variables are used to convert the centered derivative in the appropriate one-sided derivative
  ! when computing B from the curl of A
  CCTK_INT :: iplus, iminus, jplus, jminus, kplus, kminus
  CCTK_INT :: iiplus, iiminus, jjplus, jjminus, kkplus, kkminus

  CCTK_REAL :: eos_gamma, b2small_loc

  CCTK_REAL :: w, vlowx, vlowy, vlowz
  CCTK_REAL :: bst, bs2

  !!!!!!
  ! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1),one(1),dummy(1)
  CCTK_REAL :: xeos_gamma(1)
  CCTK_REAL :: Retval

  n = 1;keytemp = 0;anyerr = 0;keyerr = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
  one = 1.0d0
  ! end EOS Omni vars
  !!!!!!

  !!*  In recovery section, Con2Prim shouldn't be executed.
  IF ( Recovery_Flag == 1 ) RETURN

  idx = 1.0d0/CCTK_DELTA_SPACE(1)
  idy = 1.0d0/CCTK_DELTA_SPACE(2)
  idz = 1.0d0/CCTK_DELTA_SPACE(3)


  iplus  = 1
  iminus = 1
  jplus  = 1
  jminus = 1
  kplus  = 1
  kminus = 1


!!$  call CCTK_INFO("Conservative2Primitive_mhd")!!!!!!!!!!!!!!!!!!!!!

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)


  !trick to obtain the value of eos_ideal_fluid_gamma from thorn EOS_Ideal
  !eos_gamma = EOS_Pressure(Spritz_eos_handle,1.0d0,1.0d0)+1.0d0

  !!!!!!
  call EOS_Omni_press(Spritz_polytrope_handle,keytemp,Eos_rf_prec,n,&
                 one,one,xtemp,xye,xeos_gamma,keyerr,anyerr)
  eos_gamma = xeos_gamma(1) + 1.0d0
  !!!!!!

  b2small_loc=0.0d0 !initialize b2small_loc=b^mu*b_mu


  if(Spritz_Avec_Staggered.ne.0) then
     k_lowerbound = 2  !whisky_stencil + 1
     k_upperbound = nz !nz - 1 !nz - whisky_stencil

     j_lowerbound = 2  !whisky_stencil + 1
     j_upperbound = ny !ny - 1 !ny - whisky_stencil

     i_lowerbound = 2  !whisky_stencil + 1
     i_upperbound = nx !nx - 1 !nx - whisky_stencil
  else
     k_lowerbound = 1  !2 !Spritz_stencil + 1
     k_upperbound = nz !nz-1 !nz - Spritz_stencil

     j_lowerbound = 1  !2 !Spritz_stencil + 1
     j_upperbound = ny !ny-1 !ny - Spritz_stencil

     i_lowerbound = 1  !2 !Spritz_stencil + 1
     i_upperbound = nx !nx-1 !nx - Spritz_stencil
  end if

  !set MHD variables to be zero in order to check that they are correctly recomputed
  dBx = 0.0d0
  dBy = 0.0d0
  dBz = 0.0d0
  Bx = 0.0d0
  By = 0.0d0
  Bz = 0.0d0

  !pmin=EOS_Pressure(Spritz_polytrope_handle, Spritz_rho_min, 1.0d0)

  !!!!!!
  xrho = Spritz_rho_min
  xtemp= Spritz_hot_atmo_temp
  xye  = Spritz_hot_atmo_Y_e
  keytemp = 1
  call EOS_Omni_press(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                 xrho,one,xtemp,xye,xpress,keyerr,anyerr)
  keytemp = 0
  pmin = xpress(1)
  !!!!!!

  !$OMP PARALLEL DO PRIVATE(k, j, i, det, &
  !$OMP uxx, uxy, uxz, uyy, uyz, uzz, psi4pt, &
  !$OMP iminus, idx, jminus, idy, kminus, idz, iplus, jplus, kplus, &
  !$OMP iiplus, iiminus, jjplus, jjminus, kkplus, kkminus, &
  !$OMP epsnegative, c2perror, outside_ah, &
  !$OMP divB2_temp, Bxlow, Bylow, Bzlow, B2, err, b2small_loc, warnline, &
  !$OMP keytemp, anyerr, keyerr, xpress, xtemp, xye, xeps, xrho, &
  !$OMP xeos_gamma, Retval, infoline, &
  !$OMP w, vlowx, vlowy, vlowz, bst, bs2)
  do k = k_lowerbound, k_upperbound  !!!1, nz
     do j = j_lowerbound, j_upperbound !!!1, ny
        do i = i_lowerbound, i_upperbound !!!1, nx

           !write(infoline,'(a25,g20.7,a8,g20.7,a8,g20.7)') ' C2P input values: dBx = ',dBx(i,j,k), &
           !                                                ', dBy = ',dBy(i,j,k),', dBz = ',dBz(i,j,k)
           !call CCTK_INFO(infoline)

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


           idx = 1.0d0/CCTK_DELTA_SPACE(1)
           idy = 1.0d0/CCTK_DELTA_SPACE(2)
           idz = 1.0d0/CCTK_DELTA_SPACE(3)

           ! The following block is for setting: !!!
           ! dBx(i,1,k) = dBx(i,2,k)             !!!
           ! dBy(i,j,1) = dBy(i,j,2)             !!!
           ! dBz(1,j,k) = dBz(2,j,k)             !!!
                                                 !!!
           if (i .eq. 1) then                    !!!
              iiminus = 0                        !!!
              iiplus = 1                         !!!
           else                                  !!!
              iiminus = 1                        !!!
              iiplus = 0                         !!!
           end if                                !!!
                                                 !!!
           if (j .eq. 1) then                    !!!
              jjminus = 0                        !!!
              jjplus = 1                         !!!
           else                                  !!!
              jjminus = 1                        !!!
              jjplus = 0                         !!!
           end if                                !!!
                                                 !!!
           if (k .eq. 1) then                    !!!
              kkminus = 0                        !!!
              kkplus = 1                         !!!
           else                                  !!!
              kkminus = 1                        !!!
              kkplus = 0                         !!!
           end if                                !!!
                                                 !!!
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           if(Spritz_Avec_Staggered.ne.0) then
              !dBx is 0.5*(curl(A) at i+1/2 + curl(A) at i-1/2)
              dBx(i,j,k) = 0.5d0 * ( (idy*(Avec_z(i+iiplus,j+jjplus,k+kkplus)-Avec_z(i+iiplus,j-jjminus,k+kkplus)) - &
                                      idz*(Avec_y(i+iiplus,j+jjplus,k+kkplus)-Avec_y(i+iiplus,j+jjplus,k-kkminus))) + &
                                     (idy*(Avec_z(i-iiminus,j+jjplus,k+kkplus)-Avec_z(i-iiminus,j-jjminus,k+kkplus)) - &
                                      idz*(Avec_y(i-iiminus,j+jjplus,k+kkplus)-Avec_y(i-iiminus,j+jjplus,k-kkminus))) )
              !dBy is 0.5*(curl(A) at j+1/2 + curl(A) at j-1/2)
              dBy(i,j,k) = 0.5d0 * ( (idz*(Avec_x(i+iiplus,j+jjplus,k+kkplus)-Avec_x(i+iiplus,j+jjplus,k-kkminus)) - &
                                      idx*(Avec_z(i+iiplus,j+jjplus,k+kkplus)-Avec_z(i-iiminus,j+jjplus,k+kkplus))) + &
                                     (idz*(Avec_x(i+iiplus,j-jjminus,k+kkplus)-Avec_x(i+iiplus,j-jjminus,k-kkminus)) - &
                                      idx*(Avec_z(i+iiplus,j-jjminus,k+kkplus)-Avec_z(i-iiminus,j-jjminus,k+kkplus))) )
              !dBz is 0.5*(curl(A) at k+1/2 + curl(A) at k-1/2)
              dBz(i,j,k) = 0.5d0 * ( (idx*(Avec_y(i+iiplus,j+jjplus,k+kkplus)-Avec_y(i-iiminus,j+jjplus,k+kkplus)) - &
                                      idy*(Avec_x(i+iiplus,j+jjplus,k+kkplus)-Avec_x(i+iiplus,j-jjminus,k+kkplus))) + &
                                     (idx*(Avec_y(i+iiplus,j+jjplus,k-kkminus)-Avec_y(i-iiminus,j+jjplus,k-kkminus)) - &
                                      idy*(Avec_x(i+iiplus,j+jjplus,k-kkminus)-Avec_x(i+iiplus,j-jjminus,k-kkminus))) )

              !if ((x(i,j,k) .eq. 0.0d0) .and. (y(i,j,k) .eq. 0.0d0) .and. (z(i,j,k) .eq. 0.0d0)) then
              !  write(infoline,'(a25,g20.7,a8,g20.7,a8,g20.7)') ' C2P compu values: dBx = ',dBx(i,j,k), &
              !                                                  ', dBy = ',dBy(i,j,k),', dBz = ',dBz(i,j,k)
              !  call CCTK_INFO(infoline)
              !end if

           else
              !reset and initialize the constants before starting
              !(these are private variables in OpenMP, so we need to make sure that
              ! they are correctly initialized inside the do loop. Otherwise you may get segfaults
              ! when using optimization, because of out of boundaries errors)
              ! (declaring them as firstprivate would probably also work)
              idx = 1.0d0/CCTK_DELTA_SPACE(1)
              idy = 1.0d0/CCTK_DELTA_SPACE(2)
              idz = 1.0d0/CCTK_DELTA_SPACE(3)
              iplus  = 1
              iminus = 1
              jplus  = 1
              jminus = 1
              kplus  = 1
              kminus = 1

              !check if we have enough points to use the standard second order centered derivative
              ! or if we need to use a one-sided derivative
              if (i==1) then
                 iminus=0
                 idx = 2.0*idx
              end if
              if (j==1) then
                 jminus=0
                 idy = 2.0*idy
              end if
              if (k==1) then
                 kminus=0
                 idz = 2.0*idz
              end if
              if (i==nx) then
                 iplus=0
                 idx = 2.0*idx
              end if
              if (j==ny) then
                 jplus=0
                 idy = 2.0*idy
              end if
              if (k==nz) then
                 kplus=0
                 idz = 2.0*idz
              end if
              !Now we compute sqrt(gamma)*B from the curl of the vector potential A
              !dBx is curl(A)
              dBx(i,j,k) = 0.5d0*(idy*(Avec_z(i,j+jplus,k)-Avec_z(i,j-jminus,k)) - idz*(Avec_y(i,j,k+kplus)-Avec_y(i,j,k-kminus)))
              !dBy is curl(A)
              dBy(i,j,k) = 0.5d0*(idz*(Avec_x(i,j,k+kplus)-Avec_x(i,j,k-kminus)) - idx*(Avec_z(i+iplus,j,k)-Avec_z(i-iminus,j,k)))
              !dBz is curl(A)
              dBz(i,j,k) = 0.5d0*(idx*(Avec_y(i+iplus,j,k)-Avec_y(i-iminus,j,k)) - idy*(Avec_x(i,j+jplus,k)-Avec_x(i,j-jminus,k)))
           end if

           !write(infoline,'(a25,g20.7,a8,g20.7,a8,g20.7)') ' C2P compu values: dBx = ',dBx(i,j,k), &
           !                                                ', dBy = ',dBy(i,j,k),', dBz = ',dBz(i,j,k)
           !call CCTK_INFO(infoline)

           epsnegative = .false.
           c2perror = .false.
           outside_ah = .false.

           !Initialize the mask to zero
           !and set it to 1 when c2p fails to converge or produces NaNs
           mhd_c2p_mask(i,j,k)=0

           !B^j
           Bx(i,j,k) = dBx(i,j,k)/sqrt(det)
           By(i,j,k) = dBy(i,j,k)/sqrt(det)
           Bz(i,j,k) = dBz(i,j,k)/sqrt(det)


           !if (Spritz_mhd_HypDiv_on .ne. 0) psi_mhd(i,j,k) = dpsi_mhd(i,j,k)/sqrt(det)

           !B_j
           Bxlow = Bx(i,j,k)*gxx(i,j,k)+By(i,j,k)*gxy(i,j,k)+Bz(i,j,k)*gxz(i,j,k)
           Bylow = Bx(i,j,k)*gxy(i,j,k)+By(i,j,k)*gyy(i,j,k)+Bz(i,j,k)*gyz(i,j,k)
           Bzlow = Bx(i,j,k)*gxz(i,j,k)+By(i,j,k)*gyz(i,j,k)+Bz(i,j,k)*gzz(i,j,k)

           !B^j B_j
           B2 = gxx(i,j,k) * Bx(i,j,k)**2 + &
                gyy(i,j,k) * By(i,j,k)**2 + &
                gzz(i,j,k) * Bz(i,j,k)**2 + &
                2.d0 * gxy(i,j,k) * Bx(i,j,k) * By(i,j,k) + &
                2.d0 * gxz(i,j,k) * Bx(i,j,k) * Bz(i,j,k) + &
                2.d0 * gyz(i,j,k) * By(i,j,k) * Bz(i,j,k)


          if(evolve_Y_e.ne.0) then
             Y_e(i,j,k) = max(min(Y_e_con(i,j,k) / dens(i,j,k),Spritz_Y_e_max),&
                Spritz_Y_e_min)
          endif

           !!CHECK ON THE ATMO DONE BEFORE THE C2P Netwon-Raphson routine
           !!THIS IS DONE ONLY FOR COMPUTATIONAL EFFICIENCY
           if(dens(i,j,k)<=(sqrt(det)*(Spritz_rho_min*(1.d0+Spritz_atmo_tolerance))).or.(alp(i,j,k)<Spritz_alp_BH)) then
              !atmosphere
              rho(i,j,k) = Spritz_rho_min
              dens(i,j,k) = sqrt(det)*rho(i,j,k)
              !Since I don't know what to do in this case for eps and press
              !I set them using a polytropic EOS which means that the atmosphere is always cold
              !press(i,j,k) = EOS_Pressure(Spritz_polytrope_handle, rho(i,j,k), eps(i,j,k))
              !eps(i,j,k) = EOS_SpecificIntEnergy(Spritz_polytrope_handle,rho(i,j,k),press(i,j,k))

              !!!!!!
              if(GetTemp_FromThermo.ne.0) then
              ! set hot atmosphere values
                temperature(i,j,k) = Spritz_hot_atmo_temp
                y_e_con(i,j,k) = y_e(i,j,k) * dens(i,j,k)
                keytemp = 1 ! don't evolve temperature
                call EOS_Omni_EOS_short(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                    rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
                    press(i,j,k),entropy(i,j,k),dummy,dummy,dummy,dummy,dummy,&
                    keyerr,anyerr) !this just changes press, eps and entropy
                keytemp = 0
              else
                keytemp = 1 !don't evolve temperature, we're setting up atmosphere
                call EOS_Omni_press(Spritz_polytrope_handle,keytemp,Eos_rf_prec,n,&
                   rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),keyerr,anyerr)
                call EOS_Omni_EpsFromPress(Spritz_polytrope_handle,keytemp,Eos_rf_prec,n,&
                   rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),eps(i,j,k),keyerr,anyerr)
                keytemp = 0
              endif

              scon(i,j,k,1) = 0.0d0
              scon(i,j,k,2) = 0.0d0
              scon(i,j,k,3) = 0.0d0
              vel(i,j,k,1) = 0.0d0
              vel(i,j,k,2) = 0.0d0
              vel(i,j,k,3) = 0.0d0
              w_lorentz(i,j,k) = 1.0d0
              !idem here
              tau(i,j,k)  = sqrt(det) * ((rho(i,j,k)+rho(i,j,k)*eps(i,j,k) + &
                   press(i,j,k)+B2)*w_lorentz(i,j,k)*w_lorentz(i,j,k) - (press(i,j,k)+0.5d0*B2)) - dens(i,j,k)
              cycle
           end if !ended setting up atmo values

           if (CCTK_EQUALS(con2prim_mhd_method,"3eqs")) then
              !If the conserved density is far away from the atmo
              ! then check if this point is at least inside an AH
              ! otherwise stop the code in case of errors
              outside_ah = .true.
              call Con2Prim_pt_mhd(Spritz_eos_handle, dens(i,j,k),scon(i,j,k,1),scon(i,j,k,2), &
                   scon(i,j,k,3),dBx(i,j,k),dBy(i,j,k),dBz(i,j,k),&
                   tau(i,j,k),rho(i,j,k),vel(i,j,k,1),vel(i,j,k,2), &
                   vel(i,j,k,3),eps(i,j,k),press(i,j,k),w_lorentz(i,j,k), &
                   uxx,uxy,uxz,uyy,uyz,uzz,det,x(i,j,k),y(i,j,k), &
                   z(i,j,k),r(i,j,k),epsnegative,Spritz_rho_min, &
                   Bx(i,j,k),By(i,j,k),Bz(i,j,k),B2,Bxlow,Bylow,Bzlow,Spritz_polytrope_handle,pmin, &
                   Spritz_reflevel, c2perror, outside_ah,mhd_c2p_mask(i,j,k))
           elseif (CCTK_EQUALS(con2prim_mhd_method,"Noble")) then
              !write(infoline,'(a43)') ' Calling Noble routine...here input values '
              !call CCTK_INFO(infoline)
              !write(infoline,'(a11,i3,a3,i3,a3,i3,a2)') ' point = ( ',i,' , ',j,' , ',k,' )'
              !call CCTK_INFO(infoline)
              !write(infoline,'(a13,g20.7,a3,g20.7,a3,g20.7,a2)') ' x_input = ( ',x(i,j,k),' , ',y(i,j,k),' , ',z(i,j,k),' )'
              !call CCTK_INFO(infoline)
              !write(infoline,'(a13,g20.7,a3,g20.7,a3,g20.7,a2)') ' B_input = ( ',dBx(i,j,k),' , ',dBy(i,j,k),' , ',dBz(i,j,k),' )'
              !call CCTK_INFO(infoline)
              !write(infoline,'(a13,g20.7,a3,g20.7,a3,g20.7,a2)') ' v_input = ( ',vel(i,j,k,1),' , ',vel(i,j,k,2),' , ',vel(i,j,k,3),' )'
              !call CCTK_INFO(infoline)
              call Spritz_Con2Prim_Noble(Spritz_eos_handle,Spritz_reflevel, i, j, k, x(i,j,k), y(i,j,k), &
                   z(i,j,k), keytemp, Eos_rf_prec, Con2Prim_perc_ptol, eos_gamma, dens(i,j,k), &
                   scon(i,j,k,1),scon(i,j,k,2),scon(i,j,k,3), tau(i,j,k), &
                   dBx(i,j,k), dBy(i,j,k), dBz(i,j,k), xye, xtemp, rho(i,j,k), &
                   vel(i,j,k,1),vel(i,j,k,2),vel(i,j,k,3),eps(i,j,k),press(i,j,k),&
                   Bx(i,j,k), By(i,j,k), Bz(i,j,k), b2small_loc, w_lorentz(i,j,k),&
                   gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k), &
                   uxx,uxy,uxz,uyy,uyz,uzz,det, &
                   epsnegative,Spritz_rho_min,Retval)
                   mhd_c2p_mask(i,j,k) = Retval
              if(abs(mhd_c2p_mask(i,j,k)) .gt. 1.0d-8) then
                               c2perror=.true. !mhd_c2p_mask(i,j,k)
                               !$OMP CRITICAL
                               write(warnline,'(a55,g20.7)') 'ATTENTION: Problems in C2P_Nobel - returned Retval = ',Retval
                               call CCTK_WARN(2,warnline)
                               !$OMP END CRITICAL
                               !set mask to failure value as in 3eqs method
                               mhd_c2p_mask(i,j,k) = 1
              end if
           elseif (CCTK_EQUALS(con2prim_mhd_method,"Palenzuela")) then
              if(Constrain_Temp_now.eq.0) then
                call Con2Prim_pt_hot_Palenzuela(Spritz_eos_handle, &
                     dens(i,j,k),scon(i,j,k,1),scon(i,j,k,2), &
                     scon(i,j,k,3),dBx(i,j,k),dBy(i,j,k),dBz(i,j,k),&
                     tau(i,j,k),rho(i,j,k),vel(i,j,k,1),vel(i,j,k,2), &
                     vel(i,j,k,3),eps(i,j,k),press(i,j,k),w_lorentz(i,j,k), &
                     y_e(i,j,k),temperature(i,j,k),entropy(i,j,k), &
                     uxx,uxy,uxz,uyy,uyz,uzz,det,x(i,j,k),y(i,j,k), &
                     z(i,j,k),r(i,j,k),Spritz_rho_min, &
                     Bx(i,j,k),By(i,j,k),Bz(i,j,k),B2,Bxlow,Bylow,Bzlow, &
                     Spritz_reflevel)
              else
                call Con2Prim_pt_cold_ConstrainT(Spritz_eos_handle, &
                     dens(i,j,k),scon(i,j,k,1),scon(i,j,k,2), &
                     scon(i,j,k,3),dBx(i,j,k),dBy(i,j,k),dBz(i,j,k),&
                     tau(i,j,k),rho(i,j,k),vel(i,j,k,1),vel(i,j,k,2), &
                     vel(i,j,k,3),eps(i,j,k),press(i,j,k),w_lorentz(i,j,k), &
                     y_e(i,j,k),temperature(i,j,k),entropy(i,j,k), &
                     uxx,uxy,uxz,uyy,uyz,uzz,det,x(i,j,k),y(i,j,k), &
                     z(i,j,k),r(i,j,k),Spritz_rho_min, &
                     Bx(i,j,k),By(i,j,k),Bz(i,j,k),B2,Bxlow,Bylow,Bzlow, &
                     Spritz_reflevel)
              end if
           else
              !$OMP CRITICAL
              call CCTK_WARN(0,"con2prim_mhd_method not valid")
              !$OMP END CRITICAL
           end if

           if(temperature(i,j,k).gt.Spritz_max_temp .and. CCTK_EQUALS(Spritz_eos_table,"nuc_eos")) then
              !$OMP CRITICAL
              call CCTK_WARN(1,"C2P: Temperature too high")
              write(warnline,"(i8)") cctk_iteration
              call CCTK_WARN(1,warnline)
              write(warnline,"(4i5,1P10E15.6)") Spritz_Reflevel,i,j,k,x(i,j,k),y(i,j,k),z(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,"(1P10E15.6)") dens(i,j,k),scon(i,j,k,1:3),&
                   tau(i,j,k),w_lorentz(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,"(1P10E15.6)") rho(i,j,k),eps(i,j,k),&
                   temperature(i,j,k),Y_e(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,"(A7,i8)") "code: ",keyerr(1)
              call CCTK_WARN(1,warnline)
              write(warnline,"(A10,i5)") "reflevel: ", Spritz_Reflevel
              call CCTK_WARN(1,warnline)
              call CCTK_ERROR("Aborting!!!")
              STOP
              !$OMP END CRITICAL
           endif

           if ((epsnegative.or.c2perror).and..not.CCTK_EQUALS(Spritz_eos_table,"nuc_eos")) then
!!$           call CCTK_WARN(2,"Not able to converge using the Ideal Fluid Con2Prim, Setting the point to atmo")
!!$           !atmosphere
!!$           rho(i,j,k) = Spritz_rho_min
!!$           dens(i,j,k) = sqrt(det)*rho(i,j,k)
!!$           !Since I don't know what to do in this case for eps and press
!!$           !I set them using a polytropic EOS which means that the atmosphere is always cold
!!$           press(i,j,k) = EOS_Pressure(Spritz_polytrope_handle, rho(i,j,k), eps(i,j,k))
!!$           eps(i,j,k) = EOS_SpecificIntEnergy(Spritz_polytrope_handle,rho(i,j,k),press(i,j,k))
!!$           scon(i,j,k,1) = 0.0d0
!!$           scon(i,j,k,2) = 0.0d0
!!$           scon(i,j,k,3) = 0.0d0
!!$           vel(i,j,k,1) = 0.0d0
!!$           vel(i,j,k,2) = 0.0d0
!!$           vel(i,j,k,3) = 0.0d0
!!$           w_lorentz(i,j,k) = 1.0d0
!!$           if (zero_B_in_atmo .eq. 1) then
!!$              dBx(i,j,k) = 0.0d0
!!$              dBy(i,j,k) = 0.0d0
!!$              dBz(i,j,k) = 0.0d0
!!$              Bx(i,j,k) = 0.0d0
!!$              By(i,j,k) = 0.0d0
!!$              Bz(i,j,k) = 0.0d0
!!$              B2 = 0.0d0
!!$           end if
!!$           !idem here
!!$           tau(i,j,k)  = sqrt(det) * ((rho(i,j,k)+rho(i,j,k)*eps(i,j,k) + &
!!$                press(i,j,k)+B2)*w_lorentz(i,j,k)*w_lorentz(i,j,k) - (press(i,j,k)+0.5d0*B2)) - dens(i,j,k)
!!$           call CCTK_WARN(1,"====BEGIN MESSAGE====")
!!$           write(warnline,'(a20,3g20.7)') ' (ATMO) x,y,z: ',x(i,j,k),y(i,j,k),z(i,j,k)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) y: ',y(i,j,k)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) z: ',z(i,j,k)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) r: ',r(i,j,k)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) det: ',det
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) rho: ',rho(i,j,k)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) dens: ',dens(i,j,k)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) press: ',press(i,j,k)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) eps: ',eps(i,j,k)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) sconx: ',scon(i,j,k,1)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) scony: ',scon(i,j,k,2)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) sconz: ',scon(i,j,k,3)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) velx: ',vel(i,j,k,1)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) vely: ',vel(i,j,k,2)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) velz: ',vel(i,j,k,3)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) w_lorentz: ',w_lorentz(i,j,k)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) tau: ',tau(i,j,k)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) Bx: ',Bx(i,j,k)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) By: ',By(i,j,k)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) Bz: ',Bz(i,j,k)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) dBx: ',dBx(i,j,k)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) dBy: ',dBy(i,j,k)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) dBz: ',dBz(i,j,k)
!!$           call CCTK_WARN(1,warnline)
!!$           write(warnline,'(a20,g20.7)') ' (ATMO) B2: ',B2
!!$           call CCTK_WARN(1,warnline)
!!$           call CCTK_WARN(1,"====END MESSAGE====")
!!$           cycle
              if (not_set_to_atmo_if_c2pIF_fails .ne. 0) then
                 mhd_c2p_mask(i,j,k) = 3
                 cycle
              end if
              !$OMP CRITICAL
              call CCTK_WARN(2,"Not able to converge using the Ideal Fluid Con2Prim, I will try to use the Polytropic version now")
              !$OMP END CRITICAL
              mhd_c2p_mask(i,j,k) = 2
              call Con2Prim_ptPolytype_mhd(Spritz_polytrope_handle, &
                   dens(i,j,k), scon(i,j,k,1), scon(i,j,k,2), scon(i,j,k,3), &
                   dBx(i,j,k), dBy(i,j,k), dBz(i,j,k), &
                   tau(i,j,k), rho(i,j,k), vel(i,j,k,1), vel(i,j,k,2), &
                   vel(i,j,k,3), eps(i,j,k), press(i,j,k), w_lorentz(i,j,k), &
                   uxx, uxy, uxz, uyy, uyz, uzz, det, x(i,j,k), y(i,j,k), &
                   z(i,j,k), r(i,j,k),Spritz_rho_min, &
                   Bx(i,j,k),By(i,j,k),Bz(i,j,k),B2,Bxlow,Bylow,Bzlow, Spritz_reflevel, c2perror, mhd_c2p_mask(i,j,k))
           end if

           if(eps(i,j,k) .lt. 0.0d0 .and..not.CCTK_EQUALS(Spritz_eos_table,"nuc_eos")) then !due to binding energy, energy in nuc_eos can be negative
              if (Spritz_reflevel.ge.Spritz_c2p_warn_from_reflevel) then
                 !$OMP CRITICAL
                 call CCTK_WARN(1,'Con2Prim_MHD: stopping the code.')
                 call CCTK_WARN(1, '   specific internal energy just went below 0! ')
                 write(warnline,'(a28,i2)') 'on carpet reflevel: ',Spritz_reflevel
                 call CCTK_WARN(1,warnline)
                 write(warnline,'(a20,3g16.7)') 'xyz location: ',&
                      x(i,j,k),y(i,j,k),z(i,j,k)
                 call CCTK_WARN(1,warnline)
                 write(warnline,'(a20,g16.7)') 'radius: ',r(i,j,k)
                 call CCTK_WARN(1,warnline)
                 write(warnline,'(a20,3g16.7)') 'velocities: ',&
                      vel(i,j,k,1),vel(i,j,k,2),vel(i,j,k,3)
                 call CCTK_WARN(1,warnline)
                 write(warnline,'(a20,3g16.7)') 'magnetic field: ',&
                      Bx(i,j,k),By(i,j,k),Bz(i,j,k)
                 call CCTK_WARN(1,warnline)
                 !!!!!!
                 write(warnline,'(a20,3i5)') 'i,j,k: ',i,j,k
                 call CCTK_WARN(1,warnline)
                 write(warnline,'(a20,g16.7)') 'rho: ',rho(i,j,k)
                 call CCTK_WARN(1,warnline)
                 write(warnline,'(a20,g16.7)') 'pressure: ',press(i,j,k)
                 call CCTK_WARN(1,warnline)
                 write(warnline,'(a20,g16.7)') 'eps: ',eps(i,j,k)
                 call CCTK_WARN(1,warnline)
                 !!!!!!
                 call CCTK_WARN(Spritz_c2p_warnlevel, " - CASE A - Specific internal energy negative")
                 !$OMP END CRITICAL
              else
                 !$OMP CRITICAL
                 write(warnline,'(a60,i2)') 'Con2Prim: eps negative, but I was told to ignore level: ',Spritz_reflevel
                 call CCTK_WARN(2,warnline)
                 !$OMP END CRITICAL
                 !set this mask to 1 to identify the failure
                 !The above message is set to warning level 2 and so it would not be shown usually
                 mhd_c2p_mask(i,j,k)=1
              endif
           endif


           if((rho(i,j,k).ne.rho(i,j,k)).or. &
                (eps(i,j,k).ne.eps(i,j,k)).or. &
                (vel(i,j,k,1).ne.vel(i,j,k,1)).or. &
                (vel(i,j,k,2).ne.vel(i,j,k,2)).or. &
                (vel(i,j,k,3).ne.vel(i,j,k,3)).or. &
                (Bx(i,j,k).ne.Bx(i,j,k)).or. &
                (By(i,j,k).ne.By(i,j,k)).or. &
                (Bz(i,j,k).ne.Bz(i,j,k)).or. &
                c2perror) then
              !$OMP CRITICAL
              call CCTK_WARN(2,"Failure or NaN in C2P_MHD Ideal Fluid")
              write(warnline,'(a20,g20.7)') ' rho is NaN: ',rho(i,j,k)
              call CCTK_WARN(2,warnline)
              write(warnline,'(a20,g20.7)') ' eps is NaN: ',eps(i,j,k)
              call CCTK_WARN(2,warnline)
              write(warnline,'(a20,g20.7)') 'velx is NaN: ',vel(i,j,k,1)
              call CCTK_WARN(2,warnline)
              write(warnline,'(a20,g20.7)') 'vely is NaN: ',vel(i,j,k,2)
              call CCTK_WARN(2,warnline)
              write(warnline,'(a20,g20.7)') 'velz is NaN: ',vel(i,j,k,3)
              call CCTK_WARN(2,warnline)
              write(warnline,'(a20,g20.7)') '  Bx is NaN: ',Bx(i,j,k)
              call CCTK_WARN(2,warnline)
              write(warnline,'(a20,g20.7)') '  By is NaN: ',By(i,j,k)
              call CCTK_WARN(2,warnline)
              write(warnline,'(a20,g20.7)') '  Bz is NaN: ',Bz(i,j,k)
              call CCTK_WARN(2,warnline)
              write(warnline,'(a20,g20.7)') 'x=',x(i,j,k)
              call CCTK_WARN(2,warnline)
              write(warnline,'(a20,g20.7)') 'y=',y(i,j,k)
              call CCTK_WARN(2,warnline)
              write(warnline,'(a20,g20.7)') 'z=',z(i,j,k)
              call CCTK_WARN(2,warnline)
              call CCTK_WARN(2, 'setting the point to atmosphere')
              !$OMP END CRITICAL

              !set this mask to 1 to identify the failure
              !The above message is set to warning level 2 and so it would not be shown usually
              mhd_c2p_mask(i,j,k)=1

              !atmosphere
              rho(i,j,k) = Spritz_rho_min
              dens(i,j,k) = sqrt(det)*rho(i,j,k)
              !Since I don't know what to do in this case for eps and press
              !I set them using a polytropic EOS which means that the atmosphere is always cold
              !press(i,j,k) = EOS_Pressure(Spritz_polytrope_handle, rho(i,j,k), eps(i,j,k))
              !eps(i,j,k) = EOS_SpecificIntEnergy(Spritz_polytrope_handle,rho(i,j,k),press(i,j,k))

              if(GetTemp_FromThermo.ne.0) then
              ! set hot atmosphere values
                temperature(i,j,k) = Spritz_hot_atmo_temp
                ! y_e(i,j,k) = Spritz_hot_atmo_Y_e
                y_e_con(i,j,k) = y_e(i,j,k) * dens(i,j,k)
                keytemp = 1
                call EOS_Omni_press(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                   rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
                   press(i,j,k),keyerr,anyerr)
              else
                keytemp = 0
                call EOS_Omni_press(Spritz_polytrope_handle,keytemp,Eos_rf_prec,n,&
                   rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),keyerr,anyerr)
                call EOS_Omni_EpsFromPress(Spritz_polytrope_handle,keytemp,Eos_rf_prec,n,&
                   rho(i,j,k),eps(i,j,k),xtemp,xye,press(i,j,k),eps(i,j,k),keyerr,anyerr)
              endif

              scon(i,j,k,1) = 0.0d0
              scon(i,j,k,2) = 0.0d0
              scon(i,j,k,3) = 0.0d0
              vel(i,j,k,1) = 0.0d0
              vel(i,j,k,2) = 0.0d0
              vel(i,j,k,3) = 0.0d0
              w_lorentz(i,j,k) = 1.0d0
              !idem here
              tau(i,j,k)  = sqrt(det) * ((rho(i,j,k)+rho(i,j,k)*eps(i,j,k) + &
                   press(i,j,k)+B2)*w_lorentz(i,j,k)*w_lorentz(i,j,k) - (press(i,j,k)+0.5d0*B2)) - dens(i,j,k)
           end if

           if (Constrain_Temp_now .ne. 0) then
             !Be sure that tau is computed by definition in the case of scheme switching
             !if ((i == i_lowerbound).and.&
             !    (j == j_lowerbound).and.&
             !    (k == k_lowerbound)) then
             !   !Print INFO only one time per iteration
             !   call CCTK_INFO("Recomputing tau, due to T constraining...")
             !end if

             w = w_lorentz(i,j,k)

             !v_j
             vlowx = gxx(i,j,k)*vel(i,j,k,1) + gxy(i,j,k)*vel(i,j,k,2) + gxz(i,j,k)*vel(i,j,k,3)
             vlowy = gxy(i,j,k)*vel(i,j,k,1) + gyy(i,j,k)*vel(i,j,k,2) + gyz(i,j,k)*vel(i,j,k,3)
             vlowz = gxz(i,j,k)*vel(i,j,k,1) + gyz(i,j,k)*vel(i,j,k,2) + gzz(i,j,k)*vel(i,j,k,3)

             !b^t
             bst = w*(Bx(i,j,k)*vlowx + By(i,j,k)*vlowy + Bz(i,j,k)*vlowz) !this is b^0 * alp

             !b^mu b_mu
             bs2 = (Bx(i,j,k)*Bxlow + By(i,j,k)*Bylow + Bz(i,j,k)*Bzlow + bst**2)/w**2  
             tau(i,j,k)  = sqrt(det) * ((rho(i,j,k)*(1.0d0+eps(i,j,k))+press(i,j,k)+bs2)*w*w -&
                           (press(i,j,k)+0.5d0*bs2) - bst**2) - dens(i,j,k)
           
           end if


        end do
     end do
  end do
  !$OMP END PARALLEL DO
end subroutine Conservative2Primitive_mhd

/*@@
@routine    Con2Prim_pt_mhd
@date
@author     Bruno Giacomazzo
@desc
Given all the appropriate data, recover the primitive variables at
a single point.
@enddesc
 @calls
@calledby
@history

@endhistory

@@*/

subroutine Con2Prim_pt_mhd(handle, dens, sx, sy, sz, dBx, dBy, dBz, &
     tau, rho, velx, vely, &
     velz, eps, press, w_lorentz, uxx, uxy, uxz, uyy, &
     uyz, uzz, det, x, y, z, r, epsnegative,Spritz_rho_min, &
     Bx, By, Bz, B2, Bxlow, Bylow, Bzlow, Spritz_polytrope_handle, pmin, &
     Spritz_reflevel, c2perror, outside_ah, mhd_c2p_mask_local)

!!$  use nrmodule
  use Interfaces, only: jacobian_mhd
  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL:: dens, sx, sy, sz, dBx, dBy, dBz, tau, rho, velx, vely, velz, eps, &
       press, uxx, uxy, uxz, uyy, uyz, uzz, det, w_lorentz, x, &
       y, z, r, Spritz_rho_min, Bx, By, Bz, B2, Bxlow, Bylow, Bzlow, &
       udens, usx, usy, usz, utau, s2, BS, pold, rhoold, w_lorentzold, w, ww2, pmin

  CCTK_INT:: i, handle, Spritz_polytrope_handle, Spritz_reflevel, error_level
  CCTK_REAL, optional:: mhd_c2p_mask_local

  CCTK_REAL :: temp,ye !this won't actually be used, it is just for uding SetToAtmo
  CCTK_REAL:: vxlow, vylow, vzlow, relax,fjac_det
  CCTK_REAL, dimension(2):: unk, fvec, q  !, tmp
  CCTK_REAL, dimension(2,2):: fjac,fjac_inv

  CCTK_INT:: count, ierr

  CCTK_REAL independent(2), dependent(3)

  character(len=300) warnline
  logical epsnegative, c2perror, outside_ah, try_safe_guess

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
  !!!!!!

  c2perror=.false.

  try_safe_guess=.true. !if Newton-Raphson fails to converge, try to start again using the "safe" guess values

  count = 0 !how many time relax is changed
  relax=1.0D0 !1.0d0 !relaxation factor


!!$  Undensitize the variables

  udens = dens /sqrt(det)
  usx = sx /sqrt(det)
  usy = sy /sqrt(det)
  usz = sz /sqrt(det)
  utau = tau /sqrt(det)

  !S^j S_j
  s2 = usx*usx*uxx + usy*usy*uyy + usz*usz*uzz + 2.*usx*usy*uxy + &
       2.*usx*usz*uxz + 2.*usy*usz*uyz
  !B^j S_j
  BS = Bx*usx + By*usy + Bz*usz


!!$  Set initial guess
  !pold = max(pmin/100,EOS_Pressure(handle, rho, eps))

  !!!!!!
  xrho = rho
  xeps = eps
  call EOS_Omni_press(handle,keytemp,Eos_rf_prec,n,&
                xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
  pold = max(pmin/100,xpress(1))
  !!!!!!

  rhoold = rho
  w_lorentzold = max(w_lorentz,1.0d0)
  w_lorentz=w_lorentzold

!!$  Compute the Jacobian
  unk(1) = pold
  unk(2) = w_lorentzold

123 call jacobian_mhd(handle,udens,s2,utau,B2,BS,unk,fvec,fjac)


  do i=1,Spritz_countmax

!!$     write(*,*) "unk=",unk
     relax=1.0d0

!!$ The subroutine jacobian_mhd has supplied function values at unk
!!$ in fvec and Jacobian matrix in fjac.

!!$  Newton-Raphson routine

     fjac_det = fjac(1,1)*fjac(2,2)-fjac(1,2)*fjac(2,1)

     fjac_inv(1,1) = fjac(2,2)/fjac_det
     fjac_inv(1,2) = -fjac(1,2)/fjac_det
     fjac_inv(2,1) = -fjac(2,1)/fjac_det
     fjac_inv(2,2) = fjac(1,1)/fjac_det

     q(1) = -fjac_inv(1,1)*fvec(1) - fjac_inv(1,2)*fvec(2)
     q(2) = -fjac_inv(2,1)*fvec(1) - fjac_inv(2,2)*fvec(2)

!!$     write(*,*)
!!$     write(*,*) 'x,y,z=', x,y,z
!!$     write(*,*) 'unk=',unk
!!$     write(*,*) 'i=',i,'fvec=', fvec
!!$     write(*,*) 'q=',q
!!$     write(*,*)

!!$20   tmp=unk+relax*q

!!$     !we use a relaxation factor in order to avoid unphysical values
!!$     ! for the guess on the pressure and on the Lorentz factor,
!!$     ! but at the same time we don't want it to be smaller than 1.0e-3
!!$     ! to avoid too many iterations. If relax<1.0e-2 then we will simply use
!!$     ! the new guess and hope that everything will work well
!!$
!!$     if(((tmp(1)<0.0d0).or.(tmp(2)<0.8d0)).and.(relax>1.0e-2)) then
!!$        relax=0.5d0*relax
!!$        count=count+1
!!$        if (relax<=1.0e-2) then !1.0e-13 is too low, if q=1e-2 then q*relax=1.0e-15
!!$        write(*,*) "count=",count
!!$        write(*,*) "relax=",relax
!!$           write(*,*) 'tmp=',tmp
!!$           write(*,*) 'unk=',unk
!!$           write(*,*) 'fvec=',fvec
!!$           write(*,*) 'fjac(1,:)=',fjac(1,:)
!!$           write(*,*) 'fjac(2,:)=',fjac(2,:)
!!$           write(*,*) 'q=',q
!!$           write(*,*) 'x=',x,'y=',y,'z=',z
!!$           write(*,*) 'relax<=1.0e-2'
!!$           exit
!!$        end if
!!$        goto 20
!!$     end if

     unk=unk+relax*q !Update solution.

     !unk(1)=max(unk(1),0.0d0)
     !unk(1)=max(unk(1),pmin/100)!we avoid pressure smaller than the pressure in the atmosphere pmin (divided by 100 just to be a bit more flexible)
     unk(2)=max(unk(2),1.0d0)!avoid w_lorentz<1.0

!!$     if ((x==1.6d0).and.(y==0.0d0).and.(z==0.0d0)) then
!!$        call CCTK_WARN(1,"---START MESSAGE---")
!!$        write(warnline,'(a20,2g16.7)') 'unk (p and w)', unk(1), unk(2)
!!$        call CCTK_WARN(1,warnline)
!!$        write(warnline,'(a20,2g16.7)') 'q (p and w)', q(1), q(2)
!!$        call CCTK_WARN(1,warnline)
!!$        write(warnline,'(a20,2g16.7)') 'fvec', fvec(1), fvec(2)
!!$        call CCTK_WARN(1,warnline)
!!$        write(warnline,'(a20,2g16.7)') 'fjac1', fjac_inv(1,1), fjac_inv(1,2)
!!$        call CCTK_WARN(1,warnline)
!!$        write(warnline,'(a20,2g16.7)') 'fjac2', fjac_inv(2,1), fjac_inv(2,2)
!!$        call CCTK_WARN(1,warnline)
!!$        call CCTK_WARN(1,"---END MESSAGE---")
!!$     end if

     !Check if the correction q is either a NaN or an Inf
     !It also checks if unk is a NaN or has an unphysical value
     !in this case the NR is considered as not converged and we exit from the do cycle
     if ((q(1).ne.q(1)).or.(q(2).ne.q(2)).or.((q(1)*0.0d0).ne.0.0d0).or.((q(2)*0.0d0).ne.0.0d0).or.&
          (unk(1).ne.unk(1)).or.(unk(2).ne.unk(2)).or.(unk(1)<0.0d0).or.(unk(2)<1.0d0)) then
        !$OMP CRITICAL
        write(warnline,'(a20,2g16.7)') 'q (p and w)', q(1), q(2)
        call CCTK_WARN(2,warnline)
        write(warnline,'(a20,2g16.7)') 'unk (p and w)', unk(1), unk(2)
        call CCTK_WARN(2,warnline)
        call CCTK_WARN(2,"NaNs, Infs or unphysical values in the Newton-Raphson routine")
        !$OMP END CRITICAL
        exit !Not able to converge because of NaNs or Infs produced. Exit from the do loop.
     else
        if (max(abs(q(1)),abs(q(2))) <= Spritz_mhd_tolf) goto 1000 !Check root convergence.
     end if

!!$  Compute fvec and its Jacobian fjac with the new guess
     call jacobian_mhd(handle,udens,s2,utau,B2,BS,unk,fvec,fjac)

  end do

  if(try_safe_guess) then
     !using "safe" guess values suggested in Cerda-Duran et al 2008, A&A 492, 937-953
     B2 = (Bx*Bxlow + By*Bylow + Bz*Bzlow)!B^i*B_i
     !unk(1) = EOS_Pressure(handle, udens, (utau-0.5d0*B2)/udens)!upper limit for the pressure

     !!!!!!
     xrho(1) = udens
     xeps(1) = (utau-0.5d0*B2)/udens
     call EOS_Omni_press(handle,keytemp,Eos_rf_prec,n,&
                   xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
     unk(1) = xpress(1)!upper limit for the pressure
     !!!!!!

     unk(2) = 1.0d5 !upper limit for w_lorentz
     try_safe_guess=.false.
     goto 123
  end if

  !$OMP CRITICAL
  call CCTK_WARN(2,"Not able to converge. Let's see if we are into the atmo")
  !$OMP END CRITICAL

  !I was not able to converge so I check if I'm in the atmosphere
  !simply looking at the conserved density
  !(for efficiency reason this part of the code should be executed
  ! at the beginning, but I'm not sure of what I'm doing here to the energy)
  if(dens<=(sqrt(det)*(Spritz_rho_min*(1.d0+Spritz_atmo_tolerance)))) then
    call SetToAtmo(rho,dens,det,w_lorentz,press,&
              sx,sy,sz,s2,usx,usy,usz,tau,&
              velx,vely,velz,eps,temp,ye,handle,&
              BS,B2,Bx,By,Bz,Bxlow,Bylow,Bzlow)
     return
  end if

  !I'm sorry but I'm not in the atmo
  !Con2Prim failed its mission

  if (Spritz_reflevel.ge.Spritz_c2p_warn_from_reflevel) then
     !$OMP CRITICAL
     call CCTK_WARN(2,"Not able to converge. I also tried to see if we are into the atmo, but it does not seem so.")
     call CCTK_WARN(2, 'Con2Prim_pt_mhd: error: did not converge in ')
     write(warnline,'(a20,g16.7,a10)') '              ',&
          Spritz_countmax,' steps'
     call CCTK_WARN(2,warnline)
     write(warnline,'(a28,i2)') 'on carpet reflevel: ',Spritz_reflevel
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'relaxation factor: ',relax
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'counter: ',count
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g16.7)') 'de is ',dens
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g16.7)') 'sx is ',sx
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g16.7)') 'sy is ',sy
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g16.7)') 'sz is ',sz
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g16.7)') 'en is ',tau
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g16.7)') 'S^2 is ',s2
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g16.7)') 'B^2 is ',B2
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g16.7)') 'B^j*S_j is ',BS
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,2g16.7)') 'func values: ',fvec
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,2g16.7)') 'Jacobian(1,:): ',fjac(1,:)
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,2g16.7)') 'Jacobian(2,:): ',fjac(2,:)
     call CCTK_WARN(2,warnline)
     write(warnline,'(a30,2g16.7)') 'q (i.e. the correction)',q
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,2g16.7)') 'pnew, Wnew: ',unk
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,2g16.7)') 'pold, Wold: ',pold, w_lorentzold
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g16.7)') 'epsilon: ',eps
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g16.7)') 'Bx: ',Bx
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g16.7)') 'By: ',By
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g16.7)') 'Bz: ',Bz
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,3g16.7)') 'xyz location: ',x,y,z
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g16.7)') 'radius: ',r
     call CCTK_WARN(2,warnline)
     call CCTK_WARN(2, "Did not converge")
     !$OMP END CRITICAL
     !call CCTK_WARN(Spritz_c2p_warnlevel, "Did not converge")
     c2perror=.true.
     if(present(mhd_c2p_mask_local)) then
        !set this mask to 1 to identify the failure
        mhd_c2p_mask_local=1
     end if
     return
     !$OMP CRITICAL
     call CCTK_WARN(0,"Why am I here? Return statement did not work.") !!ONLY FOR DEBUG
     !$OMP END CRITICAL
  else
     !$OMP CRITICAL
     write(warnline,'(a70,i2)') 'Con2Prim: not able to converge, but I was told to ignore level: ',Spritz_reflevel
     call CCTK_WARN(2,warnline)!!!!!!!!!!
     call CCTK_WARN(2,'For security reason I reset the point to be atmosphere to avoid NaNs')
     !$OMP END CRITICAL
     if(present(mhd_c2p_mask_local)) then
        !set this mask to 1 to identify the failure
        mhd_c2p_mask_local=1
     end if
     !atmosphere
     call SetToAtmo(rho,dens,det,w_lorentz,press,&
               sx,sy,sz,s2,usx,usy,usz,tau,&
               velx,vely,velz,eps,temp,ye,handle,&
               BS,B2,Bx,By,Bz,Bxlow,Bylow,Bzlow)
     return
  endif

1000 rho = udens/unk(2)
  press = unk(1)
  w_lorentz = unk(2)

  !The Newton-Raphson should never return a negative pressure, but it's better to check
  if (press < 0.0d0) then
     !$OMP CRITICAL
     write(warnline,'(a20,g16.7)') 'rho       : ',rho
     call CCTK_WARN(1,warnline)
     write(warnline,'(a20,g16.7)') 'press     : ',press
     call CCTK_WARN(1,warnline)
     write(warnline,'(a20,g16.7)') 'w_lorentz : ',w_lorentz
     call CCTK_WARN(1,warnline)
     write(warnline,'(a20,3g16.7)') 'x,y,z    : ',x,y,z
     call CCTK_WARN(1,warnline)
     write(warnline,'(a20,g16.7)')  'r        : ',r
     call CCTK_WARN(1,warnline)
     call CCTK_WARN(0,"negative pressure")
     !$OMP END CRITICAL
     !I set the point to use a polytropic eos
     !(it works only if the previous warning level is set to a value larger than 0)
     !rho = Spritz_rho_min
     !press = EOS_Pressure(Spritz_polytrope_handle, rho, eps)
     !eps = EOS_SpecificIntEnergy(Spritz_polytrope_handle,rho,press)

     !!!!!!
     xrho = rho
     xeps = eps
     call EOS_Omni_press(Spritz_polytrope_handle,keytemp,Eos_rf_prec,n,&
                    xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
     press = xpress(1)
     call EOS_Omni_EpsFromPress(Spritz_polytrope_handle,keytemp,Eos_rf_prec,n,&
                           xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     eps = xeps(1)
     !!!!!!

  end if

!!$  if (count>0) then
!!$     write(*,*) 'number of times relax was changed',count
!!$  end if

!!$  Compute velocities (these are the v_j)
  !relativistic enthalpy


  !eps = EOS_SpecificIntEnergy(handle,rho,press)

  !!!!!!
  xrho = rho
  xpress = press
  call EOS_Omni_EpsFromPress(handle,keytemp,Eos_rf_prec,n,&
                        xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
  eps = xeps(1)
  !!!!!!

  w = rho + rho*eps + press

  ww2 = w*w_lorentz**2 !this is Z = rho*h*W_lorentz**2 (Anton et al. 2005)

  vxlow = Bxlow*BS/(ww2*(ww2 + B2)) + usx/(ww2 + B2)
  vylow = Bylow*BS/(ww2*(ww2 + B2)) + usy/(ww2 + B2)
  vzlow = Bzlow*BS/(ww2*(ww2 + B2)) + usz/(ww2 + B2)

  !and now compute the v^j
  velx = vxlow*uxx + vylow*uxy + vzlow*uxz
  vely = vxlow*uxy + vylow*uyy + vzlow*uyz
  velz = vxlow*uxz + vylow*uyz + vzlow*uzz

  if (rho .le. Spritz_rho_min*(1.d0+Spritz_atmo_tolerance)) then
     !atmosphere
     call SetToAtmo(rho,dens,det,w_lorentz,press,&
               sx,sy,sz,s2,usx,usy,usz,tau,&
               velx,vely,velz,eps,temp,ye,handle,&
               BS,B2,Bx,By,Bz,Bxlow,Bylow,Bzlow)
     end if


!!$If all else fails, use the polytropic EoS

  if(eps .lt. 0.0d0) then
     epsnegative = .true.
     if(present(mhd_c2p_mask_local)) then
        !set this mask to 1 to identify the failure
        mhd_c2p_mask_local=1
     end if
     return
  endif

  !FIXME TODO check AH, etc.
  if((vxlow*velx + vylow*vely + vzlow*velz)>=1.0d0) then
     if (Spritz_reflevel.ge.Spritz_c2p_warn_from_reflevel) then
        error_level=2

        !If the conserved density is far away from the atmo
        ! then check if this point is at least inside an AH
        ! otherwise stop the code
!!$        if (udens.ge.1000.0d0*Spritz_rho_min) then
!!$           if (outside_ah == .true.) error_level=1
!!$        end if

        !$OMP CRITICAL
        write(warnline,'(a20,3g16.7)') 'in x,y,z=',x,y,z
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a28,i2)') 'on carpet reflevel: ',Spritz_reflevel
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'w_lorentz =',w_lorentz
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'v2 =',vxlow*velx + vylow*vely + vzlow*velz
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'vxlow =',vxlow
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'vylow =',vylow
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'vzlow =',vzlow
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'velx =',velx
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'vely =',vely
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'velz =',velz
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'Bxlow =',Bxlow
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'Bylow =',Bylow
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'Bzlow =',Bzlow
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'Bx =',Bx
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'By =',By
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'Bz =',Bz
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'det(g_{ij}) =',det
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'uxx =',uxx
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'uxy =',uxy
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'uxz =',uxz
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'uyy =',uyy
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'uyz =',uyz
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'uzz =',uzz
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'udens =',udens
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'Sx =',usx
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'Sy =',usy
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'Sz =',usz
        call CCTK_WARN(error_level,warnline)
        if (error_level==1) then
           call CCTK_WARN(0,"Velocity greater than the speed of light!")
        end if
        !The following warning is called only if error_level=2
        call CCTK_WARN(error_level,"But this point is inside an AH or so close to atmo that I ignored it and set the point to atmo")
        ! check that it is really inside the atmosphere!!!
        !$OMP END CRITICAL
        if(present(mhd_c2p_mask_local)) then
           !set this mask to 1 to identify the failure
           mhd_c2p_mask_local=1
        end if

        !atmosphere
        call SetToAtmo(rho,dens,det,w_lorentz,press,&
                  sx,sy,sz,s2,usx,usy,usz,tau,&
                  velx,vely,velz,eps,temp,ye,handle,&
                  BS,B2,Bx,By,Bz,Bxlow,Bylow,Bzlow)
        return
     else
        !$OMP CRITICAL
        write(warnline,'(a60,i2)') 'Con2Prim: v^2>1, but I was told to ignore level: ',Spritz_reflevel
        call CCTK_WARN(2,warnline)
        !$OMP END CRITICAL
        if(present(mhd_c2p_mask_local)) then
           !set this mask to 1 to identify the failure
           mhd_c2p_mask_local=1
        end if
        !To Avoid NaN I set this point to be atmosphere
        ! If this region is covered by a finer grid, the point would be correcly set to its right value
        ! by the Con2Prim routine execute on that level
        !atmosphere
        call SetToAtmo(rho,dens,det,w_lorentz,press,&
                  sx,sy,sz,s2,usx,usy,usz,tau,&
                  velx,vely,velz,eps,temp,ye,handle,&
                  BS,B2,Bx,By,Bz,Bxlow,Bylow,Bzlow)
        return
     end if
  end if

  !check whether rho is NaN
  if ((rho .ne. rho).or. &
       (eps .ne. eps).or. &
       (velx .ne. velx).or. &
       (vely .ne. vely).or. &
       (velz .ne. velz).or. &
       (Bx .ne. Bx).or. &
       (By .ne. By).or. &
       (Bz .ne. Bz)) then
     !$OMP CRITICAL
     write(warnline,'(a20,g20.7)') ' rho is NaN: ',rho
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') ' eps is NaN: ',eps
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') 'velx is NaN: ',velx
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') 'vely is NaN: ',vely
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') 'velz is NaN: ',velz
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') '  Bx is NaN: ',Bx
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') '  By is NaN: ',By
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') '  Bz is NaN: ',Bz
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') 'x=',x
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') 'y=',y
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') 'z=',z
     call CCTK_WARN(2,warnline)
     call CCTK_WARN(2, 'setting the point to atmosphere')
     !$OMP END CRITICAL
     if(present(mhd_c2p_mask_local)) then
        !set this mask to 1 to identify the failure
        !The above message is set to warning level 2 and so it would not be shown usually
        mhd_c2p_mask_local=1
     end if

     !atmosphere
     call SetToAtmo(rho,dens,det,w_lorentz,press,&
               sx,sy,sz,s2,usx,usy,usz,tau,&
               velx,vely,velz,eps,temp,ye,handle,&
               BS,B2,Bx,By,Bz,Bxlow,Bylow,Bzlow)
     return
  end if

end subroutine Con2Prim_pt_mhd


/*@@
@routine    Con2Prim_pt_hot_Palenzuela
@date
@author     Lorenzo Sala
@desc
C2P described by Siegel et al 2018, uses Brent Method in
one dimension for recovery. (First described by Palenzuela, GRMHD 2015,
used for first time by Galeazzi/Kastaun 2013, GRHD)
@enddesc
 @calls
@calledby
@history

@endhistory

@@*/

subroutine Con2Prim_pt_hot_Palenzuela(handle, dens, sx, sy, sz, dBx, dBy, dBz, &
     tau, rho, velx, vely, velz, &
     eps, press, w_lorentz, ye, temp, entropy, uxx, uxy, uxz, uyy, &
     uyz, uzz, det, x, y, z, r, Spritz_rho_min, &
     Bx, By, Bz, B2, Bxlow, Bylow, Bzlow, &
     Spritz_reflevel)

!!$  use nrmodule
  use Interfaces, only: BrentMethod
  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL:: dens, sx, sy, sz, dBx, dBy, dBz, tau, rho, velx, vely, velz, eps, &
       press, uxx, uxy, uxz, uyy, uyz, uzz, det, w_lorentz, x, &
       y, z, r, Spritz_rho_min, Bx, By, Bz, B2, Bxlow, Bylow, Bzlow, &
       udens, usx, usy, usz, utau, s2, BS, ye, temp, entropy

  CCTK_REAL:: vxlow, vylow, vzlow, relax
  CCTK_REAL :: pale_q, pale_r, pale_s, pale_t, Wminus2, pale_Z, pale_x
  CCTK_REAL :: xlow,xup,flow,fup,dx

  CCTK_INT:: i, handle, Spritz_reflevel, error_level


  character(len=300) warnline
  logical c2perror, outside_ah

  !!!!!!
  ! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1),xent(1),dummy(1)

  n = 1;keytemp = 0;anyerr = 0;keyerr = 0
  xpress = 0.0d0;xeps = 0.0d0;xent = 0.0d0;xtemp = temp; xye = ye
  ! end EOS Omni vars
  !!!!!!

  if(.not.CCTK_EQUALS(Spritz_eos_table,"nuc_eos")) call CCTK_WARN(0,"Palenzuela C2P not supported with non-tabulated EoS.")

!!$  Undensitize the variables
  udens = dens /sqrt(det)
  usx = sx /sqrt(det)
  usy = sy /sqrt(det)
  usz = sz /sqrt(det)
  utau = tau /sqrt(det)

  !S^j S_j
  s2 = usx*usx*uxx + usy*usy*uyy + usz*usz*uzz + 2.*usx*usy*uxy + &
       2.*usx*usz*uxz + 2.*usy*usz*uyz
  !B^j S_j
  BS = Bx*usx + By*usy + Bz*usz

  pale_q = utau / udens
  pale_r = S2 / udens**2
  pale_s = B2/udens
  pale_t = BS / udens**1.5d0

  !first, let's check bounds
  xlow = max(1.0d-10,1.0d0 + pale_q - pale_s) !this is for sure a lower bound
  xup  = 2.0d0 + 2.0d0*pale_q - pale_s !this corresponds to the dominant energy condition (rho*(1+eps)>p)
  call f_Palenzuela_value(xlow,pale_q,pale_r,pale_s,pale_t,udens,ye,temp,flow)
  call f_Palenzuela_value(xup ,pale_q,pale_r,pale_s,pale_t,udens,ye,temp,fup )

  !just to quickly check that dominant energy condition has not been violated
  if(flow*fup.gt.0.0d0) then
    xup = 3.0d0 + 3.0d0*pale_q - 1.5d0*pale_s !this corresponds to condition rho*(1+eps)>p/2
    call f_Palenzuela_value(xup ,pale_q,pale_r,pale_s,pale_t,udens,ye,temp,fup )
    if(flow*fup.le.0.0d0) then
      write(warnline,'(a70,3g16.7,i2)') 'Dominant energy condition has been violated! x,y,z,reflevel',x,y,z,Spritz_reflevel
      call CCTK_WARN(2,warnline)
      write(warnline,'(a20,5g16.7)') 'q,r,s,t,D',pale_q,pale_r,pale_s,pale_t,udens
      call CCTK_WARN(2,warnline)
    end if
  end if

  !atmo if there's no chance of finding root
  if(flow*fup.gt.0.0d0) then
    write(warnline,'(a50,3g16.7,i2)') 'C2P_hot setting to atmosphere. x,y,z,reflevel',x,y,z,Spritz_reflevel
    call CCTK_WARN(2,warnline)
    call SetToAtmo(rho,dens,det,w_lorentz,press,&
              sx,sy,sz,s2,usx,usy,usz,tau,&
              velx,vely,velz,eps,temp,ye,handle,&
              BS,B2,Bx,By,Bz,Bxlow,Bylow,Bzlow)
    goto 1234
  end if

  !Let's compute iteratively x, using Brent's Method
  call BrentMethod(xlow,xup,pale_q,pale_r,pale_s,pale_t,udens,ye,temp,3.0d-15,Spritz_mhd_tolf, pale_x)

  !now compute back primitive from x,q,r,s,t
  Wminus2 = 1.0d0 - ( pale_x*pale_x*pale_r + (2.0d0*pale_x+pale_s)*pale_t*pale_t ) / ( pale_x*pale_x*(pale_x+pale_s)*(pale_x+pale_s) )
  Wminus2 = min(max(Wminus2 ,1.0d-13 ), 1.0d0-1d-10)
  w_lorentz = Wminus2**(-0.5d0)

  rho = max(udens/w_lorentz,rho_abs_eos_min)

  !find eps and press
  eps = w_lorentz - 1.0d0 + (1.0d0-w_lorentz**2)*pale_x/w_lorentz + w_lorentz*(pale_q - pale_s + pale_t**2/(2.0d0*pale_x**2) + pale_s/(2.0d0*w_lorentz**2))

  !from eps and rho, find temperature and pressure
  xrho = rho
  xeps = eps
  xtemp = temp
  xye = ye
  keytemp = 0
  call EOS_Omni_press(handle,keytemp,Eos_rf_prec,n,&
  xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
  press = xpress(1)
  temp   = xtemp(1)
  !done.

  !revert to cold_constrainT C2P if temperature is too low. If any other errors, abort.
  if(anyerr.ne.0) then
    if(keyerr(1).eq.668) then
       !TODO: if this point is covered by finer grids, don't care and go on (saves time).
       !now, constrain temperature to atmo and try with the other c2p scheme.
       if(Spritz_reflevel.ge.Spritz_c2p_warn_from_reflevel) then
         write(warnline,"(1P10E15.6)") "Energy is too low. Resetting temp and using constrC2P. reflevel:",Spritz_reflevel
         call CCTK_WARN(2, warnline)
         temp = Spritz_hot_atmo_temp
         call Con2Prim_pt_cold_ConstrainT(handle, dens, sx, sy, sz, dBx, dBy, dBz, &
                  tau, rho, velx, vely, velz, &
                  eps, press, w_lorentz, ye, temp, entropy, uxx, uxy, uxz, uyy, &
                  uyz, uzz, det, x, y, z, r, Spritz_rho_min, &
                  Bx, By, Bz, B2, Bxlow, Bylow, Bzlow, &
                  Spritz_reflevel)
                  ! TODO should I recompute tau?
         call CCTK_WARN(2, "Should I recompute tau?")
       else
         call SetToAtmo(rho,dens,det,w_lorentz,press,&
                   sx,sy,sz,s2,usx,usy,usz,tau,&
                   velx,vely,velz,eps,temp,ye,handle,&
                   BS,B2,Bx,By,Bz,Bxlow,Bylow,Bzlow)
       end if
       goto 1234
     elseif (keyerr(1).eq.106) then ! if rho < rho_abs_eos_min, should never happen
       call SetToAtmo(rho,dens,det,w_lorentz,press,&
                 sx,sy,sz,s2,usx,usy,usz,tau,&
                 velx,vely,velz,eps,temp,ye,handle,&
                 BS,B2,Bx,By,Bz,Bxlow,Bylow,Bzlow)
       goto 1234
     else
       !$OMP CRITICAL
	   write(warnline,'(a30,3g16.7,i2)') 'Eos error in x y z, reflevel',x,y,z,Spritz_reflevel
       call CCTK_WARN(1,warnline)
	   write(warnline,'(a40,5g16.7)') 'udens, utau, w_lorentz, eps, temp',udens,utau,w_lorentz,eps,temp
       call CCTK_WARN(1,warnline)       
       write(warnline,'(a10,i4)') 'keyerr:', keyerr(1)
       call CCTK_WARN(1,warnline) 
       call CCTK_WARN(0,"EOS_Omni error in Palenzuela C2P.")
       !$OMP END CRITICAL
     end if
  end if

  pale_Z = pale_x * rho * w_lorentz

  !now covariant velocities
  vxlow = Bxlow*BS/(pale_Z*(pale_Z + B2)) + usx/(pale_Z + B2)
  vylow = Bylow*BS/(pale_Z*(pale_Z + B2)) + usy/(pale_Z + B2)
  vzlow = Bzlow*BS/(pale_Z*(pale_Z + B2)) + usz/(pale_Z + B2)

  !and now compute the v^j
  velx = vxlow*uxx + vylow*uxy + vzlow*uxz
  vely = vxlow*uxy + vylow*uyy + vzlow*uyz
  velz = vxlow*uxz + vylow*uyz + vzlow*uzz

  if(rho<=Spritz_rho_min*(1.d0+Spritz_atmo_tolerance)) then
    call SetToAtmo(rho,dens,det,w_lorentz,press,&
              sx,sy,sz,s2,usx,usy,usz,tau,&
              velx,vely,velz,eps,temp,ye,handle,&
              BS,B2,Bx,By,Bz,Bxlow,Bylow,Bzlow)
  end if

  !Pressure has been inferred from EoS, should not have problems. but we can check.
  if (press < 0.0d0) then
     !$OMP CRITICAL
     write(warnline,'(a20,g16.7)') 'rho       : ',rho
     call CCTK_WARN(1,warnline)
     write(warnline,'(a20,g16.7)') 'press     : ',press
     call CCTK_WARN(1,warnline)
     write(warnline,'(a20,g16.7)') 'w_lorentz : ',w_lorentz
     call CCTK_WARN(1,warnline)
     write(warnline,'(a20,3g16.7)') 'x,y,z    : ',x,y,z
     call CCTK_WARN(1,warnline)
     write(warnline,'(a20,g16.7)')  'r        : ',r
     call CCTK_WARN(1,warnline)
     call CCTK_WARN(0,"negative pressure")
     !$OMP END CRITICAL
  end if

  !FIXME FIXME qui si dice che vengono fatte cose che non sono fatte, come controllare se si è in AH
  if((vxlow*velx + vylow*vely + vzlow*velz)>=1.0d0) then
     if (Spritz_reflevel.ge.Spritz_c2p_warn_from_reflevel) then
        error_level=2
        !If the conserved density is far away from the atmo
        ! then check if this point is at least inside an AH
        ! otherwise stop the code
!!$        if (udens.ge.1000.0d0*Spritz_rho_min) then
!!$           if (outside_ah == .true.) error_level=1
!!$        end if

        !$OMP CRITICAL
        write(warnline,'(a20,3g16.7)') 'in x,y,z=',x,y,z
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a28,i2)') 'on carpet reflevel: ',Spritz_reflevel
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'w_lorentz =',w_lorentz
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'v2 =',vxlow*velx + vylow*vely + vzlow*velz
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'vxlow =',vxlow
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'vylow =',vylow
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'vzlow =',vzlow
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'velx =',velx
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'vely =',vely
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'velz =',velz
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'Bxlow =',Bxlow
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'Bylow =',Bylow
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'Bzlow =',Bzlow
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'Bx =',Bx
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'By =',By
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'Bz =',Bz
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'det(g_{ij}) =',det
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'uxx =',uxx
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'uxy =',uxy
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'uxz =',uxz
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'uyy =',uyy
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'uyz =',uyz
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'uzz =',uzz
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'udens =',udens
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'Sx =',usx
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'Sy =',usy
        call CCTK_WARN(error_level,warnline)
        write(warnline,'(a20,1g16.7)') 'Sz =',usz
        call CCTK_WARN(error_level,warnline)
        if (error_level==1) then
           call CCTK_WARN(0,"Velocity greater than the speed of light!")
        end if
        !The following warning is called only if error_level=2
        call CCTK_WARN(error_level,"v>c. But this point is inside an AH or so close to atmo that I ignored it and set the point to atmo.")
        !$OMP END CRITICAL
        call SetToAtmo(rho,dens,det,w_lorentz,press,&
                    sx,sy,sz,s2,usx,usy,usz,tau,&
                    velx,vely,velz,eps,temp,ye,handle,&
                    BS,B2,Bx,By,Bz,Bxlow,Bylow,Bzlow)
        goto 1234
     else
        !$OMP CRITICAL
        write(warnline,'(a60,i2)') 'Con2Prim: v^2>1, but I was told to ignore level: ',Spritz_reflevel
        call CCTK_WARN(2,warnline)
        !$OMP END CRITICAL
        !To Avoid NaN I set this point to be atmosphere
        ! If this region is covered by a finer grid, the point would be correcly set to its right value
        ! by the Con2Prim routine execute on that level
        !atmosphere
        call SetToAtmo(rho,dens,det,w_lorentz,press,&
                    sx,sy,sz,s2,usx,usy,usz,tau,&
                    velx,vely,velz,eps,temp,ye,handle,&
                    BS,B2,Bx,By,Bz,Bxlow,Bylow,Bzlow)
        goto 1234
     end if
  end if

  !check NaNs
  if ((rho .ne. rho).or. &
       (eps .ne. eps).or. &
       (velx .ne. velx).or. &
       (vely .ne. vely).or. &
       (velz .ne. velz).or. &
       (Bx .ne. Bx).or. &
       (By .ne. By).or. &
       (Bz .ne. Bz)) then
     !$OMP CRITICAL
     write(warnline,'(a20,g20.7)') ' rho is NaN: ',rho
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') ' eps is NaN: ',eps
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') 'velx is NaN: ',velx
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') 'vely is NaN: ',vely
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') 'velz is NaN: ',velz
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') '  Bx is NaN: ',Bx
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') '  By is NaN: ',By
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') '  Bz is NaN: ',Bz
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') 'x=',x
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') 'y=',y
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') 'z=',z
     call CCTK_WARN(2,warnline)
     call CCTK_WARN(2, 'setting the point to atmosphere')
     !$OMP END CRITICAL
     ! if(present(mhd_c2p_mask_local)) then
     !    !set this mask to 1 to identify the failure
     !    !The above message is set to warning level 2 and so it would not be shown usually
     !    mhd_c2p_mask_local=1
     ! end if

     call SetToAtmo(rho,dens,det,w_lorentz,press,&
                 sx,sy,sz,s2,usx,usy,usz,tau,&
                 velx,vely,velz,eps,temp,ye,handle,&
                 BS,B2,Bx,By,Bz,Bxlow,Bylow,Bzlow)
     goto 1234
  end if

  !finally, let's find entropy
1234 xrho = rho
  xtemp = temp
  xye = ye
  keytemp = 1
  call EOS_Omni_EOS_short(handle,keytemp,Eos_rf_prec,n,&
       xrho,xeps,xtemp,xye,xpress,xent,dummy,dummy,dummy,dummy,dummy,&
       keyerr,anyerr)
  entropy = xent(1)
end subroutine Con2Prim_pt_hot_Palenzuela

subroutine BrentMethod ( a, b, cons_q, cons_r, cons_s, cons_t, cons_udens, ye, temp, machep, t, x)
!*****************************************************************************
!! ZERO seeks the root of a function F(X) in an interval [A,B].
!
!  Discussion:
!    The interval [A,B] must be a change of sign interval for F.
!    That is, F(A) and F(B) must be of opposite signs.  Then
!    assuming that F is continuous implies the existence of at least
!    one value C between A and B for which F(C) = 0.
!
!    The location of the zero is determined to within an accuracy
!    of 6 * MACHEPS * abs ( C ) + 2 * T.
!
!    Thanks to Thomas Secretin for pointing out a transcription error in the
!    setting of the value of P, 11 February 2013.
!
!  Licensing: This code is distributed under the GNU LGPL license.
!
!  Modified: 11 February 2013
!
!  Author:
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!    Input, real ( kind = 8 ) A, B, the endpoints of the change of sign interval.
!    Input, real ( kind = 8 ) MACHEP, an estimate for the relative machine precision.
!    Input, real ( kind = 8 ) T, a positive error tolerance.
!    function whose zero is being sought.
!    Output, real ( kind = 8 ) ZERO, the estimated value of a zero of
!    the function F.
!
  implicit none

  CCTK_REAL, intent(in)  :: a, b, cons_q, cons_r, cons_s, cons_t, cons_udens, ye, machep, t
  CCTK_REAL, intent(out) :: x
  CCTK_REAL, intent(inout) :: temp

  CCTK_REAL :: c,d,e,f,fa,fb,fc,m,&
               p,q,r,s,sa,sb,tol

!
!  Make local copies of A and B.
!
  sa = a
  sb = b
  call f_Palenzuela_value(sa, cons_q, cons_r, cons_s, cons_t, cons_udens, ye, temp, fa)
  call f_Palenzuela_value(sb, cons_q, cons_r, cons_s, cons_t, cons_udens, ye, temp, fb)

  c = sa
  fc = fa
  e = sb - sa
  d = e

  do
    if ( abs ( fc ) < abs ( fb ) ) then
      sa = sb
      sb = c
      c = sa
      fa = fb
      fb = fc
      fc = fa
    end if

    tol = 2.0D+00 * machep * abs ( sb ) + t
    m = 0.5D+00 * ( c - sb )

    if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
      exit
    end if

    if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then
      e = m
      d = e
    else
      s = fb / fa
      if ( sa == c ) then
        p = 2.0D+00 * m * s
        q = 1.0D+00 - s
      else
        q = fa / fc
        r = fb / fc
        p = s * ( 2.0D+00 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
        q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )
      end if
      if ( 0.0D+00 < p ) then
        q = - q
      else
        p = - p
      end if
      s = e
      e = d
      if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. &
        p < abs ( 0.5D+00 * s * q ) ) then
        d = p / q
      else
        e = m
        d = e
      end if
    end if
    sa = sb
    fa = fb
    if ( tol < abs ( d ) ) then
      sb = sb + d
    else if ( 0.0D+00 < m ) then
      sb = sb + tol
    else
      sb = sb - tol
    end if
    call f_Palenzuela_value(sb, cons_q, cons_r, cons_s, cons_t, cons_udens, ye, temp, fb)
    if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
         ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
      c = sa
      fc = fa
      e = sb - sa
      d = e
    end if
  end do
  x = sb !this returns x
end subroutine BrentMethod

subroutine f_Palenzuela_value ( x, q, r, s, t, udens, ye, temp, fval )
  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, intent(in)    :: x, q, r, s, t, udens, ye
  CCTK_REAL, intent(inout) :: temp
  CCTK_REAL, intent(out)   :: fval
  CCTK_REAL :: Wminus2, W, rho, eps, press, h

  !!!!!!
  ! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1),handle
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1)

  n = 1;keytemp = 0;anyerr = 0;keyerr = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = temp; xye = ye
  ! end EOS Omni vars
  !!!!!!

  Wminus2 = 1.0d0 - ( x*x*r + (2.0d0*x+s)*t*t ) / ( x*x*(x+s)*(x+s) )
  Wminus2 = min(max(Wminus2 ,1.0d-13 ), 1.0d0-1d-10)
  W       = Wminus2**(-0.5d0)
  rho = udens/W
  eps = W - 1.0d0 + (1.0d0-W**2)*x/W + W*(q - s + t**2/(2.0d0*x**2) + s/(2.0d0*W**2))

  !find temp from initial guess and pressure
  xrho = max(rho,rho_abs_eos_min) !so that EOS_Omni does not fail
  xeps = eps
  xtemp = temp
  xye = ye
  handle  = 4
  keytemp = 0
  call EOS_Omni_press(handle,keytemp,Eos_rf_prec,n,&
  xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
  press = xpress(1)
  !ok

  h = 1.0d0 + eps + press/rho
  fval = x - h*W
end subroutine f_Palenzuela_value

/*@@
@routine    Con2Prim_pt_cold_ConstrainT
@date
@author     Lorenzo Sala
@desc
Constrain temperature, solving in one dimension for w_lorentz
@enddesc
 @calls
@calledby
@history

@endhistory

@@*/

subroutine Con2Prim_pt_cold_ConstrainT(handle, dens, sx, sy, sz, dBx, dBy, dBz, &
     tau, rho, velx, vely, velz, &
     eps, press, w_lorentz, ye, temp, entropy, uxx, uxy, uxz, uyy, &
     uyz, uzz, det, x, y, z, r, Spritz_rho_min, &
     Bx, By, Bz, B2, Bxlow, Bylow, Bzlow, &
     Spritz_reflevel)

!!$  use nrmodule
  use Interfaces, only: BrentMethodCold
  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL:: dens, sx, sy, sz, dBx, dBy, dBz, tau, rho, velx, vely, velz, eps, &
       press, uxx, uxy, uxz, uyy, uyz, uzz, det, w_lorentz, x, &
       y, z, r, Spritz_rho_min, Bx, By, Bz, B2, Bxlow, Bylow, Bzlow, &
       udens, usx, usy, usz, s2, BS, ye, temp, entropy

  CCTK_REAL:: vxlow, vylow, vzlow, relax
  CCTK_REAL :: Z_val
  CCTK_REAL :: wlow,wup,flow,fup

  CCTK_INT:: i, handle, Spritz_reflevel, error_level


  character(len=300) warnline
  logical c2perror, outside_ah

  !!!!!!
  ! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1),xent(1),dummy(1)

  n = 1;keytemp = 0;anyerr = 0;keyerr = 0
  xpress = 0.0d0;xeps = 0.0d0;xent = 0.0d0;xtemp = temp; xye = ye
  ! end EOS Omni vars
  !!!!!!

  if(.not.CCTK_EQUALS(Spritz_eos_table,"nuc_eos")) call CCTK_WARN(0,"Constrained temp C2P not supported with non-tabulated EoS.")

!!$  Undensitize the variables
  udens = dens /sqrt(det)
  usx = sx /sqrt(det)
  usy = sy /sqrt(det)
  usz = sz /sqrt(det)

  !S^j S_j
  s2 = usx*usx*uxx + usy*usy*uyy + usz*usz*uzz + 2.*usx*usy*uxy + &
       2.*usx*usz*uxz + 2.*usy*usz*uyz
  !B^j S_j
  BS = Bx*usx + By*usy + Bz*usz

  !first, let's guess bounds
  wlow = 1.0d0
  wup  = 1.5d0
  call f_cold_eos_mhd_value(wlow,udens,S2,B2,BS,ye,temp,flow)
  call f_cold_eos_mhd_value(wup ,udens,S2,B2,BS,ye,temp,fup )

  do i = 1, Spritz_countmax
    if(flow*fup.le.0.0d0) exit

    write(warnline,'(a40,2i2,3g16.7)') 'C2P constT not converging, i,rl,x,y,z',i,Spritz_reflevel,x,y,z
    call CCTK_WARN(2,warnline)
    write(warnline,'(a20,3g16.7)') 'wup,flow,fup:',wup,flow,fup
    call CCTK_WARN(2,warnline)

    wup  = wup * 1.2d0
    call f_cold_eos_mhd_value(wup ,udens,S2,B2,BS,ye,temp,fup )
  end do

  !Let's compute iteratively x, using Brent's Method
  if(flow*fup.gt.0.0d0) then
    write(warnline,'(a40,3g16.7,i2)') 'Setting to atmosphere. x,y,z,reflevel',x,y,z,Spritz_reflevel
    call CCTK_WARN(2,warnline)
    call SetToAtmo(rho,dens,det,w_lorentz,press,&
              sx,sy,sz,s2,usx,usy,usz,tau,&
              velx,vely,velz,eps,temp,ye,handle,&
              BS,B2,Bx,By,Bz,Bxlow,Bylow,Bzlow)
    return
  end if

  if(flow.eq.0.0d0) then !checking can save time.
    w_lorentz = wlow
  else
    call BrentMethodCold(wlow,wup,udens,s2,B2,BS,ye,temp,3.0d-15,Spritz_mhd_tolf,w_lorentz)
  end if

  !now I have w_lorentz
  rho = max(udens/w_lorentz, rho_abs_eos_min)

  xrho = rho
  xtemp = temp
  xye = ye
  handle  = 4
  keytemp = 1
  call EOS_Omni_EOS_short(handle,keytemp,Eos_rf_prec,n,&
       xrho,xeps,xtemp,xye,xpress,xent,dummy,dummy,dummy,dummy,dummy,&
       keyerr,anyerr)
  press   = xpress(1)
  eps     = xeps(1)
  entropy = xent(1)

  if(anyerr.ne.0) then
      !$OMP CRITICAL
      write(*,*) "Eos error in x y z",x,y,z
      write(*,*) "udens,w_lorentz,eps,temp",udens,w_lorentz,eps,temp
      write(*,*) "keyerr:",keyerr(1)
      call CCTK_WARN(0, "EOS_Omni error in cold eos mhd C2P.")
      !$OMP END CRITICAL
  end if

  Z_val = (rho+rho*eps+press)*w_lorentz*w_lorentz
  vxlow = Bxlow*BS/(Z_val*(Z_val + B2)) + usx/(Z_val + B2)
  vylow = Bylow*BS/(Z_val*(Z_val + B2)) + usy/(Z_val + B2)
  vzlow = Bzlow*BS/(Z_val*(Z_val + B2)) + usz/(Z_val + B2)

  !and now compute the v^j
  velx = vxlow*uxx + vylow*uxy + vzlow*uxz
  vely = vxlow*uxy + vylow*uyy + vzlow*uyz
  velz = vxlow*uxz + vylow*uyz + vzlow*uzz

  if(rho<=Spritz_rho_min*(1.d0+Spritz_atmo_tolerance)) then
    call SetToAtmo(rho,dens,det,w_lorentz,press,&
              sx,sy,sz,s2,usx,usy,usz,tau,&
              velx,vely,velz,eps,temp,ye,handle,&
              BS,B2,Bx,By,Bz,Bxlow,Bylow,Bzlow)
    goto 12345
  end if

  !Pressure has been inferred from EoS, should not have problems. but we can check.
  if (press < 0.0d0) then
     !$OMP CRITICAL
     write(warnline,'(a20,g16.7)') 'rho       : ',rho
     call CCTK_WARN(1,warnline)
     write(warnline,'(a20,g16.7)') 'press     : ',press
     call CCTK_WARN(1,warnline)
     write(warnline,'(a20,g16.7)') 'w_lorentz : ',w_lorentz
     call CCTK_WARN(1,warnline)
     write(warnline,'(a20,3g16.7)') 'x,y,z    : ',x,y,z
     call CCTK_WARN(1,warnline)
     write(warnline,'(a20,g16.7)')  'r        : ',r
     call CCTK_WARN(1,warnline)
     call CCTK_WARN(0,"negative pressure")
     !$OMP END CRITICAL
  end if

  if((vxlow*velx + vylow*vely + vzlow*velz)>=1.0d0) then
    !call CCTK_WARN(2,warnline)
    call CCTK_WARN(2, 'Speed is greater than c. Setting to atmo.')
    call SetToAtmo(rho,dens,det,w_lorentz,press,&
                    sx,sy,sz,s2,usx,usy,usz,tau,&
                    velx,vely,velz,eps,temp,ye,handle,&
                    BS,B2,Bx,By,Bz,Bxlow,Bylow,Bzlow)
    goto 12345
  end if

  !check NaNs
  if ((rho .ne. rho).or. &
       (eps .ne. eps).or. &
       (velx .ne. velx).or. &
       (vely .ne. vely).or. &
       (velz .ne. velz).or. &
       (Bx .ne. Bx).or. &
       (By .ne. By).or. &
       (Bz .ne. Bz)) then
     !$OMP CRITICAL
     write(warnline,'(a20,g20.7)') ' rho is NaN: ',rho
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') ' eps is NaN: ',eps
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') 'velx is NaN: ',velx
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') 'vely is NaN: ',vely
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') 'velz is NaN: ',velz
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') '  Bx is NaN: ',Bx
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') '  By is NaN: ',By
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') '  Bz is NaN: ',Bz
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') 'x=',x
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') 'y=',y
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g20.7)') 'z=',z
     call CCTK_WARN(2,warnline)
     call CCTK_WARN(2, 'setting the point to atmosphere')
     !$OMP END CRITICAL
     ! if(present(mhd_c2p_mask_local)) then
     !    !set this mask to 1 to identify the failure
     !    !The above message is set to warning level 2 and so it would not be shown usually
     !    mhd_c2p_mask_local=1
     ! end if

     call SetToAtmo(rho,dens,det,w_lorentz,press,&
                 sx,sy,sz,s2,usx,usy,usz,tau,&
                 velx,vely,velz,eps,temp,ye,handle,&
                 BS,B2,Bx,By,Bz,Bxlow,Bylow,Bzlow)
     goto 12345
  end if

  !finally, let's find entropy
12345 xrho = rho
  xtemp = temp
  xye = ye
  keytemp = 1
  call EOS_Omni_EOS_short(handle,keytemp,Eos_rf_prec,n,&
       xrho,xeps,xtemp,xye,xpress,xent,dummy,dummy,dummy,dummy,dummy,&
       keyerr,anyerr)
  entropy = xent(1)

end subroutine Con2Prim_pt_cold_ConstrainT

subroutine BrentMethodCold ( a, b, udens, s2, B2, BS, ye, temp, machep, t, x)
  implicit none

  CCTK_REAL, intent(in)  :: a, b, udens, s2, B2, BS, ye, temp, machep, t
  CCTK_REAL, intent(out) :: x

  CCTK_REAL :: c,d,e,f,fa,fb,fc,m,&
               p,q,r,s,sa,sb,tol

  sa = a
  sb = b
  call f_cold_eos_mhd_value(sa,udens,S2,B2,BS,ye,temp,fa)
  call f_cold_eos_mhd_value(sb,udens,S2,B2,BS,ye,temp,fb)

  c = sa
  fc = fa
  e = sb - sa
  d = e

  do
    if ( abs ( fc ) < abs ( fb ) ) then
      sa = sb
      sb = c
      c = sa
      fa = fb
      fb = fc
      fc = fa
    end if

    tol = 2.0D+00 * machep * abs ( sb ) + t
    m = 0.5D+00 * ( c - sb )

    if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
      exit
    end if

    if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then
      e = m
      d = e
    else
      s = fb / fa
      if ( sa == c ) then
        p = 2.0D+00 * m * s
        q = 1.0D+00 - s
      else
        q = fa / fc
        r = fb / fc
        p = s * ( 2.0D+00 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
        q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )
      end if
      if ( 0.0D+00 < p ) then
        q = - q
      else
        p = - p
      end if
      s = e
      e = d
      if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. &
        p < abs ( 0.5D+00 * s * q ) ) then
        d = p / q
      else
        e = m
        d = e
      end if
    end if
    sa = sb
    fa = fb
    if ( tol < abs ( d ) ) then
      sb = sb + d
    else if ( 0.0D+00 < m ) then
      sb = sb + tol
    else
      sb = sb - tol
    end if
    call f_cold_eos_mhd_value(sb,udens,S2,B2,BS,ye,temp,fb)
    if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
         ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
      c = sa
      fc = fa
      e = sb - sa
      d = e
    end if
  end do
  x = sb !this returns x
end subroutine BrentMethodCold

subroutine f_cold_eos_mhd_value ( w_lorentz, udens, s2, B2, BS, ye, temp, fval )
  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, intent(inout)    :: w_lorentz
  CCTK_REAL, intent(in)       :: udens, s2, B2, BS, ye, temp
  CCTK_REAL, intent(out)   :: fval
  CCTK_REAL :: Z, rho, eps, press

  character(len=100) warnline

  !!!!!!
  ! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1),handle
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1)

  n = 1;keytemp = 0;anyerr = 0;keyerr = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = temp; xye = ye
  ! end EOS Omni vars
  !!!!!!

  w_lorentz = max(w_lorentz,1.0d0) !do not allow for unphysical W, would bring to negative density.
  rho = max(udens/w_lorentz,rho_abs_eos_min)

  xrho = rho
  xtemp = temp
  xye = ye
  handle  = 4
  keytemp = 1
  call EOS_Omni_press(handle,keytemp,Eos_rf_prec,n,&
                        xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
  press   = xpress(1)
  eps     = xeps(1)

  if(anyerr.ne.0) then
    !$OMP CRITICAL
    write(warnline,'(a20,g16.7)') 'udens= ',udens
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'w_lorentz= ',w_lorentz
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'eps= ',eps
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'temp= ',temp
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'ye= ',ye
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'keyerr= ',keyerr(1)
    call CCTK_WARN(1,warnline)
    call CCTK_WARN(0,"EOS_Omni error in f_cold_eos_mhd_value")
    !$OMP END CRITICAL
  end if

  Z    = (rho+rho*eps+press)*w_lorentz*w_lorentz
  fval = s2 - ((Z+B2)**2*(w_lorentz**2-1.0d0)/w_lorentz**2 - (2.0d0*Z+B2)*BS**2/Z**2)
end subroutine f_cold_eos_mhd_value



 /*@@
   @routine    Con2PrimPolytype_mhd
   @date
   @author     Bruno Giacomazzo
   @desc
   All routines below are identical to those above, just
   specialised from polytropic EOS.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/


subroutine Con2PrimPolytype_mhd(CCTK_ARGUMENTS)

  use Interfaces,only: Con2Prim_ptPolytype_mhd

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

!#ifdef _EOS_BASE_INC_
!#undef _EOS_BASE_INC_
!#endif
!#include "EOS_Base.inc"


  integer :: i, j, k, nx, ny, nz
  CCTK_REAL :: uxx, uxy, uxz, uyy, uyz, uzz, det, psi4pt
  CCTK_REAL :: Bxlow, Bylow, Bzlow, B2

  CCTK_INT :: k_lowerbound, k_upperbound, j_lowerbound, j_upperbound, i_lowerbound, i_upperbound
  !the following variables are used to convert the centered derivative in the appropriate one-sided derivative
  ! when computing B from the curl of A
  CCTK_INT :: iplus, iminus, jplus, jminus, kplus, kminus
  CCTK_INT :: iiplus, iiminus, jjplus, jjminus, kkplus, kkminus

  CCTK_REAL :: idx, idy, idz

  character(len=100) warnline

  CCTK_INT :: epsnegative
  CCTK_REAL :: eos_gamma, b2small_loc, eos_kappa, sc

  !!!!!!
  ! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1),one(1),ten(1)
  CCTK_REAL :: retvalcheck

  n = 1;keytemp = 0;anyerr = 0;keyerr = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
  one = 1.0d0;ten = 10.0d0
  retvalcheck = 0.0d0
  ! end EOS Omni vars
  !!!!!!

  !!*  In recovery section, Con2Prim shouldn't be executed.
  IF ( Recovery_Flag == 1 ) RETURN

  !trick to obtain the value of eos_ideal_fluid_gamma from thorn EOS_Ideal
  !eos_kappa = EOS_Pressure(Spritz_eos_handle,1.0d0,1.0d0)
  !p=k*rho**eos_gamma
  !if rho=10 then log10(k*rho**eos_gamma/k)=log10(10**eos_gamma)=eos_gamma
  !eos_gamma = log10(EOS_Pressure(Spritz_eos_handle,10.0d0,1.0d0)/eos_kappa)

  !!!!!!
  call EOS_Omni_press(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                 one,one,xtemp,xye,xpress,keyerr,anyerr)
  eos_kappa = xpress(1)
  call EOS_Omni_press(Spritz_eos_handle,keytemp,Eos_rf_prec,n,&
                 ten,one,xtemp,xye,xpress,keyerr,anyerr)
  eos_gamma = log10(xpress(1)/eos_kappa)
  !!!!!!

!!$  write(*,*) "eos_kappa=", eos_kappa
!!$  write(*,*) "eos_gamma=", eos_gamma

  b2small_loc=0.0d0 !initialize b2small_loc=b^mu*b_mu
  epsnegative=0

!!$  write(*,*) "Con2PrimPolytype_mhd"

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  idx = 1.0d0/CCTK_DELTA_SPACE(1)
  idy = 1.0d0/CCTK_DELTA_SPACE(2)
  idz = 1.0d0/CCTK_DELTA_SPACE(3)

  iplus  = 1
  iminus = 1
  jplus  = 1
  jminus = 1
  kplus  = 1
  kminus = 1

  !Initialize the mask to zero
  mhd_c2p_mask=0


  if(Spritz_Avec_Staggered.ne.0) then
     k_lowerbound = 2 !whisky_stencil + 1
     k_upperbound = nz !nz - whisky_stencil

     j_lowerbound = 2 !whisky_stencil + 1
     j_upperbound = ny !ny - whisky_stencil

     i_lowerbound = 2 !whisky_stencil + 1
     i_upperbound = nx !nx - whisky_stencil
  else
     k_lowerbound = 1 !2 !Spritz_stencil + 1
     k_upperbound = nz !nz-1 !nz - Spritz_stencil

     j_lowerbound = 1 !2 !Spritz_stencil + 1
     j_upperbound = ny !ny-1 !ny - Spritz_stencil

     i_lowerbound = 1 !2 !Spritz_stencil + 1
     i_upperbound = nx !nx-1 !nx - Spritz_stencil
  end if


  !set MHD variables to be zero in order to check that they are correctly recomputed
  dBx = 0.0d0
  dBy = 0.0d0
  dBz = 0.0d0
  Bx = 0.0d0
  By = 0.0d0
  Bz = 0.0d0





  !$OMP PARALLEL DO PRIVATE(k, j, i, iminus, idx, jminus, idy, kminus, idz, iplus, jplus, kplus, &
  !$OMP iiplus, iiminus, jjplus, jjminus, kkplus, kkminus, &
  !$OMP warnline, uxx, uxy, uxz, uyy, uyz, uzz, det, psi4pt, &
  !$OMP  Bxlow, Bylow, Bzlow, B2, b2small_loc, retvalcheck)
  do k = k_lowerbound, k_upperbound  !!!1, nz
     do j = j_lowerbound, j_upperbound !!!1, ny
        do i = i_lowerbound, i_upperbound !!!1, nx

           !write(infoline,'(a25,g20.7,a8,g20.7,a8,g20.7)') ' C2P input values: dBx = ',dBx(i,j,k), &
           !                                                ', dBy = ',dBy(i,j,k),', dBz = ',dBz(i,j,k)
           !call CCTK_INFO(infoline)
           
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

	
           idx = 1.0d0/CCTK_DELTA_SPACE(1)	
           idy = 1.0d0/CCTK_DELTA_SPACE(2)
           idz = 1.0d0/CCTK_DELTA_SPACE(3)


           ! The following block is for setting: !!!
           ! dBx(i,1,k) = dBx(i,2,k)             !!!
           ! dBy(i,j,1) = dBy(i,j,2)             !!!
           ! dBz(1,j,k) = dBz(2,j,k)             !!!
                                                 !!!
           if (i .eq. 1) then                    !!!
              iiminus = 0                        !!!
              iiplus = 1                         !!!
           else                                  !!!
              iiminus = 1                        !!!
              iiplus = 0                         !!!
           end if                                !!!
                                                 !!!
           if (j .eq. 1) then                    !!!
              jjminus = 0                        !!!
              jjplus = 1                         !!!
           else                                  !!!
              jjminus = 1                        !!!
              jjplus = 0                         !!!
           end if                                !!!
                                                 !!!
           if (k .eq. 1) then                    !!!
              kkminus = 0                        !!!
              kkplus = 1                         !!!
           else                                  !!!
              kkminus = 1                        !!!
              kkplus = 0                         !!!
           end if                                !!!
                                                 !!!
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           if(Spritz_Avec_Staggered.ne.0) then
              !dBx is 0.5*(curl(A) at i+1/2 + curl(A) at i-1/2)
              dBx(i,j,k) = 0.5d0 * ( (idy*(Avec_z(i+iiplus,j+jjplus,k+kkplus)-Avec_z(i+iiplus,j-jjminus,k+kkplus)) - &
                                      idz*(Avec_y(i+iiplus,j+jjplus,k+kkplus)-Avec_y(i+iiplus,j+jjplus,k-kkminus))) + &
                                     (idy*(Avec_z(i-iiminus,j+jjplus,k+kkplus)-Avec_z(i-iiminus,j-jjminus,k+kkplus)) - &
                                      idz*(Avec_y(i-iiminus,j+jjplus,k+kkplus)-Avec_y(i-iiminus,j+jjplus,k-kkminus))) )
              !dBy is 0.5*(curl(A) at j+1/2 + curl(A) at j-1/2)
              dBy(i,j,k) = 0.5d0 * ( (idz*(Avec_x(i+iiplus,j+jjplus,k+kkplus)-Avec_x(i+iiplus,j+jjplus,k-kkminus)) - &
                                      idx*(Avec_z(i+iiplus,j+jjplus,k+kkplus)-Avec_z(i-iiminus,j+jjplus,k+kkplus))) + &
                                     (idz*(Avec_x(i+iiplus,j-jjminus,k+kkplus)-Avec_x(i+iiplus,j-jjminus,k-kkminus)) - &
                                      idx*(Avec_z(i+iiplus,j-jjminus,k+kkplus)-Avec_z(i-iiminus,j-jjminus,k+kkplus))) )
              !dBz is 0.5*(curl(A) at k+1/2 + curl(A) at k-1/2)
              dBz(i,j,k) = 0.5d0 * ( (idx*(Avec_y(i+iiplus,j+jjplus,k+kkplus)-Avec_y(i-iiminus,j+jjplus,k+kkplus)) - &
                                      idy*(Avec_x(i+iiplus,j+jjplus,k+kkplus)-Avec_x(i+iiplus,j-jjminus,k+kkplus))) + &
                                     (idx*(Avec_y(i+iiplus,j+jjplus,k-kkminus)-Avec_y(i-iiminus,j+jjplus,k-kkminus)) - &
                                      idy*(Avec_x(i+iiplus,j+jjplus,k-kkminus)-Avec_x(i+iiplus,j-jjminus,k-kkminus))) )
           else
              !reset and initialize the constants before starting
              !(these are private variables in OpenMP, so we need to make sure that
              ! they are correctly initialized inside the do loop. Otherwise you may get segfaults
              ! when using optimization, because of out of boundaries errors)
              ! (declaring them as firstprivate would probably also work)
              idx = 1.0d0/CCTK_DELTA_SPACE(1)
              idy = 1.0d0/CCTK_DELTA_SPACE(2)
              idz = 1.0d0/CCTK_DELTA_SPACE(3)
              iplus  = 1
              iminus = 1
              jplus  = 1
              jminus = 1
              kplus  = 1
              kminus = 1

              !check if we have enough point to use the standard second order centered derivative or if we need to use a one-sided derivative
              if (i==1) then
                 iminus=0
                 idx = 2.0*idx
              end if
              if (j==1) then
                 jminus=0
                 idy = 2.0*idy
              end if
              if (k==1) then
                 kminus=0
                 idz = 2.0*idz
              end if
              if (i==nx) then
                 iplus=0
                 idx = 2.0*idx
              end if
              if (j==ny) then
                 jplus=0
                 idy = 2.0*idy
              end if
              if (k==nz) then
                 kplus=0
                 idz = 2.0*idz
              end if

              !Now we compute sqrt(gamma)*B from the curl of the vector potential A
              !dBx is curl(A)
              dBx(i,j,k) = 0.5d0*(idy*(Avec_z(i,j+jplus,k)-Avec_z(i,j-jminus,k)) - idz*(Avec_y(i,j,k+kplus)-Avec_y(i,j,k-kminus)))
              !dBy is curl(A)
              dBy(i,j,k) = 0.5d0*(idz*(Avec_x(i,j,k+kplus)-Avec_x(i,j,k-kminus)) - idx*(Avec_z(i+iplus,j,k)-Avec_z(i-iminus,j,k)))
              !dBz is curl(A)
              dBz(i,j,k) = 0.5d0*(idx*(Avec_y(i+iplus,j,k)-Avec_y(i-iminus,j,k)) - idy*(Avec_x(i,j+jplus,k)-Avec_x(i,j-jminus,k)))
           end if

           !BRUNO DEBUG STUFF
           if((dBx(i,j,k).ne.dBx(i,j,k)).or.(dBy(i,j,k).ne.dBy(i,j,k)).or.(dBz(i,j,k).ne.dBz(i,j,k))) then
              !$OMP CRITICAL
              write(warnline,'(a20,g16.7)') 'dBx= ',dBx(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'dBy= ',dBy(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'dBz= ',dBz(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'Avec_x= ',Avec_x(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'Avec_y= ',Avec_y(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'Avec_z= ',Avec_z(i,j,k)
              call CCTK_WARN(1,warnline)

              write(warnline,'(a20,g16.7)') 'Avec_x(j+jplus)= ',Avec_x(i,j+jplus,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'Avec_x(j-jminus)= ',Avec_x(i,j-jminus,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'Avec_x(k+kplus)= ',Avec_x(i,j,k+kplus)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'Avec_x(k-kminus)= ',Avec_x(i,j,k-kminus)
              call CCTK_WARN(1,warnline)

              write(warnline,'(a20,g16.7)') 'Avec_z(j+jplus)= ',Avec_z(i,j+jplus,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'Avec_z(j-jminus)= ',Avec_z(i,j-jminus,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'Avec_z(i+iplus)= ',Avec_z(i+iplus,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'Avec_z(i-iminus)= ',Avec_z(i-iminus,j,k)
              call CCTK_WARN(1,warnline)

              write(warnline,'(a20,g16.7)') 'Avec_y(i+iplus)= ',Avec_y(i+iplus,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'Avec_y(i-iminus)= ',Avec_y(i-iminus,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'Avec_y(k+kplus)= ',Avec_y(i,j,k+kplus)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'Avec_y(k-kminus)= ',Avec_y(i,j,k-kminus)
              call CCTK_WARN(1,warnline)

              write(warnline,'(a20,g16.7)') 'iplus= ',iplus
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'iminus= ',iminus
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'jplus= ',jplus
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'jminus= ',jminus
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'kplus= ',kplus
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'kminus= ',kminus
              call CCTK_WARN(1,warnline)

              write(warnline,'(a20,g16.7)') 'x= ',x(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'y= ',y(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'z= ',z(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'r= ',r(i,j,k)
              call CCTK_WARN(1,warnline)

              call CCTK_WARN(0,"The magnetic field is a NaN in C2P Poly")
              !$OMP END CRITICAL
           end if


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



           !BRUNO DEBUG STUFF
           if((det.ne.det).or.(det<0.0d0)) then
              !$OMP CRITICAL
              write(warnline,'(a20,g16.7)') 'det= ',det
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'x= ',x(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'y= ',y(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'z= ',z(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g16.7)') 'r= ',r(i,j,k)
              call CCTK_WARN(1,warnline)
              call CCTK_WARN(0,"The determinant of the 3-metric is NaN or negative in C2P poly")
              !$OMP END CRITICAL
           end if


           !B^j
           Bx(i,j,k) = dBx(i,j,k)/sqrt(det)
           By(i,j,k) = dBy(i,j,k)/sqrt(det)
           Bz(i,j,k) = dBz(i,j,k)/sqrt(det)

           !if (Spritz_mhd_HypDiv_on .ne. 0) psi_mhd(i,j,k) = dpsi_mhd(i,j,k)/sqrt(det)

           !B_j
           Bxlow = Bx(i,j,k)*gxx(i,j,k)+By(i,j,k)*gxy(i,j,k)+Bz(i,j,k)*gxz(i,j,k)
           Bylow = Bx(i,j,k)*gxy(i,j,k)+By(i,j,k)*gyy(i,j,k)+Bz(i,j,k)*gyz(i,j,k)
           Bzlow = Bx(i,j,k)*gxz(i,j,k)+By(i,j,k)*gyz(i,j,k)+Bz(i,j,k)*gzz(i,j,k)

           !B^j B_j
           B2 = gxx(i,j,k) * Bx(i,j,k)**2 + &
                gyy(i,j,k) * By(i,j,k)**2 + &
                gzz(i,j,k) * Bz(i,j,k)**2 + &
                2.d0 * gxy(i,j,k) * Bx(i,j,k) * By(i,j,k) + &
                2.d0 * gxz(i,j,k) * Bx(i,j,k) * Bz(i,j,k) + &
                2.d0 * gyz(i,j,k) * By(i,j,k) * Bz(i,j,k)



                      !!CHECK ON THE ATMO DONE BEFORE THE C2P Netwon-Raphson routine
           !!THIS IS DONE ONLY FOR COMPUTATIONAL EFFICIENCY
           if(dens(i,j,k)<=(sqrt(det)*(Spritz_rho_min*(1.d0+Spritz_atmo_tolerance))).or.(alp(i,j,k)<Spritz_alp_BH)) then
              !atmosphere
              rho(i,j,k) = Spritz_rho_min
              dens(i,j,k) = sqrt(det)*rho(i,j,k)
              !Since I don't know what to do in this case for eps and press
              !I set them using a polytropic EOS which means that the atmosphere is always cold
              
              !!!!!!
              xrho = rho(i,j,k)
              xeps = eps(i,j,k)
              call EOS_Omni_press(Spritz_polytrope_handle,keytemp,Eos_rf_prec,n,&
                             xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
              press(i,j,k) = xpress(1)
              call EOS_Omni_EpsFromPress(Spritz_polytrope_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
              eps(i,j,k) = xeps(1)
              !!!!!!

              scon(i,j,k,1) = 0.0d0
              scon(i,j,k,2) = 0.0d0
              scon(i,j,k,3) = 0.0d0
              vel(i,j,k,1) = 0.0d0
              vel(i,j,k,2) = 0.0d0
              vel(i,j,k,3) = 0.0d0
              w_lorentz(i,j,k) = 1.0d0
              !idem here
              tau(i,j,k)  = sqrt(det) * ((rho(i,j,k)+rho(i,j,k)*eps(i,j,k) + &
                   press(i,j,k)+B2)*w_lorentz(i,j,k)*w_lorentz(i,j,k) - (press(i,j,k)+0.5d0*B2)) - dens(i,j,k)
              cycle
           end if


           if (CCTK_EQUALS(con2prim_mhd_method,"3eqs")) then
              call Con2Prim_ptPolytype_mhd(Spritz_eos_handle, dens(i,j,k),&
                   scon(i,j,k,1),scon(i,j,k,2), &
                   scon(i,j,k,3), dBx(i,j,k), dBy(i,j,k), dBz(i,j,k), &
                   tau(i,j,k),rho(i,j,k),vel(i,j,k,1),vel(i,j,k,2), &
                   vel(i,j,k,3),eps(i,j,k),press(i,j,k),w_lorentz(i,j,k), &
                   uxx,uxy,uxz,uyy,uyz,uzz,det,x(i,j,k),y(i,j,k), &
                   z(i,j,k),r(i,j,k),Spritz_rho_min, &
                   Bx(i,j,k), By(i,j,k), Bz(i,j,k), B2, Bxlow, Bylow, Bzlow, Spritz_reflevel)
           elseif (CCTK_EQUALS(con2prim_mhd_method,"Noble")) then
              !! In principle the same subroutine in Spritz_Con2PrimM_Noble_pt_EOSOmni.c
              !! allows the implementation of both polytropic and general tabulated EOSs...
              call Spritz_Con2Prim_Noble(Spritz_eos_handle,Spritz_reflevel, i, j, k, x(i,j,k), y(i,j,k), &
                   z(i,j,k), keytemp, Eos_rf_prec, Con2Prim_perc_ptol, eos_gamma, dens(i,j,k), &
                   scon(i,j,k,1),scon(i,j,k,2),scon(i,j,k,3), tau(i,j,k), &
                   dBx(i,j,k), dBy(i,j,k), dBz(i,j,k), xye, xtemp, rho(i,j,k),&
                   vel(i,j,k,1),vel(i,j,k,2),vel(i,j,k,3),eps(i,j,k),press(i,j,k),&
                   Bx(i,j,k), By(i,j,k), Bz(i,j,k), b2small_loc, w_lorentz(i,j,k),&
                   gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k), &
                   uxx,uxy,uxz,uyy,uyz,uzz,det, &
                   epsnegative,Spritz_rho_min,retvalcheck)        !mhd_c2p_mask(i,j,k)
              if(abs(retvalcheck) .gt. 1.0d-8 ) then
                               !$OMP CRITICAL
                               write(warnline,'(a80,g20.7)') 'ATTENTION: Problems in C2P_Nobel, polytropic case - returned Retval = ',retvalcheck
                               call CCTK_WARN(2,warnline)
                               !$OMP END CRITICAL
              end if
           else
              !$OMP CRITICAL
              call CCTK_WARN(0,"MHD C2P Polytyope allows only con2prim_mhd_method equal to 3eqs or Noble")
              !$OMP END CRITICAL
           end if

           if((rho(i,j,k).ne.rho(i,j,k)).or. &
                (eps(i,j,k).ne.eps(i,j,k)).or. &
                (vel(i,j,k,1).ne.vel(i,j,k,1)).or. &
                (vel(i,j,k,2).ne.vel(i,j,k,2)).or. &
                (vel(i,j,k,3).ne.vel(i,j,k,3)).or. &
                (Bx(i,j,k).ne.Bx(i,j,k)).or. &
                (By(i,j,k).ne.By(i,j,k)).or. &
                (Bz(i,j,k).ne.Bz(i,j,k))) then  !.or.c2perror
              !$OMP CRITICAL
              call CCTK_WARN(2,"Failure or NaN in C2P_MHD Ideal Fluid")
              write(warnline,'(a20,g20.7)') ' rho is NaN: ',rho(i,j,k)
              call CCTK_WARN(2,warnline)
              write(warnline,'(a20,g20.7)') ' eps is NaN: ',eps(i,j,k)
              call CCTK_WARN(2,warnline)
              write(warnline,'(a20,g20.7)') 'velx is NaN: ',vel(i,j,k,1)
              call CCTK_WARN(2,warnline)
              write(warnline,'(a20,g20.7)') 'vely is NaN: ',vel(i,j,k,2)
              call CCTK_WARN(2,warnline)
              write(warnline,'(a20,g20.7)') 'velz is NaN: ',vel(i,j,k,3)
              call CCTK_WARN(2,warnline)
              write(warnline,'(a20,g20.7)') '  Bx is NaN: ',Bx(i,j,k)
              call CCTK_WARN(2,warnline)
              write(warnline,'(a20,g20.7)') '  By is NaN: ',By(i,j,k)
              call CCTK_WARN(2,warnline)
              write(warnline,'(a20,g20.7)') '  Bz is NaN: ',Bz(i,j,k)
              call CCTK_WARN(2,warnline)
              write(warnline,'(a20,g20.7)') 'x=',x(i,j,k)
              call CCTK_WARN(2,warnline)
              write(warnline,'(a20,g20.7)') 'y=',y(i,j,k)
              call CCTK_WARN(2,warnline)
              write(warnline,'(a20,g20.7)') 'z=',z(i,j,k)
              call CCTK_WARN(2,warnline)
              call CCTK_WARN(2, 'setting the point to atmosphere')
              !$OMP END CRITICAL
              
              !set this mask to 1 to identify the failure
              !The above message is set to warning level 2 and so it would not be shown usually
              mhd_c2p_mask(i,j,k)=1

              !atmosphere
              rho(i,j,k) = Spritz_rho_min
              dens(i,j,k) = sqrt(det)*rho(i,j,k)
              !Since I don't know what to do in this case for eps and press
              !I set them using a polytropic EOS which means that the atmosphere is always cold
              !press(i,j,k) = EOS_Pressure(Spritz_polytrope_handle, rho(i,j,k), eps(i,j,k))
              !eps(i,j,k) = EOS_SpecificIntEnergy(Spritz_polytrope_handle,rho(i,j,k),press(i,j,k))

              !!!!!!
              xrho = rho(i,j,k)
              xeps = eps(i,j,k)
              call EOS_Omni_press(Spritz_polytrope_handle,keytemp,Eos_rf_prec,n,&
                             xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
              press(i,j,k) = xpress(1)
              call EOS_Omni_EpsFromPress(Spritz_polytrope_handle,keytemp,Eos_rf_prec,n,&
                                    xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
              eps(i,j,k) = xeps(1)
              !!!!!!

              scon(i,j,k,1) = 0.0d0
              scon(i,j,k,2) = 0.0d0
              scon(i,j,k,3) = 0.0d0
              vel(i,j,k,1) = 0.0d0
              vel(i,j,k,2) = 0.0d0
              vel(i,j,k,3) = 0.0d0
              w_lorentz(i,j,k) = 1.0d0
              !idem here
              tau(i,j,k)  = sqrt(det) * ((rho(i,j,k)+rho(i,j,k)*eps(i,j,k) + &
                   press(i,j,k)+B2)*w_lorentz(i,j,k)*w_lorentz(i,j,k) - (press(i,j,k)+0.5d0*B2)) - dens(i,j,k)
           end if

        end do
     end do
  end do
  !$OMP END PARALLEL DO

!!$write(*,*) 'Con2PrimPolytype_mhd OK'

end subroutine Con2PrimPolytype_mhd

/*@@
@routine    Con2Prim_ptPolytype_mhd
@date
@author    Bruno Giacomazzo
@desc
Given all the appropriate data, recover the primitive variables at
a single point.
@enddesc
 @calls
@calledby
@history

@endhistory

@@*/

subroutine Con2Prim_ptPolytype_mhd(                       &
     & handle,                                            &
     & dens, sx, sy, sz,                                  &
     & dBx, dBy, dBz,                                     &
     & tau, rho, velx, vely, velz, eps, press, w_lorentz, &
     & uxx, uxy, uxz, uyy, uyz, uzz, det,                 &
     & x, y, z, r, Spritz_rho_min,                        &
     & Bx, By, Bz, B2, Bxlow, Bylow, Bzlow,               &
     & Spritz_reflevel, c2perror, mhd_c2p_mask_local )

  use Interfaces, only: jacobian_mhd_polytype
  use nrmodule
  implicit none

  DECLARE_CCTK_PARAMETERS

  CCTK_REAL:: dens, sx, sy, sz, dBx, dBy, dBz, tau, rho, velx, vely, velz, eps, &
       press, uxx, uxy, uxz, uyy, uyz, uzz, det, w_lorentz, x, &
       y, z, r, Spritz_rho_min, Bx, By, Bz, B2, Bxlow, Bylow, Bzlow, &
       udens, usx, usy, usz, s2, w_lorentzold, BS, bst, bs2
  CCTK_INT :: i, handle, Spritz_reflevel, error_level
  CCTK_REAL:: vxlow, vylow, vzlow, w, ww2
  CCTK_REAL, optional:: mhd_c2p_mask_local
  logical, optional::c2perror

  CCTK_REAL :: w_lorentz_min, w_lorentz_max, fl, fh, f, df, dfl, xl, xh,&
       dx, dxold, temp

  character(len=100) warnline

!#ifdef _EOS_BASE_INC_
!#undef _EOS_BASE_INC_
!#endif
!#include "EOS_Base.inc"

  !!!!!!
  ! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1),zero(1)

  n = 1;keytemp = 0;anyerr = 0;keyerr = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
  zero = 0.0d0
  ! end EOS Omni vars
  !!!!!

  if(present(c2perror)) c2perror=.false.


!!$  Undensitize the variables
  udens = dens /sqrt(det)
  usx = sx /sqrt(det)
  usy = sy /sqrt(det)
  usz = sz /sqrt(det)
  s2 = usx*usx*uxx + usy*usy*uyy + usz*usz*uzz + 2.d0*usx*usy*uxy + &
       2.d0*usx*usz*uxz + 2.d0*usy*usz*uyz

  !B^j S_j
  BS = Bx*usx + By*usy + Bz*usz


  !!$ if dens is already below the atmosphere level, then no Newton-Raphson
  if(dens<=(sqrt(det)*Spritz_rho_min*(1.d0+Spritz_atmo_tolerance))) then !NB in Spritz_Dev there is no factor sqrt(det) here
     !atmosphere
     rho = Spritz_rho_min
     w_lorentz = 1.0d0
     sx = 0.d0
     sy = 0.d0
     sz = 0.d0
     s2 = 0.d0
     usx = 0.d0
     usy = 0.d0
     usz = 0.d0
     BS=0.0d0
     goto 1001
  end if


!!$  Set initial guess

!!$  rhoold = max(Spritz_rho_min,rho)
  w_lorentzold = max(1.0d0,w_lorentz)

  w_lorentz = w_lorentzold

  w_lorentz_min = 1.0d0
  w_lorentz_max = 1.0d4

!!$  !USED FOR DEBUG
!!$  B2=10.0d0
!!$  s2=1.0d4
!!$  BS=1.0d3
!!$  udens=1.0d0
!!$  open(unit=11,file='func_1eq.dat',status='unknown')
!!$  w_lorentz = w_lorentz_min
!!$  do i=1,100
!!$     call  jacobian_mhd_polytype(handle,udens,s2,B2,BS,w_lorentz,fl,dfl)
!!$     write(*,*) w_lorentz,fl
!!$     write(11,FMT="(3g17.4)") w_lorentz,fl,dfl
!!$     w_lorentz = w_lorentz + 1
!!$  end do
!!$  close(11)
!!$  call CCTK_WARN(0,"results in func_1eq.dat")


  call jacobian_mhd_polytype(handle,udens,s2,B2,BS,w_lorentz_min,fl,dfl)
  if ((fl == 0.0).and.(abs(fl/dfl)<=Spritz_mhd_tolf)) then
     w_lorentz=w_lorentz_min
     !$OMP CRITICAL
     write(warnline,'(a20,g16.7)') 'w_lorentz_min ',w_lorentz_min
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g16.7)') 'fl ',fl
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g16.7)') 'df ',dfl
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g16.7)') 'fl/df ',fl/dfl
     call CCTK_WARN(2,warnline)
     call CCTK_WARN(2,"Using w_lorentz_min as root")
     !$OMP END CRITICAL
     goto 2000
  end if

5 call jacobian_mhd_polytype(handle,udens,s2,B2,BS,w_lorentz_max,fh,df)
  if ((fh == 0.0).and.(abs(fh/df)<=Spritz_mhd_tolf)) then
     !$OMP CRITICAL
     write(warnline,'(a20,g16.7)') 'w_lorentz_max ',w_lorentz_max
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g16.7)') 'fh ',fh
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g16.7)') 'df ',df
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,g16.7)') 'fh/df ',fh/df
     call CCTK_WARN(2,warnline)
     call CCTK_WARN(2,"Using w_lorentz_max as root")
     !$OMP END CRITICAL
     w_lorentz=w_lorentz_max
     goto 2000
  end if

  if ((fl > 0.0 .and. fh > 0.0) .or.(fl < 0.0 .and. fh < 0.0)) then
     if (w_lorentz_max > 1.0d20) then
        !$OMP CRITICAL
        write(warnline,'(a3,g16.7,a3,g16.7,a3,g16.7,a3,g16.7)') 'x=',x,'y=',y,'z=',z,'r=',r
        call CCTK_WARN(1,warnline)
        write(warnline,'(a15,g16.7,a15,g16.7)') 'w_lorentz_min=',w_lorentz_min,'w_lorentz_max=',w_lorentz_max
        call CCTK_WARN(1,warnline)
        write(warnline,'(a4,g16.7,a4,g16.7)') 'fl=',fl,'fh=',fh
        call CCTK_WARN(1,warnline)
        write(warnline,'(a7,g16.7,a7,g16.7)') 'fl/df=',fl/dfl,'fh/df=',fh/df
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g16.7)') 'de is ',dens
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g16.7)') 'sx is ',sx
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g16.7)') 'sy is ',sy
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g16.7)') 'sz is ',sz
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g16.7)') 'en is ',tau
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g16.7)') 'Bx: ',Bx
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g16.7)') 'By: ',By
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g16.7)') 'Bz: ',Bz
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g16.7)') 's2: ',s2
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g16.7)') 'B2: ',B2
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g16.7)') 'udens: ',udens
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g16.7)') 'BS: ',BS
        call CCTK_WARN(1,warnline)
        call CCTK_WARN(0,"root not bracketed in rtsafe!")
        !$OMP END CRITICAL
     end if
     w_lorentz_max=w_lorentz_max*1.5d0
     goto 5
  end if

  if (fl < 0.0) then !Orient the search so that f(xl)<0.
     xl=w_lorentz_min
     xh=w_lorentz_max
  else
     xh=w_lorentz_min
     xl=w_lorentz_max
  end if
!!$  w_lorentz=0.5d0*(w_lorentz_min+w_lorentz_max) !Initialize the guess for root,
  dxold=abs(w_lorentz_max-w_lorentz_min) !the "stepsize before last",
  dx=dxold !and the last step.
  call jacobian_mhd_polytype(handle,udens,s2,B2,BS,w_lorentz,f,df)
  do i=1, Spritz_countmax !Loop over allowed iterations.
     if (((w_lorentz-xh)*df-f)*((w_lorentz-xl)*df-f) > 0.0 .or. &
          abs(2.0d0*f) > abs(dxold*df) ) then
        !Bisect if Newton out of range, or not decreasing fast enough.
        dxold=dx
        dx=0.5d0*(xh-xl)
        w_lorentz=xl+dx
        if (xl == w_lorentz) goto 2000 !Change in root is negligible.
     else !Newton step acceptable. Take it.
        dxold=dx
        dx=f/df
        temp=w_lorentz
        w_lorentz=w_lorentz-dx
        if (temp == w_lorentz) goto 2000
     end if
     if (abs(dx/w_lorentz) < Spritz_mhd_tolf) goto 2000 !Convergence criterion. (RELATIVE ERROR)
     call jacobian_mhd_polytype(handle,udens,s2,B2,BS,w_lorentz,f,df)
     !if (f==0.0d0) goto 2000 !Convergence criterion.

!!$     if ((x==15.36d0).and.(y==-3.84d0).and.(z==-7.68d0)) then
!!$        write(warnline,'(a4,g16.7,a4,g16.7,a4,g16.7,a4,g16.7,a4,g16.7)') 'i=',i,'W=',w_lorentz,'f=',f,'df=',df,'dx=',dx!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$        call CCTK_WARN(1,warnline)
!!$     end if

     !One new function evaluation per iteration.
     if (f < 0.0) then !Maintain the bracket on the root.
        xl=w_lorentz
     else
        xh=w_lorentz
     end if
  end do

  ! print warnings only on a level for which this is required
  if (Spritz_reflevel.ge.Spritz_c2p_warn_from_reflevel) then
     error_level=1 !set warning level to 1 to print message to screen
     if(present(c2perror)) error_level=2
     !$OMP CRITICAL
     call CCTK_WARN(error_level, 'Con2Prim_ptPolytype_mhd: error: did not converge in ')
     write(warnline,'(a20,g16.7,a10)') '              ',&
          Spritz_countmax,' steps'
     call CCTK_WARN(error_level,warnline)
     write(warnline,'(a28,i2)') 'on carpet reflevel: ',Spritz_reflevel
     call CCTK_WARN(error_level,warnline)
     write(warnline,'(a20,g16.7)') 'de is ',dens
     call CCTK_WARN(error_level,warnline)
     write(warnline,'(a20,g16.7)') 'sx is ',sx
     call CCTK_WARN(error_level,warnline)
     write(warnline,'(a20,g16.7)') 'sy is ',sy
     call CCTK_WARN(error_level,warnline)
     write(warnline,'(a20,g16.7)') 'sz is ',sz
     call CCTK_WARN(error_level,warnline)
     write(warnline,'(a20,g16.7)') 'en is ',tau
     call CCTK_WARN(error_level,warnline)
     write(warnline,'(a20,g16.7)') 'func value: ',f
     call CCTK_WARN(error_level,warnline)
     write(warnline,'(a20,g16.7)') 'Jacobian: ',df
     call CCTK_WARN(error_level,warnline)
     write(warnline,'(a20,g16.7)') 'Wnew: ',w_lorentz
     call CCTK_WARN(error_level,warnline)
     write(warnline,'(a20,g16.7)') 'Wold: ',w_lorentzold
     call CCTK_WARN(error_level,warnline)
     write(warnline,'(a20,g16.7)') 'correction: ',dx
     call CCTK_WARN(error_level,warnline)
     write(warnline,'(a20,g16.7)') 'epsilon: ',eps
     call CCTK_WARN(error_level,warnline)
     write(warnline,'(a20,g16.7)') 'Bx: ',Bx
     call CCTK_WARN(error_level,warnline)
     write(warnline,'(a20,g16.7)') 'By: ',By
     call CCTK_WARN(error_level,warnline)
     write(warnline,'(a20,g16.7)') 'Bz: ',Bz
     call CCTK_WARN(error_level,warnline)
     write(warnline,'(a20,g16.7)') 's2: ',s2
     call CCTK_WARN(error_level,warnline)
     write(warnline,'(a20,g16.7)') 'B2: ',B2
     call CCTK_WARN(error_level,warnline)
     write(warnline,'(a20,g16.7)') 'udens: ',udens
     call CCTK_WARN(error_level,warnline)
     write(warnline,'(a20,g16.7)') 'BS: ',BS
     call CCTK_WARN(error_level,warnline)
     write(warnline,'(a20,3g16.7)') 'xyz location: ',x,y,z
     call CCTK_WARN(error_level,warnline)
     write(warnline,'(a20,g16.7)') 'radius: ',r
     call CCTK_WARN(error_level,warnline)
     !$OMP END CRITICAL
     if(present(c2perror)) then
        c2perror=.true.
        return
     end if
     !$OMP CRITICAL
     call CCTK_WARN(Spritz_c2p_warnlevel, "Did not converge")
     !$OMP END CRITICAL
  else
     !$OMP CRITICAL
     write(warnline,'(a60,i2)') 'Con2Prim: did not converge, but I was told to ignore level: ',Spritz_reflevel
     call CCTK_WARN(2,warnline)!!!!!!!!!!!!!!!!!!
     call CCTK_WARN(2,'For security reason I reset the point to be atmosphere to avoid NaNs')
     !$OMP END CRITICAL
     !atmosphere
     rho = Spritz_rho_min
     udens = rho
     dens = sqrt(det)*rho
     w_lorentz = 1.0d0
     sx = 0.d0
     sy = 0.d0
     sz = 0.d0
     s2 = 0.d0
     usx = 0.d0
     usy = 0.d0
     usz = 0.d0
     velx = 0.0d0
     vely = 0.0d0
     velz = 0.0d0
     BS=0.0d0
     B2 = (Bx*Bxlow + By*Bylow + Bz*Bzlow)!this is true only if v=0
     !Since I don't know what to do in this case for eps and rho
     !I set them using a polytropic EOS...
     !(I don't know if this can create problems...
     !...but look here if you search for a bug in con2prim)
     !press= EOS_Pressure(handle, rho, eps)
     !eps = EOS_SpecificIntEnergy(handle,rho,press)

     !!!!!!
     xrho = rho
     xeps = eps
     call EOS_Omni_press(handle,keytemp,Eos_rf_prec,n,&
                    xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
     press = xpress(1)
     call EOS_Omni_EpsFromPress(handle,keytemp,Eos_rf_prec,n,&
                           xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     eps = xeps(1)
     !!!!!!

     !idem here
     tau  = sqrt(det) * ((rho+rho*eps+press+B2)*w_lorentz*w_lorentz - (press+0.5d0*B2)) - dens
     return
  end if

2000 rho = udens/w_lorentz
  !rho = unk(1)
  !w_lorentz = unk(2)


  if (rho .le. Spritz_rho_min*(1.d0+Spritz_atmo_tolerance)) then
     !atmosphere
     rho = Spritz_rho_min
     !eps = EOS_SpecificIntEnergy(handle,rho,0.0d0)
     !press = EOS_Pressure(handle,rho,0.0d0)
     !velx = 0.0d0
     !vely = 0.0d0
     !velz = 0.0d0
     w_lorentz = 1.0d0
     sx = 0.d0
     sy = 0.d0
     sz = 0.d0
     s2 = 0.d0
     usx = 0.d0
     usy = 0.d0
     usz = 0.d0
     BS=0.0d0
  end if



!!$  Compute velocities (these are the v_j)
  !relativistic enthalpy
!1001 press = EOS_Pressure(handle,rho,0.0d0)
!  eps = EOS_SpecificIntEnergy(handle,rho,0.0d0)

   !!!!!!
1001 CONTINUE
   xrho = rho
   call EOS_Omni_press(handle,keytemp,Eos_rf_prec,n,&
                  xrho,zero,xtemp,xye,xpress,keyerr,anyerr)
   press = xpress(1)
   call EOS_Omni_EpsFromPress(handle,keytemp,Eos_rf_prec,n,&
                         xrho,xeps,xtemp,xye,zero,xeps,keyerr,anyerr)
   eps = xeps(1)
   !!!!!!

  w = rho + rho*eps + press

  ww2 = w*w_lorentz**2 !this is Z = rho*h*W_lorentz**2 (Anton et al. 2005)

  vxlow = Bxlow*BS/(ww2*(ww2 + B2)) + usx/(ww2 + B2)
  vylow = Bylow*BS/(ww2*(ww2 + B2)) + usy/(ww2 + B2)
  vzlow = Bzlow*BS/(ww2*(ww2 + B2)) + usz/(ww2 + B2)


  !and now compute the v^j
  velx = vxlow*uxx + vylow*uxy + vzlow*uxz
  vely = vxlow*uxy + vylow*uyy + vzlow*uyz
  velz = vxlow*uxz + vylow*uyz + vzlow*uzz



  !We need to compute tau since it's not evolved in the polytropic case
  !If tau is not computed than you may experience very high mass losses
  !This is due to the UpdateMask routine which checks if tau+taurhs<0 and
  ! if so sets the point in which it happens to be atmosphere.
  ! Tau is not evolved when using a polytropic EOS,
  ! but taurhs is anyway computed (the routine
  ! was written to work with any EOS). So if you don't update tau then you
  ! can get tau+taurhs<0 even in a point which is not supposed to be
  ! atmosphere.

  bst = w_lorentz*(Bx*vxlow + By*vylow + Bz*vzlow) !this is b^0*alp
  bs2 = (Bx*Bxlow + By*Bylow + Bz*Bzlow + bst**2)/w_lorentz**2

  tau  = sqrt(det) * ((w+bs2)*w_lorentz*w_lorentz - (press+0.5d0*bs2) - &
       bst**2) - dens



  if((vxlow*velx + vylow*vely + vzlow*velz)>=1.0d0) then
     !$OMP CRITICAL
     write(warnline,'(a20,3g16.7)') 'in x,y,z=',x,y,z
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'w_lorentz =',w_lorentz
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'vxlow =',vxlow
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'vylow =',vylow
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'vzlow =',vzlow
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'velx =',velx
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'vely =',vely
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'velz =',velz
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'Bxlow =',Bxlow
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'Bylow =',Bylow
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'Bzlow =',Bzlow
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'Bx =',Bx
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'By =',By
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'Bz =',Bz
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'det(g_{ij}) =',det
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'uxx =',uxx
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'uxy =',uxy
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'uxz =',uxz
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'uyy =',uyy
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'uyz =',uyz
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'uzz =',uzz
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'udens =',udens
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'Sx =',usx
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'Sy =',usy
     call CCTK_WARN(2,warnline)
     write(warnline,'(a20,1g16.7)') 'Sz =',usz
     call CCTK_WARN(2,warnline)

!!$     do w_lorentz=1.0d0,100.0d0,1.0d0
!!$        call jacobian_mhd_polytype(handle,udens,s2,B2,BS,w_lorentz,f,df)
!!$        write(warnline,'(a20,3g16.7)') 'W, f, df =', w_lorentz,f,df
!!$        call CCTK_WARN(1,warnline)
!!$      end do

     call CCTK_WARN(2,"Velocity greater than the speed of light!")
     call CCTK_WARN(2,"(C2P poly) But this point is inside an AH or so close to atmo that I ignored it and set the point to atmo")
     !$OMP END CRITICAL
     ! check that it is really inside the atmosphere!!!

     if(present(mhd_c2p_mask_local)) then
        !set this mask to 1 to identify the failure
        mhd_c2p_mask_local=1
      end if

     !atmosphere
        rho = Spritz_rho_min
        udens = rho
        dens = sqrt(det)*rho
        w_lorentz = 1.0d0
        sx = 0.d0
        sy = 0.d0
        sz = 0.d0
        s2 = 0.d0
        usx = 0.d0
        usy = 0.d0
        usz = 0.d0
        velx = 0.0d0
        vely = 0.0d0
        velz = 0.0d0
        BS = 0.0d0
        B2 = (Bx*Bxlow + By*Bylow + Bz*Bzlow)!this is true only if v=0
        !Since I don't know what to do in this case for eps and rho
        !I set them using a polytropic EOS...
        !(I don't know if this can create problems...
        !...but look here if you search for a bug in con2prim)
        !press= EOS_Pressure(handle, rho, eps)
        !eps = EOS_SpecificIntEnergy(handle,rho,press)

        !!!!!!
        xrho = rho
        xeps = eps
        call EOS_Omni_press(handle,keytemp,Eos_rf_prec,n,&
                       xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
        press = xpress(1)
        call EOS_Omni_EpsFromPress(handle,keytemp,Eos_rf_prec,n,&
                              xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
        eps = xeps(1)
        !!!!!!

        !idem here
        tau  = sqrt(det) * ((rho+rho*eps+press+B2)*w_lorentz*w_lorentz - (press+0.5d0*B2)) - dens
        return
  end if

  if(eps .lt. 0.0d0) then
     if (Spritz_reflevel.ge.Spritz_c2p_warn_from_reflevel) then
        !$OMP CRITICAL
        call CCTK_WARN(1,'Con2prim_ptPolytype: stopping the code.')
        call CCTK_WARN(1, '   specific internal energy just went below 0! ')
        write(warnline,'(a28,i2)') 'on carpet reflevel: ',Spritz_reflevel
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,3g16.7)') 'xyz location: ',x,y,z
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g16.7)') 'radius: ',r
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,3g16.7)') 'velocities: ',velx,vely,velz
        call CCTK_WARN(1,warnline)
        call CCTK_WARN(Spritz_c2p_warnlevel, " - CASE B - Specific internal energy negative")
        !$OMP END CRITICAL
     else
        !$OMP CRITICAL
        write(warnline,'(a60,i2)') 'Con2Prim: eps negative, but I was told to ignore level: ',Spritz_reflevel
        call CCTK_WARN(1,warnline)
        call CCTK_WARN(2,'For security reason I reset the point to be atmosphere to avoid NaNs')
        !$OMP END CRITICAL
        !atmosphere
        rho = Spritz_rho_min
        udens = rho
        dens = sqrt(det)*rho
        w_lorentz = 1.0d0
        sx = 0.d0
        sy = 0.d0
        sz = 0.d0
        s2 = 0.d0
        usx = 0.d0
        usy = 0.d0
        usz = 0.d0
        velx = 0.0d0
        vely = 0.0d0
        velz = 0.0d0
        BS=0.0d0
        B2 = (Bx*Bxlow + By*Bylow + Bz*Bzlow)!this is true only if v=0
        !Since I don't know what to do in this case for eps and rho
        !I set them using a polytropic EOS...
        !(I don't know if this can create problems...
        !...but look here if you search for a bug in con2prim)
        !press= EOS_Pressure(handle, rho, eps)
        !eps = EOS_SpecificIntEnergy(handle,rho,press)

        !!!!!!
        xrho = rho
        xeps = eps
        call EOS_Omni_press(handle,keytemp,Eos_rf_prec,n,&
                       xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
        press = xpress(1)
        call EOS_Omni_EpsFromPress(handle,keytemp,Eos_rf_prec,n,&
                              xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
        eps = xeps(1)
        !!!!!!

        !idem here
        tau  = sqrt(det) * ((rho+rho*eps+press+B2)*w_lorentz*w_lorentz - (press+0.5d0*B2)) - dens
     endif
  end if


!write(*,*) 'Con2Prim_ptPolytype_mhd OK'

end subroutine Con2Prim_ptPolytype_mhd





 /*@@
   @routine    jacobian_mhd
   @date
   @author     Bruno Giacomazzo
   @desc
   Compute the Jacobian and the values of the equation used to get the values
   of rho, p and w_lorentz from the conserved variables.
   The Jacobian is computed analytically
   @enddesc
   @calls
   @calledby
   @history
   @endhistory

@@*/

subroutine jacobian_mhd(handle,udens,s2,utau,B2,BS,unk,fvec,fjac)
  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

!#ifdef _EOS_BASE_INC_
!#undef _EOS_BASE_INC_
!#endif
!#include "EOS_Base.inc"

  CCTK_INT :: handle
  CCTK_REAL :: udens,s2,utau,B2,BS
  CCTK_REAL,dimension(2) :: unk,fvec
  CCTK_REAL,dimension(2,2) :: fjac

  CCTK_REAL :: Z,rho,p,w_lorentz,w,eps,depsdP,depsdrho,dZdP,dZdW

  CCTK_REAL independent(2), dependent(3), epsdependent(3)
  CCTK_REAL :: eos_gamma
  CCTK_INT :: ierr

  !!!!!!
  ! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1),one(1)
  CCTK_REAL :: xdpdeps(1),xdpdrho(1)

  n = 1;keytemp = 0;anyerr = 0;keyerr = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
  one = 1.0d0
  ! end EOS Omni vars
  !!!!!

  rho =  udens/unk(2)

  p = unk(1)
  w_lorentz = unk(2)

  if (use_PPEOS .eq. 0) then
     !this is only valid for Ideal Fluid p=rho*eps*(gamma-1)
     !trick to obtain the value of eos_ideal_fluid_gamma from thorn EOS_Ideal
     !eos_gamma = EOS_Pressure(handle,1.0d0,1.0d0) +1.0d0

     !!!!!!
     call EOS_Omni_press(handle,keytemp,Eos_rf_prec,n,&
                    one,one,xtemp,xye,xpress,keyerr,anyerr)
     eos_gamma = xpress(1) + 1.0d0
     !!!!!!

     w = rho + eos_gamma*p/(eos_gamma-1.0d0)
  else
     !eps = EOS_SpecificIntEnergy(handle,rho,p)

     !!!!!!
     xrho = rho
     xpress = p
     call EOS_Omni_EpsFromPress(handle,keytemp,Eos_rf_prec,n,&
                           xrho,xeps,xtemp,xye,xpress,xeps,keyerr,anyerr)
     eps = xeps(1)
     !!!!!!

     w = rho+rho*eps+p
  end if


  Z = w*w_lorentz**2

  fvec(1) = s2 - ((Z+B2)**2*(w_lorentz**2-1.0d0)/w_lorentz**2 - (2.0d0*Z+B2)*BS**2/Z**2)
  fvec(2) = utau - (Z + B2 - p - B2/(2.0d0*w_lorentz**2) - BS**2/(2.0d0*Z**2) - udens)

  !Note that rho is only a function of w_lorentz becuase rho=D/w_lorentz
  !while eps depends both on P and rho=rho(w_lorentz)
  !NOTE: Z=(rho+rho*eps+P)*w_lorentz**2 = (D/w_lorentz + eps*D/w_lorentz + p)*w_lorentz**2
  if (use_PPEOS .eq. 0) then
     !this is only valid for Ideal Fluid p=rho*eps*(gamma-1)
     dZdp = w_lorentz**2*eos_gamma/(eos_gamma-1.0d0)
     dZdW = udens + 2.0d0*w_lorentz*p*eos_gamma/(eos_gamma-1.0d0)
  else
     !depsdP = 1.0d0/EOS_DPressByDEps(handle,rho,eps)
     !depsdrho = EOS_DPressByDRho(handle,rho,eps)/EOS_DPressByDEps(handle,rho,eps)

     !!!!!!
     xrho = rho
     xeps = eps
     call EOS_Omni_DPressByDEps(handle,keytemp,Eos_rf_prec,n,&
                           xrho,xeps,xtemp,xye,xdpdeps,keyerr,anyerr)
     call EOS_Omni_DPressByDRho(handle,keytemp,Eos_rf_prec,n,&
                           xrho,xeps,xtemp,xye,xdpdrho,keyerr,anyerr)
     depsdP = 1.0d0/xdpdeps(1)
     depsdrho = xdpdrho(1)/xdpdeps(1)
     !!!!!!

     dZdp = w_lorentz**2*(rho*depsdP+1)
     dZdW = udens*(1+eps-udens/w_lorentz*depsdrho)+2.0*p*w_lorentz
  end if

  fjac(1,1) = (2.0d0*B2*(-1.0d0+1.0d0/w_lorentz**2)-2.0d0*B2*BS**2/Z**3 - &
       2.0d0*BS**2/Z**2+(-2.0d0+2.0d0/w_lorentz**2)*Z)*dZdp

  fjac(1,2) = 1.0d0/w_lorentz**3*(-2.0d0*(B2+Z)**2-1.0d0/Z**3*&
       (2.0d0*w_lorentz*(B2*BS**2*w_lorentz**2 + BS**2*w_lorentz**2*Z + &
       B2*(-1.0d0+w_lorentz**2)*Z**3+(w_lorentz**2-1)*Z**4)*dZdW))

  fjac(2,1) = 1.0d0+(-1.0d0-BS**2/Z**3)*dZdp
  fjac(2,2) = -B2/w_lorentz**3 + (-1.0d0-BS**2/Z**3)*dZdW
end subroutine jacobian_mhd

 /*@@
   @routine    jacobian_mhd_1eq
   @date
   @author     Bruno Giacomazzo
   @desc
   Compute the Jacobian for the Neilsen et al 2006 method
   The Jacobian is computed analytically
   @enddesc
   @calls
   @calledby
   @history
   @endhistory

@@*/
subroutine jacobian_mhd_1eq(handle,udens,s2,utau,B2,BS,unk,fvec,fjac)
  implicit none

  DECLARE_CCTK_PARAMETERS

!#ifdef _EOS_BASE_INC_
!#undef _EOS_BASE_INC_
!#endif
!#include "EOS_Base.inc"

  CCTK_INT :: handle
  CCTK_REAL :: udens,s2,utau,B2,BS,unk,fvec,fjac
  CCTK_REAL :: rho,p,w_lorentz
  CCTK_REAL :: eos_gamma

  !!!!!!
  ! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),one(1)

  n = 1;keytemp = 0;anyerr = 0;keyerr = 0
  xpress = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
  one = 1.0d0
  ! end EOS Omni vars
  !!!!!

  !trick to obtain the value of eos_ideal_fluid_gamma from thorn EOS_Ideal
  !eos_gamma = EOS_Pressure(handle,1.0d0,1.0d0) + 1.0d0

  !!!!!!
  call EOS_Omni_press(handle,keytemp,Eos_rf_prec,n,&
                 one,one,xtemp,xye,xpress,keyerr,anyerr)
  eos_gamma = xpress(1) + 1.0d0
  !!!!!!

  w_lorentz = sqrt((unk**2*(unk+B2)**2) / (unk**2*(unk+B2)**2-unk**2*s2-(2.0d0*unk+B2)*BS**2))

  if ((unk**2*(unk+B2)**2-unk**2*s2-(2.0d0*unk+B2)*BS**2).le.0.0d0) then
     write(*,*) 'unk**2*(unk+B2)**2-unk**2*s2-(2.0d0*unk+B2)*BS**2',unk**2*(unk+B2)**2-unk**2*s2-(2.0d0*unk+B2)*BS**2
     write(*,*) 'unk=',unk
     write(*,*) 'B2=',B2
     write(*,*) 's2=',s2
     write(*,*) 'BS=',BS
     !$OMP CRITICAL
     call CCTK_WARN(0,"w_lorentz is NaN in jacobian_mhd_1eq")
     !$OMP END CRITICAL
  end if

  rho =  udens/w_lorentz

  p = (unk/w_lorentz**2-udens/w_lorentz)*(eos_gamma-1.0d0)/eos_gamma


  fvec = unk - p -0.5d0*((BS**2)/(unk**2)+B2/w_lorentz**2) + B2 - utau - udens

  fjac = 1 + BS**2/unk**3 + (B2*(-(S2*unk**2) + unk**2*(B2 + unk)**2 - BS**2*(B2 + 2*unk)))/(unk**2*(B2 + unk)**3) + &
       (B2*(-(S2*unk**2) + unk**2*(B2 + unk)**2 - BS**2*(B2 + 2*unk)))/(unk**3*(B2 + unk)**2) - &
       (B2*(-BS**2 + unk*(-S2 + (B2 + unk)*(B2 + 2*unk))))/(unk**2*(B2 + unk)**2) + &
       ((-1 + eos_gamma)*(2*unk**2*(-(S2*unk**2) + unk**2*(B2 + unk)**2 - BS**2*(B2 + 2*unk)) + &
       unk*(B2 + unk)*(-(S2*unk**2) + unk**2*(B2 + unk)**2 - BS**2*(B2 + 2*unk)) + &
       udens*Sqrt((unk**2*(B2 + unk)**2)/(-(B2*BS**2) - 2*BS**2*unk + (B2**2 - S2)*unk**2 + &
       2*B2*unk**3 + unk**4))*(B2**2*BS**2 + 3*B2*BS**2*unk + unk**2*(3*BS**2 + S2*unk)) - &
       unk**2*(B2 + unk)*(-2*BS**2 + 2*unk*(-S2 + (B2 + unk)*(B2 + 2*unk)))))/(eos_gamma*unk**3*(B2 + unk)**3)

end subroutine jacobian_mhd_1eq


 /*@@
   @routine    jacobian_mhd_polytype
   @date
   @author     Bruno Giacomazzo
   @desc
   Compute the Jacobian and the values of the equation used to get the values
   of rho, p and w_lorentz from the conserved variables.
   The Jacobian is computed analytically.
   A polytropyc equation of state is used.
   @enddesc
   @calls
   @calledby
   @history
   @endhistory

@@*/

subroutine jacobian_mhd_polytype(handle,udens,s2,B2,BS,unk,fvec,fjac)
  implicit none

  DECLARE_CCTK_PARAMETERS

!#ifdef _EOS_BASE_INC_
!#undef _EOS_BASE_INC_
!#endif
!#include "EOS_Base.inc"

  CCTK_INT :: handle
  CCTK_REAL :: udens,s2,B2,BS
  CCTK_REAL :: unk,fvec
  CCTK_REAL :: fjac

  CCTK_REAL :: Z,rho,p,w_lorentz,w,eps,dZdW
  CCTK_REAL :: eos_gamma, eos_kappa

  !!!!!!
  ! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1),one(1),zero(1)

  n = 1;keytemp = 0;anyerr = 0;keyerr = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = 0.0d0;xye = 0.0d0
  zero = 0.0d0;one = 1.0d0
  ! end EOS Omni vars
  !!!!!

  rho = udens/unk
  w_lorentz = unk

  !P = EOS_Pressure(handle,rho,0.0d0)
  !eps = EOS_SpecificIntEnergy(handle,rho,0.0d0)

  !!!!!!
  xrho = rho
  call EOS_Omni_press(handle,keytemp,Eos_rf_prec,n,&
                 xrho,zero,xtemp,xye,xpress,keyerr,anyerr)
  P = xpress(1)
  call EOS_Omni_EpsFromPress(handle,keytemp,Eos_rf_prec,n,&
                        xrho,xeps,xtemp,xye,zero,xeps,keyerr,anyerr)
  eps = xeps(1)
  !!!!!!

  w = rho*(1.0+eps)+P

  Z = w*w_lorentz**2

  fvec = s2 - ((Z+B2)**2*(w_lorentz**2-1)/w_lorentz**2 - (2.0*Z+B2)*BS**2/Z**2)

  !and now the Jacobian matrix computed analytically
  !eos_kappa = EOS_Pressure(handle,1.0d0,0.0d0)
  !eos_gamma = eos_kappa/EOS_SpecificIntEnergy(handle,1.0d0,0.0d0) + 1.0d0

  !!!!!!
  call EOS_Omni_press(handle,keytemp,Eos_rf_prec,n,&
                 one,zero,xtemp,xye,xpress,keyerr,anyerr)
  eos_kappa = xpress(1)
  call EOS_Omni_EpsFromPress(handle,keytemp,Eos_rf_prec,n,&
                        one,xeps,xtemp,xye,zero,xeps,keyerr,anyerr)
  eos_gamma = eos_kappa/xeps(1) + 1.0d0
  !!!!!!

!!$  dZdW = 2.0d0*(rho + rho**(eos_gamma-1)*eos_gamma*eos_kappa/(eos_gamma-1.0d0))*w_lorentz !!WRONG
  dZdW =  ( rho + rho**(eos_gamma)*eos_gamma*(2.0d0-eos_gamma)*eos_kappa/(eos_gamma-1.0d0) )*w_lorentz

!!$  fjac = -B2/w_lorentz**3+(-1.0d0-BS**2/Z**3)*dZdW !!WRONG

  fjac = 1.0d0/w_lorentz**3*(-2.0d0*(B2+Z)**2-1.0d0/Z**3*&
       (2.0d0*w_lorentz*(B2*BS**2*w_lorentz**2 + BS**2*w_lorentz**2*Z + &
       B2*(-1.0d0+w_lorentz**2)*Z**3+(w_lorentz**2-1)*Z**4)*dZdW))

end subroutine jacobian_mhd_polytype

subroutine SetToAtmo(rho,dens,det,w_lorentz,press,&
                        sx,sy,sz,s2,usx,usy,usz,tau,&
                        velx,vely,velz,eps,temp,ye,handle,&
                        BS,B2,Bx,By,Bz,Bxlow,Bylow,Bzlow) !sets point to atmosphere
  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, intent(inout) :: ye
  CCTK_REAL, intent(out)   :: rho,dens,w_lorentz,press,sx,sy,sz,s2,usx,usy,usz,tau,velx,vely,velz,BS,B2,eps,temp
  CCTK_REAL, intent(in)    :: det,Bx,By,Bz,Bxlow,Bylow,Bzlow
  CCTK_INT,  intent(in)    :: handle

  !!!!!!
  ! begin EOS Omni vars
  CCTK_INT  :: n,keytemp,anyerr,keyerr(1)
  CCTK_REAL :: xpress(1),xtemp(1),xye(1),xeps(1),xrho(1)

  n = 1;keytemp = 0;anyerr = 0;keyerr = 0
  xpress = 0.0d0;xeps = 0.0d0;xtemp = temp; xye = ye
  ! end EOS Omni vars
  !!!!!!

  rho  = rho_abs_min
  temp = Spritz_hot_atmo_temp
  dens = sqrt(det)*rho
  w_lorentz = 1.0d0
  sx = 0.d0
  sy = 0.d0
  sz = 0.d0
  s2 = 0.d0
  usx = 0.d0
  usy = 0.d0
  usz = 0.d0
  velx = 0.0d0
  vely = 0.0d0
  velz = 0.0d0
  BS=0.0d0
  B2 = (Bx*Bxlow + By*Bylow + Bz*Bzlow)

  xrho = rho
  xeps = eps !will be overwritten
  xtemp= temp
  xye = ye
  keytemp = 1 ! don't change temperature
  call EOS_Omni_press(handle,keytemp,Eos_rf_prec,n,&
                 xrho,xeps,xtemp,xye,xpress,keyerr,anyerr)
  press = xpress(1)
  eps = xeps(1)
  keytemp = 0
  tau  = sqrt(det) * ((rho+rho*eps+press+B2)*w_lorentz*w_lorentz - (press+0.5d0*B2)) - dens
end subroutine SetToAtmo


subroutine sync_Spritzmhd_c2p_mask(CCTK_ARGUMENTS)
! a dummy routine to syncronise Spritzmhd_c2p_mask
  implicit none

  DECLARE_CCTK_ARGUMENTS

  return

end subroutine sync_Spritzmhd_c2p_mask
