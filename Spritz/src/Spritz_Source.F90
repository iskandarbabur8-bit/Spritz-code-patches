 /*@@
   @file      Spritz_Source.F90
   @date      Sat Jan 26 02:03:56 2002
   @author    Ian Hawke
   @desc
   The geometric source terms for the matter evolution
   @enddesc
 @@*/

! Second order f.d.

#define DIFF_X_2(q) 0.5d0 * (q(i+1,j,k) - q(i-1,j,k)) * idx
#define DIFF_Y_2(q) 0.5d0 * (q(i,j+1,k) - q(i,j-1,k)) * idy
#define DIFF_Z_2(q) 0.5d0 * (q(i,j,k+1) - q(i,j,k-1)) * idz

! Fourth order f.d.

#define DIFF_X_4(q) (-q(i+2,j,k) + 8.d0 * q(i+1,j,k) - 8.d0 * q(i-1,j,k) + \
                      q(i-2,j,k)) / 12.d0 * idx
#define DIFF_Y_4(q) (-q(i,j+2,k) + 8.d0 * q(i,j+1,k) - 8.d0 * q(i,j-1,k) + \
                      q(i,j-2,k)) / 12.d0 * idy
#define DIFF_Z_4(q) (-q(i,j,k+2) + 8.d0 * q(i,j,k+1) - 8.d0 * q(i,j,k-1) + \
                      q(i,j,k-2)) / 12.d0 * idz

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#define velx(i,j,k) vel(i,j,k,1)
#define vely(i,j,k) vel(i,j,k,2)
#define velz(i,j,k) vel(i,j,k,3)

 /*@@
   @routine    SourceTerms
   @date       Sat Jan 26 02:04:21 2002
   @author     Ian Hawke
   @desc
   Calculate the geometric source terms and add to the update GFs
   @enddesc
   @calls
   @calledby
   @history
   Minor alterations of routine from GR3D.
   @endhistory

@@*/

subroutine SourceTerms(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: i, j, k, nx, ny, nz
  CCTK_REAL :: one, two, half
  CCTK_REAL :: t00, t0x, t0y, t0z, txx, txy, txz, tyy, tyz, tzz
  CCTK_REAL :: sqrtdet, det, uxx, uxy, uxz, uyy, uyz, uzz, rhoenthalpyW2
  CCTK_REAL :: shiftx, shifty, shiftz, velxshift, velyshift, velzshift
  CCTK_REAL :: vlowx, vlowy, vlowz
  CCTK_REAL :: dx_betax, dx_betay, dx_betaz, dy_betax, dy_betay,&
       dy_betaz, dz_betax, dz_betay, dz_betaz
  CCTK_REAL :: dx_alp, dy_alp, dz_alp
  CCTK_REAL :: tau_source, sx_source, sy_source, sz_source
  CCTK_REAL :: localgxx,localgxy,localgxz,localgyy,localgyz,localgzz
  CCTK_REAL :: dx_gxx, dx_gxy, dx_gxz, dx_gyy, dx_gyz, dx_gzz
  CCTK_REAL :: dy_gxx, dy_gxy, dy_gxz, dy_gyy, dy_gyz, dy_gzz
  CCTK_REAL :: dz_gxx, dz_gxy, dz_gxz, dz_gyy, dz_gyz, dz_gzz
  CCTK_REAL :: dx, dy, dz, idx, idy, idz
  CCTK_REAL :: psi4pt, dx_psi4, dy_psi4, dz_psi4
  CCTK_REAL :: shiftshiftk, shiftkx, shiftky, shiftkz
  CCTK_REAL :: sumTK
  CCTK_REAL :: halfshiftdgx, halfshiftdgy, halfshiftdgz
  CCTK_REAL :: halfTdgx, halfTdgy, halfTdgz
  CCTK_REAL :: Bxlow, Bylow, Bzlow, bst, bsx, bsy, bsz, bs2
  CCTK_REAL :: bsxlow, bsylow, bszlow, t0lowx, t0lowy, t0lowz, beta2
  CCTK_REAL :: betax_low, betay_low, betaz_low
  CCTK_REAL :: divB2tmp

  CCTK_REAL :: Keppens_cd_local, divBipo, divBimo, divBjpo, divBjmo, divBkpo, divBkmo, dt

  logical, allocatable, dimension (:,:,:) :: force_spatial_second_order

  character(len=300) warnline

  one = 1.0d0
  two = 2.0d0
  half = 0.5d0
  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  dx = CCTK_DELTA_SPACE(1)
  dy = CCTK_DELTA_SPACE(2)
  dz = CCTK_DELTA_SPACE(3)
  idx = 1.d0/dx
  idy = 1.d0/dy
  idz = 1.d0/dz

  dt = CCTK_DELTA_TIME

!!$  write(*,*) 'dt=',dt, ' dx=',dx

!!$  Initialize the update terms to be zero.
!!$  This will guarantee that no garbage in the boundaries is updated.

  densrhs = 0.d0
  srhs = 0.d0
  taurhs = 0.d0

  if (evolve_Y_e .ne. 0) then
    y_e_con_rhs = 0.0d0
  endif

  Bxrhs = 0.0d0
  Byrhs = 0.0d0
  Bzrhs = 0.0d0
  Avec_x_rhs=0.0d0
  Avec_y_rhs=0.0d0
  Avec_z_rhs=0.0d0

  if((Spritz_mhd_fluxCD_Lorenz.ne.0)) then
     psi_mhdrhs=0.0d0
  end if

!!$  Set up the array for checking the order. We switch to second order
!!$  differencing at boundaries.
!!$  Copied straight from BSSN.

  allocate (force_spatial_second_order(nx,ny,nz))
  force_spatial_second_order = .FALSE.

  if (spatial_order > 2) then
    do k = 1 + Spritz_stencil, nz - Spritz_stencil
      do j = 1 + Spritz_stencil, ny - Spritz_stencil
        do i = 1 + Spritz_stencil, nx - Spritz_stencil
          if ((i < 3).or.(i > cctk_lsh(1) - 2).or. &
               (j < 3).or.(j > cctk_lsh(2) - 2).or. &
               (k < 3).or.(k > cctk_lsh(3) - 2) ) then
            force_spatial_second_order(i,j,k) = .TRUE.
           end if
        end do
      end do
    end do
  end if

  !$OMP PARALLEL DO PRIVATE(i, j, k, &
  !$OMP t00, t0x, t0y, t0z, txx, txy, txz, tyy, tyz, tzz, &
  !$OMP sqrtdet, det, uxx, uxy, uxz, uyy, uyz, uzz, rhoenthalpyW2, &
  !$OMP shiftx, shifty, shiftz, velxshift, velyshift, velzshift, &
  !$OMP vlowx, vlowy, vlowz, &
  !$OMP dx_betax, dx_betay, dx_betaz, dy_betax, dy_betay,&
  !$OMP dy_betaz, dz_betax, dz_betay, dz_betaz, &
  !$OMP dx_alp, dy_alp, dz_alp, &
  !$OMP tau_source, sx_source, sy_source, sz_source, &
  !$OMP localgxx,localgxy,localgxz,localgyy,localgyz,localgzz, &
  !$OMP dx_gxx, dx_gxy, dx_gxz, dx_gyy, dx_gyz, dx_gzz, &
  !$OMP dy_gxx, dy_gxy, dy_gxz, dy_gyy, dy_gyz, dy_gzz, &
  !$OMP dz_gxx, dz_gxy, dz_gxz, dz_gyy, dz_gyz, dz_gzz, &
  !$OMP psi4pt, dx_psi4, dy_psi4, dz_psi4, &
  !$OMP shiftshiftk, shiftkx, shiftky, shiftkz, &
  !$OMP sumTK, &
  !$OMP halfshiftdgx, halfshiftdgy, halfshiftdgz, &
  !$OMP halfTdgx, halfTdgy, halfTdgz, &
  !$OMP Bxlow, Bylow, Bzlow, bst, bsx, bsy, bsz, bs2, &
  !$OMP bsxlow, bsylow, bszlow, t0lowx, t0lowy, t0lowz, beta2, &
  !$OMP betax_low, betay_low, betaz_low, &
  !$OMP divB2tmp, warnline, &
  !$OMP divBipo, divBimo, divBjpo, divBjmo, divBkpo, divBkmo, Keppens_cd_local)
  do k=1 + Spritz_stencil,nz - Spritz_stencil
    do j=1 + Spritz_stencil,ny - Spritz_stencil
      do i=1 + Spritz_stencil,nx - Spritz_stencil


        local_spatial_order = spatial_order
        if (force_spatial_second_order(i,j,k)) then
          local_spatial_order = 2
        end if

!!$        Set the metric terms.

        if (conformal_state .eq. 0) then

          localgxx = gxx(i,j,k)
          localgxy = gxy(i,j,k)
          localgxz = gxz(i,j,k)
          localgyy = gyy(i,j,k)
          localgyz = gyz(i,j,k)
          localgzz = gzz(i,j,k)

        else

          psi4pt = psi(i,j,k)**4

          localgxx = psi4pt*gxx(i,j,k)
          localgxy = psi4pt*gxy(i,j,k)
          localgxz = psi4pt*gxz(i,j,k)
          localgyy = psi4pt*gyy(i,j,k)
          localgyz = psi4pt*gyz(i,j,k)
          localgzz = psi4pt*gzz(i,j,k)

        end if


        call SpatialDeterminant(localgxx, localgxy, localgxz,&
             localgyy, localgyz, localgzz, det)
        if (det<0.0d0) then
           !$OMP CRITICAL
           write(warnline,'(a50,g20.7)') 'sqrtdet is NaN in Spritz Source, det: ',det
           call CCTK_WARN(1,warnline)
           write(warnline,'(a20,g20.7)') 'gxx=',localgxx
           call CCTK_WARN(1,warnline)
           write(warnline,'(a20,g20.7)') 'gxy=',localgxy
           call CCTK_WARN(1,warnline)
           write(warnline,'(a20,g20.7)') 'gxz=',localgxz
           call CCTK_WARN(1,warnline)
           write(warnline,'(a20,g20.7)') 'gyy=',localgyy
           call CCTK_WARN(1,warnline)
           write(warnline,'(a20,g20.7)') 'gyz=',localgyz
           call CCTK_WARN(1,warnline)
           write(warnline,'(a20,g20.7)') 'gzz=',localgzz
           call CCTK_WARN(1,warnline)
           write(warnline,'(a20,g20.7)') 'dens: ',dens(i,j,k)
           call CCTK_WARN(1,warnline)
           write(warnline,'(a20,g20.7)') ' tau: ',tau(i,j,k)
           call CCTK_WARN(1,warnline)
           write(warnline,'(a20,g20.7)') '  sx: ',scon(i,j,k,1)
           call CCTK_WARN(1,warnline)
           write(warnline,'(a20,g20.7)') '  sy: ',scon(i,j,k,2)
           call CCTK_WARN(1,warnline)
           write(warnline,'(a20,g20.7)') '  sz: ',scon(i,j,k,3)
           call CCTK_WARN(1,warnline)
           write(warnline,'(a20,g20.7)') ' dBx: ',dBx(i,j,k)
           call CCTK_WARN(1,warnline)
           write(warnline,'(a20,g20.7)') ' dBy: ',dBy(i,j,k)
           call CCTK_WARN(1,warnline)
           write(warnline,'(a20,g20.7)') ' dBz: ',dBz(i,j,k)
           call CCTK_WARN(1,warnline)
           write(warnline,'(a20,g20.7)') 'x=',x(i,j,k)
           call CCTK_WARN(1,warnline)
           write(warnline,'(a20,g20.7)') 'y=',y(i,j,k)
           call CCTK_WARN(1,warnline)
           write(warnline,'(a20,g20.7)') 'z=',z(i,j,k)
           call CCTK_WARN(1,warnline)
           write(warnline,'(a20,g20.7)') 'dx=',dx
           call CCTK_WARN(1,warnline)
           write(warnline,'(a20,g20.7)') 'dy=',dy
           call CCTK_WARN(1,warnline)
           write(warnline,'(a20,g20.7)') 'dz=',dz
           call CCTK_WARN(0,warnline)
           !$OMP END CRITICAL
        end if
        sqrtdet = sqrt(det)
        call UpperMetric(uxx, uxy, uxz, uyy, uyz, uzz, det, localgxx,&
             localgxy, localgxz, localgyy, localgyz, localgzz)

!!$        All the matter variables (except velocity) always appear
!!$        together in this form

        rhoenthalpyW2 = (rho(i,j,k)*(one + eps(i,j,k)) + press(i,j,k))*&
             w_lorentz(i,j,k)**2

        if (shift_state .ne. 0) then

          shiftx = betax(i,j,k)
          shifty = betay(i,j,k)
          shiftz = betaz(i,j,k)

          if (local_spatial_order .eq. 2) then

            dx_betax = DIFF_X_2(betax)
            dx_betay = DIFF_X_2(betay)
            dx_betaz = DIFF_X_2(betaz)

            dy_betax = DIFF_Y_2(betax)
            dy_betay = DIFF_Y_2(betay)
            dy_betaz = DIFF_Y_2(betaz)

            dz_betax = DIFF_Z_2(betax)
            dz_betay = DIFF_Z_2(betay)
            dz_betaz = DIFF_Z_2(betaz)

          else

            dx_betax = DIFF_X_4(betax)
            dx_betay = DIFF_X_4(betay)
            dx_betaz = DIFF_X_4(betaz)

            dy_betax = DIFF_Y_4(betax)
            dy_betay = DIFF_Y_4(betay)
            dy_betaz = DIFF_Y_4(betaz)

            dz_betax = DIFF_Z_4(betax)
            dz_betay = DIFF_Z_4(betay)
            dz_betaz = DIFF_Z_4(betaz)

          end if

        else

          shiftx = 0.0d0
          shifty = 0.0d0
          shiftz = 0.0d0

          dx_betax = 0.0d0
          dx_betay = 0.0d0
          dx_betaz = 0.0d0

          dy_betax = 0.0d0
          dy_betay = 0.0d0
          dy_betaz = 0.0d0

          dz_betax = 0.0d0
          dz_betay = 0.0d0
          dz_betaz = 0.0d0

        endif

        velxshift = velx(i,j,k) - shiftx/alp(i,j,k)
        velyshift = vely(i,j,k) - shifty/alp(i,j,k)
        velzshift = velz(i,j,k) - shiftz/alp(i,j,k)
        vlowx = velx(i,j,k)*localgxx + vely(i,j,k)*localgxy +&
             velz(i,j,k)*localgxz
        vlowy = velx(i,j,k)*localgxy + vely(i,j,k)*localgyy +&
             velz(i,j,k)*localgyz
        vlowz = velx(i,j,k)*localgxz + vely(i,j,k)*localgyz +&
             velz(i,j,k)*localgzz

        !B_j
        Bxlow = Bx(i,j,k)*localgxx + By(i,j,k)*localgxy +&
             Bz(i,j,k)*localgxz
        Bylow = Bx(i,j,k)*localgxy + By(i,j,k)*localgyy +&
             Bz(i,j,k)*localgyz
        Bzlow = Bx(i,j,k)*localgxz + By(i,j,k)*localgyz +&
             Bz(i,j,k)*localgzz

        !these are b^mu
        bst=w_lorentz(i,j,k)*(Bx(i,j,k)*vlowx+ &
             By(i,j,k)*vlowy+Bz(i,j,k)*vlowz)/alp(i,j,k)
        bsx=(Bx(i,j,k)+alp(i,j,k)*bst*w_lorentz(i,j,k)*velxshift)/w_lorentz(i,j,k)
        bsy=(By(i,j,k)+alp(i,j,k)*bst*w_lorentz(i,j,k)*velyshift)/w_lorentz(i,j,k)
        bsz=(Bz(i,j,k)+alp(i,j,k)*bst*w_lorentz(i,j,k)*velzshift)/w_lorentz(i,j,k)
        !b^mu b_mu
        bs2=(Bx(i,j,k)*Bxlow+By(i,j,k)*Bylow+Bz(i,j,k)*Bzlow + &
             alp(i,j,k)**2*bst**2)/w_lorentz(i,j,k)**2

!!$        For a change, these are T^{ij}
        !g^{00} = -1/alp^2
        !g^{0j} = beta^j/alp^2
        !g^{ij} = \gamma^{ij} - beta^i beta^j/alp^2

        t00 = (rhoenthalpyW2 - press(i,j,k))/(alp(i,j,k)**2)
        t0x = rhoenthalpyW2*velxshift/alp(i,j,k) +&
             press(i,j,k)*shiftx/(alp(i,j,k)**2)
        t0y = rhoenthalpyW2*velyshift/alp(i,j,k) +&
             press(i,j,k)*shifty/(alp(i,j,k)**2)
        t0z = rhoenthalpyW2*velzshift/alp(i,j,k) +&
             press(i,j,k)*shiftz/(alp(i,j,k)**2)
        txx = rhoenthalpyW2*velxshift*velxshift +&
             press(i,j,k)*(uxx - shiftx*shiftx/(alp(i,j,k)**2))
        txy = rhoenthalpyW2*velxshift*velyshift +&
             press(i,j,k)*(uxy - shiftx*shifty/(alp(i,j,k)**2))
        txz = rhoenthalpyW2*velxshift*velzshift +&
             press(i,j,k)*(uxz - shiftx*shiftz/(alp(i,j,k)**2))
        tyy = rhoenthalpyW2*velyshift*velyshift +&
             press(i,j,k)*(uyy - shifty*shifty/(alp(i,j,k)**2))
        tyz = rhoenthalpyW2*velyshift*velzshift +&
             press(i,j,k)*(uyz - shifty*shiftz/(alp(i,j,k)**2))
        tzz = rhoenthalpyW2*velzshift*velzshift +&
             press(i,j,k)*(uzz - shiftz*shiftz/(alp(i,j,k)**2))

        beta2 = localgxx*shiftx**2 + localgyy*shifty**2 + localgzz*shiftz**2 + &
             2.0*(localgxy*shiftx*shifty + localgxz*shiftx*shiftz + &
             localgyz*shifty*shiftz)

        t00 = t00 + Spritz_mhd_coeff*((w_lorentz(i,j,k)**2/alp(i,j,k)**2-&
             half/alp(i,j,k)**2)*bs2-bst*bst)
        t0x = t0x + Spritz_mhd_coeff*((w_lorentz(i,j,k)**2*velxshift/alp(i,j,k)+&
             half*shiftx/alp(i,j,k)**2)*bs2-bst*bsx)
        t0y = t0y + Spritz_mhd_coeff*((w_lorentz(i,j,k)**2*velyshift/alp(i,j,k)+&
             half*shifty/alp(i,j,k)**2)*bs2-bst*bsy)
        t0z = t0z + Spritz_mhd_coeff*((w_lorentz(i,j,k)**2*velzshift/alp(i,j,k)+&
             half*shiftz/alp(i,j,k)**2)*bs2-bst*bsz)
        txx = txx + Spritz_mhd_coeff*((w_lorentz(i,j,k)**2*velxshift**2+&
             half*(uxx-shiftx**2/alp(i,j,k)**2))*bs2-bsx**2)
        txy = txy + Spritz_mhd_coeff*((w_lorentz(i,j,k)**2*velxshift*velyshift+&
             half*(uxy-shiftx*shifty/alp(i,j,k)**2))*bs2-bsx*bsy)
        txz = txz + Spritz_mhd_coeff*((w_lorentz(i,j,k)**2*velxshift*velzshift+&
             half*(uxz-shiftx*shiftz/alp(i,j,k)**2))*bs2-bsx*bsz)
        tyy = tyy + Spritz_mhd_coeff*((w_lorentz(i,j,k)**2*velyshift**2+&
             half*(uyy-shifty**2/alp(i,j,k)**2))*bs2-bsy**2)
        tyz = tyz + Spritz_mhd_coeff*((w_lorentz(i,j,k)**2*velyshift*velzshift+&
             half*(uyz-shifty*shiftz/alp(i,j,k)**2))*bs2-bsy*bsz)
        tzz = tzz + Spritz_mhd_coeff*((w_lorentz(i,j,k)**2*velzshift**2+&
             half*(uzz-shiftz**2/alp(i,j,k)**2))*bs2-bsz**2)

        betax_low = localgxx*shiftx + localgxy*shifty + localgxz*shiftz
        betay_low = localgxy*shiftx + localgyy*shifty + localgyz*shiftz
        betaz_low = localgxz*shiftx + localgyz*shifty + localgzz*shiftz

        !b_j
        bsxlow  = betax_low*bst  + localgxx*bsx  + localgxy*bsy  + localgxz*bsz
        bsylow  = betay_low*bst  + localgxy*bsx  + localgyy*bsy  + localgyz*bsz
        bszlow  = betaz_low*bst  + localgxz*bsx  + localgyz*bsy  + localgzz*bsz


        !T_i^0
        t0lowx = (rhoenthalpyW2 + Spritz_mhd_coeff*bs2*w_lorentz(i,j,k)**2)*&
             vlowx/alp(i,j,k) - &
             Spritz_mhd_coeff*bsxlow*bst
        t0lowy = (rhoenthalpyW2 + Spritz_mhd_coeff*bs2*w_lorentz(i,j,k)**2)*&
             vlowy/alp(i,j,k) - &
             Spritz_mhd_coeff*bsylow*bst
        t0lowz = (rhoenthalpyW2 + Spritz_mhd_coeff*bs2*w_lorentz(i,j,k)**2)*&
             vlowz/alp(i,j,k) - &
             Spritz_mhd_coeff*bszlow*bst

!!$        Derivatives of the lapse, metric and shift

        if (local_spatial_order .eq. 2) then

          dx_alp = DIFF_X_2(alp)
          dy_alp = DIFF_Y_2(alp)
          dz_alp = DIFF_Z_2(alp)

        else

          dx_alp = DIFF_X_4(alp)
          dy_alp = DIFF_Y_4(alp)
          dz_alp = DIFF_Z_4(alp)

        end if

        if (conformal_state .eq. 0) then

          if (local_spatial_order .eq. 2) then

            dx_gxx = DIFF_X_2(gxx)
            dx_gxy = DIFF_X_2(gxy)
            dx_gxz = DIFF_X_2(gxz)
            dx_gyy = DIFF_X_2(gyy)
            dx_gyz = DIFF_X_2(gyz)
            dx_gzz = DIFF_X_2(gzz)
            dy_gxx = DIFF_Y_2(gxx)
            dy_gxy = DIFF_Y_2(gxy)
            dy_gxz = DIFF_Y_2(gxz)
            dy_gyy = DIFF_Y_2(gyy)
            dy_gyz = DIFF_Y_2(gyz)
            dy_gzz = DIFF_Y_2(gzz)
            dz_gxx = DIFF_Z_2(gxx)
            dz_gxy = DIFF_Z_2(gxy)
            dz_gxz = DIFF_Z_2(gxz)
            dz_gyy = DIFF_Z_2(gyy)
            dz_gyz = DIFF_Z_2(gyz)
            dz_gzz = DIFF_Z_2(gzz)

          else

            dx_gxx = DIFF_X_4(gxx)
            dx_gxy = DIFF_X_4(gxy)
            dx_gxz = DIFF_X_4(gxz)
            dx_gyy = DIFF_X_4(gyy)
            dx_gyz = DIFF_X_4(gyz)
            dx_gzz = DIFF_X_4(gzz)
            dy_gxx = DIFF_Y_4(gxx)
            dy_gxy = DIFF_Y_4(gxy)
            dy_gxz = DIFF_Y_4(gxz)
            dy_gyy = DIFF_Y_4(gyy)
            dy_gyz = DIFF_Y_4(gyz)
            dy_gzz = DIFF_Y_4(gzz)
            dz_gxx = DIFF_Z_4(gxx)
            dz_gxy = DIFF_Z_4(gxy)
            dz_gxz = DIFF_Z_4(gxz)
            dz_gyy = DIFF_Z_4(gyy)
            dz_gyz = DIFF_Z_4(gyz)
            dz_gzz = DIFF_Z_4(gzz)

          end if

        else if (conformal_state > 1) then

!!$          Using conformal factor, and the derivatives are
!!$          already calculated.

          psi4pt = psi(i,j,k)**4
          dx_psi4 = 4.d0*psi4pt*psix(i,j,k)
          dy_psi4 = 4.d0*psi4pt*psiy(i,j,k)
          dz_psi4 = 4.d0*psi4pt*psiz(i,j,k)

          if (local_spatial_order .eq. 2) then

            dx_gxx = DIFF_X_2(gxx)*psi4pt + dx_psi4*gxx(i,j,k)
            dx_gxy = DIFF_X_2(gxy)*psi4pt + dx_psi4*gxy(i,j,k)
            dx_gxz = DIFF_X_2(gxz)*psi4pt + dx_psi4*gxz(i,j,k)
            dx_gyy = DIFF_X_2(gyy)*psi4pt + dx_psi4*gyy(i,j,k)
            dx_gyz = DIFF_X_2(gyz)*psi4pt + dx_psi4*gyz(i,j,k)
            dx_gzz = DIFF_X_2(gzz)*psi4pt + dx_psi4*gzz(i,j,k)
            dy_gxx = DIFF_Y_2(gxx)*psi4pt + dy_psi4*gxx(i,j,k)
            dy_gxy = DIFF_Y_2(gxy)*psi4pt + dy_psi4*gxy(i,j,k)
            dy_gxz = DIFF_Y_2(gxz)*psi4pt + dy_psi4*gxz(i,j,k)
            dy_gyy = DIFF_Y_2(gyy)*psi4pt + dy_psi4*gyy(i,j,k)
            dy_gyz = DIFF_Y_2(gyz)*psi4pt + dy_psi4*gyz(i,j,k)
            dy_gzz = DIFF_Y_2(gzz)*psi4pt + dy_psi4*gzz(i,j,k)
            dz_gxx = DIFF_Z_2(gxx)*psi4pt + dz_psi4*gxx(i,j,k)
            dz_gxy = DIFF_Z_2(gxy)*psi4pt + dz_psi4*gxy(i,j,k)
            dz_gxz = DIFF_Z_2(gxz)*psi4pt + dz_psi4*gxz(i,j,k)
            dz_gyy = DIFF_Z_2(gyy)*psi4pt + dz_psi4*gyy(i,j,k)
            dz_gyz = DIFF_Z_2(gyz)*psi4pt + dz_psi4*gyz(i,j,k)
            dz_gzz = DIFF_Z_2(gzz)*psi4pt + dz_psi4*gzz(i,j,k)

          else

            dx_gxx = DIFF_X_4(gxx)*psi4pt + dx_psi4*gxx(i,j,k)
            dx_gxy = DIFF_X_4(gxy)*psi4pt + dx_psi4*gxy(i,j,k)
            dx_gxz = DIFF_X_4(gxz)*psi4pt + dx_psi4*gxz(i,j,k)
            dx_gyy = DIFF_X_4(gyy)*psi4pt + dx_psi4*gyy(i,j,k)
            dx_gyz = DIFF_X_4(gyz)*psi4pt + dx_psi4*gyz(i,j,k)
            dx_gzz = DIFF_X_4(gzz)*psi4pt + dx_psi4*gzz(i,j,k)
            dy_gxx = DIFF_Y_4(gxx)*psi4pt + dy_psi4*gxx(i,j,k)
            dy_gxy = DIFF_Y_4(gxy)*psi4pt + dy_psi4*gxy(i,j,k)
            dy_gxz = DIFF_Y_4(gxz)*psi4pt + dy_psi4*gxz(i,j,k)
            dy_gyy = DIFF_Y_4(gyy)*psi4pt + dy_psi4*gyy(i,j,k)
            dy_gyz = DIFF_Y_4(gyz)*psi4pt + dy_psi4*gyz(i,j,k)
            dy_gzz = DIFF_Y_4(gzz)*psi4pt + dy_psi4*gzz(i,j,k)
            dz_gxx = DIFF_Z_4(gxx)*psi4pt + dz_psi4*gxx(i,j,k)
            dz_gxy = DIFF_Z_4(gxy)*psi4pt + dz_psi4*gxy(i,j,k)
            dz_gxz = DIFF_Z_4(gxz)*psi4pt + dz_psi4*gxz(i,j,k)
            dz_gyy = DIFF_Z_4(gyy)*psi4pt + dz_psi4*gyy(i,j,k)
            dz_gyz = DIFF_Z_4(gyz)*psi4pt + dz_psi4*gyz(i,j,k)
            dz_gzz = DIFF_Z_4(gzz)*psi4pt + dz_psi4*gzz(i,j,k)

          end if

        else
!!$          Have to finite difference the conformal factor
           !$OMP CRITICAL
          call CCTK_WARN(0,"At this point in the source terms we should be differencing the conformal factor. Ian has been too lazy to put this bit of code in yet.")
          !$OMP END CRITICAL
        end if

!!$        Contract the shift with the extrinsic curvature

        shiftshiftk = shiftx*shiftx*kxx(i,j,k) + &
                      shifty*shifty*kyy(i,j,k) + &
                      shiftz*shiftz*kzz(i,j,k) + &
             two*(shiftx*shifty*kxy(i,j,k) + &
                  shiftx*shiftz*kxz(i,j,k) + &
                  shifty*shiftz*kyz(i,j,k))

        shiftkx = shiftx*kxx(i,j,k) + shifty*kxy(i,j,k) + shiftz*kxz(i,j,k)
        shiftky = shiftx*kxy(i,j,k) + shifty*kyy(i,j,k) + shiftz*kyz(i,j,k)
        shiftkz = shiftx*kxz(i,j,k) + shifty*kyz(i,j,k) + shiftz*kzz(i,j,k)

!!$        Contract the matter terms with the extrinsic curvature

        sumTK = txx*kxx(i,j,k) + tyy*kyy(i,j,k) + tzz*kzz(i,j,k) &
             + two*(txy*kxy(i,j,k) + txz*kxz(i,j,k) + tyz*kyz(i,j,k))

!!$        Update term for tau

        tau_source = t00* &
             (shiftshiftk - (shiftx*dx_alp + shifty*dy_alp + shiftz*dz_alp) )&
             + t0x*(-dx_alp + two*shiftkx) &
             + t0y*(-dy_alp + two*shiftky) &
             + t0z*(-dz_alp + two*shiftkz) &
             + sumTK

!!$        The following looks very little like the terms in the
!!$        standard papers. Take a look in the ThornGuide to see why
!!$        it is really the same thing.

!!$        Contract the shift with derivatives of the metric

        halfshiftdgx = half*(shiftx*shiftx*dx_gxx + &
             shifty*shifty*dx_gyy + shiftz*shiftz*dx_gzz) + &
             shiftx*shifty*dx_gxy + shiftx*shiftz*dx_gxz + &
             shifty*shiftz*dx_gyz
        halfshiftdgy = half*(shiftx*shiftx*dy_gxx + &
             shifty*shifty*dy_gyy + shiftz*shiftz*dy_gzz) + &
             shiftx*shifty*dy_gxy + shiftx*shiftz*dy_gxz + &
             shifty*shiftz*dy_gyz
        halfshiftdgz = half*(shiftx*shiftx*dz_gxx + &
             shifty*shifty*dz_gyy + shiftz*shiftz*dz_gzz) + &
             shiftx*shifty*dz_gxy + shiftx*shiftz*dz_gxz + &
             shifty*shiftz*dz_gyz

!!$        Contract the matter with derivatives of the metric

        halfTdgx = half*(txx*dx_gxx + tyy*dx_gyy + tzz*dx_gzz) +&
             txy*dx_gxy + txz*dx_gxz + tyz*dx_gyz
        halfTdgy = half*(txx*dy_gxx + tyy*dy_gyy + tzz*dy_gzz) +&
             txy*dy_gxy + txz*dy_gxz + tyz*dy_gyz
        halfTdgz = half*(txx*dz_gxx + tyy*dz_gyy + tzz*dz_gzz) +&
             txy*dz_gxy + txz*dz_gxz + tyz*dz_gyz

        sx_source = t00*&
             (halfshiftdgx - alp(i,j,k)*dx_alp) +&
             t0x*(shiftx*dx_gxx + shifty*dx_gxy + shiftz*dx_gxz) +&
             t0y*(shiftx*dx_gxy + shifty*dx_gyy + shiftz*dx_gyz) +&
             t0z*(shiftx*dx_gxz + shifty*dx_gyz + shiftz*dx_gzz) +&
             halfTdgx + t0lowx*dx_betax + t0lowy*dx_betay + t0lowz*dx_betaz

        sy_source = t00*&
             (halfshiftdgy - alp(i,j,k)*dy_alp) +&
             t0x*(shiftx*dy_gxx + shifty*dy_gxy + shiftz*dy_gxz) +&
             t0y*(shiftx*dy_gxy + shifty*dy_gyy + shiftz*dy_gyz) +&
             t0z*(shiftx*dy_gxz + shifty*dy_gyz + shiftz*dy_gzz) +&
             halfTdgy + t0lowx*dy_betax + t0lowy*dy_betay + t0lowz*dy_betaz

        sz_source = t00*&
             (halfshiftdgz - alp(i,j,k)*dz_alp) +&
             t0x*(shiftx*dz_gxx + shifty*dz_gxy + shiftz*dz_gxz) +&
             t0y*(shiftx*dz_gxy + shifty*dz_gyy + shiftz*dz_gyz) +&
             t0z*(shiftx*dz_gxz + shifty*dz_gyz + shiftz*dz_gzz) +&
             halfTdgz + t0lowx*dz_betax + t0lowy*dz_betay + t0lowz*dz_betaz



        densrhs(i,j,k) = 0.d0
        srhs(i,j,k,1)  = alp(i,j,k)*sqrtdet*sx_source
        srhs(i,j,k,2)  = alp(i,j,k)*sqrtdet*sy_source
        srhs(i,j,k,3)  = alp(i,j,k)*sqrtdet*sz_source
        taurhs(i,j,k) = alp(i,j,k)*sqrtdet*tau_source

        Bxrhs(i,j,k) = 0.d0
        Byrhs(i,j,k) = 0.d0
        Bzrhs(i,j,k) = 0.d0

     enddo
  enddo
enddo
  !$OMP END PARALLEL DO


  !$OMP PARALLEL DO PRIVATE(i, j, k, warnline)
  do k=1 + Spritz_stencil,nz - Spritz_stencil
    do j=1 + Spritz_stencil,ny - Spritz_stencil
      do i=1 + Spritz_stencil,nx - Spritz_stencil
           !check whether rhs is NaN
           if ((densrhs(i,j,k) .ne. densrhs(i,j,k)).or.&
                (taurhs(i,j,k) .ne. taurhs(i,j,k)).or.&
                (srhs(i,j,k,1) .ne. srhs(i,j,k,1)).or.&
                (srhs(i,j,k,2) .ne. srhs(i,j,k,2)).or.&
                (srhs(i,j,k,3) .ne. srhs(i,j,k,3)).or.&
                (Bxrhs(i,j,k) .ne. Bxrhs(i,j,k)).or.&
                (Byrhs(i,j,k) .ne. Byrhs(i,j,k)).or.&
                (Bzrhs(i,j,k) .ne. Bzrhs(i,j,k)))  then
              !$OMP CRITICAL


             ! added local calculations again, as code wasn't checking them ...
              if (conformal_state .eq. 0) then

                localgxx = gxx(i,j,k)
                localgxy = gxy(i,j,k)
                localgxz = gxz(i,j,k)
                localgyy = gyy(i,j,k)
                localgyz = gyz(i,j,k)
                localgzz = gzz(i,j,k)

              else

                psi4pt = psi(i,j,k)**4

                localgxx = psi4pt*gxx(i,j,k)
                localgxy = psi4pt*gxy(i,j,k)
                localgxz = psi4pt*gxz(i,j,k)
                localgyy = psi4pt*gyy(i,j,k)
                localgyz = psi4pt*gyz(i,j,k)
                localgzz = psi4pt*gzz(i,j,k)

              end if


              call SpatialDeterminant(localgxx, localgxy, localgxz, localgyy, localgyz, localgzz, det)

              write(warnline,'(a50,g20.7)') 'rhs is NaN in Spritz Source, densrhs: ',densrhs(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a50,g20.7)') 'rhs is NaN in Spritz Source,  taurhs: ',taurhs(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a50,g20.7)') 'rhs is NaN in Spritz Source,   sxrhs: ',srhs(i,j,k,1)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a50,g20.7)') 'rhs is NaN in Spritz Source,   syrhs: ',srhs(i,j,k,2)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a50,g20.7)') 'rhs is NaN in Spritz Source,   szrhs: ',srhs(i,j,k,3)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a50,g20.7)') 'rhs is NaN in Spritz Source,   Bxrhs: ',Bxrhs(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a50,g20.7)') 'rhs is NaN in Spritz Source,   Byrhs: ',Byrhs(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a50,g20.7)') 'rhs is NaN in Spritz Source,   Bzrhs: ',Bzrhs(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a50,g20.7)') 'rhs is NaN in Spritz Source,   Avec_x_rhs: ',Avec_x_rhs(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a50,g20.7)') 'rhs is NaN in Spritz Source,   Avec_y_rhs: ',Avec_y_rhs(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a50,g20.7)') 'rhs is NaN in Spritz Source,   Avec_z_rhs: ',Avec_z_rhs(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a50,g20.7)') 'rhs is NaN in Spritz Source,   psi_mhdrhs: ',psi_mhdrhs(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a50,g20.7)') 'det: ',det
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') 'gxx=',localgxx
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') 'gxy=',localgxy
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') 'gxz=',localgxz
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') 'gyy=',localgyy
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') 'gyz=',localgyz
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') 'gzz=',localgzz
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') 'dens: ',dens(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') ' tau: ',tau(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') '  sx: ',scon(i,j,k,1)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') '  sy: ',scon(i,j,k,2)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') '  sz: ',scon(i,j,k,3)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') ' dBx: ',dBx(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') ' dBy: ',dBy(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') ' dBz: ',dBz(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') ' rho: ',rho(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') 'press: ',press(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') ' eps: ',eps(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') 'velx: ',velx(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') 'vely: ',vely(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') 'velz: ',velz(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') '  Bx: ',Bx(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') '  By: ',By(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') '  Bz: ',Bz(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') '  Avec_x: ',Avec_x(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') '  Avec_y: ',Avec_y(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') '  Avec_z: ',Avec_z(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') '  psi_mhd: ',psi_mhd(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') 'x=',x(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') 'y=',y(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') 'z=',z(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') 'dx=',dx
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') 'dy=',dy
              call CCTK_WARN(1,warnline)
              write(warnline,'(a20,g20.7)') 'dz=',dz
              call CCTK_WARN(1,warnline)
              call CCTK_WARN(0,"NaN in Spritz_Source (i.e. the values computed by C2P were wrong)")
              !$OMP END CRITICAL
           end if
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

  deallocate(force_spatial_second_order)

end subroutine SourceTerms
