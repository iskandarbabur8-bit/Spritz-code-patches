  /*@@
   @file      Spritz_Analysis_MHD.F90
   @date      
   @author    Bruno Giacomazzo
   @desc 
   Analysis routine for the magnetic field
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

! second order centered finite differencing
#define DIFF_X(VAR,i,j,k)                                                     \
  0.5d0/CCTK_DELTA_SPACE(1)*(VAR(i+1,j,k)-VAR(i-1,j,k))
#define DIFF_Y(VAR,i,j,k)                                                     \
  0.5d0/CCTK_DELTA_SPACE(2)*(VAR(i,j+1,k)-VAR(i,j-1,k))
#define DIFF_Z(VAR,i,j,k)                                                     \
  0.5d0/CCTK_DELTA_SPACE(3)*(VAR(i,j,k+1)-VAR(i,j,k-1))


 /*@@
   @routine    Spritz_divB
   @date       
   @author     Bruno Giacomazzo
   @desc 
   Compute the divergence of the magnetic field
   @enddesc 
   @calls     
   @calledby   
@@*/


subroutine Spritz_divB(CCTK_ARGUMENTS)

  implicit none

!#include "EOS_Base.inc"

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  integer :: i, j, k
  CCTK_REAL :: idx,idy,idz

  CCTK_REAL :: bst,Blowx,Blowy,Blowz,vlowx,vlowy,vlowz,B2big

  CCTK_REAL :: betaxlow, betaylow, betazlow, beta2, bstlow, bsx, bsy, bsz, utlow, ux, uy, uz
  CCTK_REAL :: normalx, normaly, normalz, Exlocal, Eylocal, Ezlocal
  CCTK_REAL :: vtildelowx, vtildelowy, vtildelowz, det
!!$  CCTK_REAL :: velxlow, velylow, velzlow

!!$  CCTK_REAL :: velx_tilde, uxx,uxy,uxz,uyy,uyz,uzz,det,wtot,cs2,dpdrho,dpdeps,presstmp,htmp
!!$  CCTK_REAL :: vely_tilde,velz_tilde,Bv
!!$  CCTK_REAL, dimension(4)::lam

  CCTK_REAL :: Aupx, Aupy, Aupz, uxx, uxy, uxz, uyy, uyz, uzz

  CCTK_REAL :: omegax, omegay, local_angular_velocity, bsxlow, bsylow, bszlow, local_alfven_speed


  idx = 1.0d0/CCTK_DELTA_SPACE(1)
  idy = 1.0d0/CCTK_DELTA_SPACE(2)
  idz = 1.0d0/CCTK_DELTA_SPACE(3)


  !initialize the variables to zero
  divB = 0.0d0
  divB2 = 0.0d0
  divA = 0.0d0
  B_norm = 0.0d0
  A_norm = 0.0d0
  divB_normalized = 0.0d0
  divB2_normalized = 0.0d0
  b2small = 0.0d0
  mhd_press_ratio = 0.0d0
  mhd_energy_ratio = 0.0d0
  
  poynting_flux_x = 0.0d0
  poynting_flux_y = 0.0d0
  poynting_flux_z = 0.0d0
  poynting_flux_norm = 0.0d0
  
  poynting_vector_x = 0.d0
  poynting_vector_y = 0.d0
  poynting_vector_z = 0.d0
  poynting_scalar   = 0.d0

  lambda_MRI_x = 0.0d0
  lambda_MRI_y = 0.0d0
  lambda_MRI_z = 0.0d0

  !$OMP PARALLEL DO PRIVATE(k, j, i, vlowx, vlowy, vlowz, Blowx, Blowy, Blowz, B2big, &
  !$OMP det, uxx, uxy, uxz, uyy, uyz, uzz, Aupx, Aupy, Aupz, bst, betaxlow, betaylow, betazlow, &
  !$OMP beta2, ux, uy, uz, utlow, bsx, bsy, bsz, bstlow, vtildelowx, vtildelowy, vtildelowz, &
  !$OMP Exlocal, Eylocal, Ezlocal, normalx, normaly, normalz, &
  !$OMP omegax, omegay, local_angular_velocity, bsxlow, bsylow, bszlow, local_alfven_speed)
  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)


           !this part can be used to compute the eigenvalues
           !if one wants to check them during the evolution

!!$           w = w_lorentz(i,j,k)
           
           !v_j
           vlowx = gxx(i,j,k)*velx(i,j,k) + gxy(i,j,k)*vely(i,j,k) + gxz(i,j,k)*velz(i,j,k)
           vlowy = gxy(i,j,k)*velx(i,j,k) + gyy(i,j,k)*vely(i,j,k) + gyz(i,j,k)*velz(i,j,k)
           vlowz = gxz(i,j,k)*velx(i,j,k) + gyz(i,j,k)*vely(i,j,k) + gzz(i,j,k)*velz(i,j,k)
           
           
           !B_j
           Blowx = gxx(i,j,k)*Bx(i,j,k) + gxy(i,j,k)*By(i,j,k) + gxz(i,j,k)*Bz(i,j,k)
           Blowy = gxy(i,j,k)*Bx(i,j,k) + gyy(i,j,k)*By(i,j,k) + gyz(i,j,k)*Bz(i,j,k)
           Blowz = gxz(i,j,k)*Bx(i,j,k) + gyz(i,j,k)*By(i,j,k) + gzz(i,j,k)*Bz(i,j,k)

           !B^j B_j
           B2big = Bx(i,j,k)*Blowx + By(i,j,k)*Blowy + Bz(i,j,k)*Blowz

           B_norm(i,j,k) = sqrt(B2big)


           
           !compute the norm of the vector potential
              
           !compute the determinant of the three metric and the upper spatial metric
           call SpatialDeterminant(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),&
                gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det)
           call UpperMetric(uxx, uxy, uxz, uyy, uyz, uzz, &
                det, gxx(i,j,k), gxy(i,j,k), gxz(i,j,k), &
                gyy(i,j,k), gyz(i,j,k), gzz(i,j,k))
           
           !A^j              
           Aupx = uxx*Avec_x(i,j,k) + uxy*Avec_y(i,j,k) + uxz*Avec_z(i,j,k)
           Aupy = uxy*Avec_x(i,j,k) + uyy*Avec_y(i,j,k) + uyz*Avec_z(i,j,k)
           Aupz = uxz*Avec_x(i,j,k) + uyz*Avec_y(i,j,k) + uzz*Avec_z(i,j,k)
           
           !sqrt(A^j A_j)
           A_norm(i,j,k) = sqrt(Avec_x(i,j,k)*Aupx + Avec_y(i,j,k)*Aupy + Avec_z(i,j,k)*Aupz)

!!$           !B^jv_j
!!$           Bv = Bx(i,j,k)*vlowx + By(i,j,k)*vlowy + Bz(i,j,k)*vlowz

           !b^t
           bst = w_lorentz(i,j,k)*(Bx(i,j,k)*vlowx + By(i,j,k)*vlowy + Bz(i,j,k)*vlowz) !this is b^0 * alp
           !b^mu b_mu
           b2small(i,j,k) = (B2big + bst**2)/w_lorentz(i,j,k)**2
                      
           mhd_press_ratio(i,j,k) = b2small(i,j,k)/(2.0*press(i,j,k))
           mhd_energy_ratio(i,j,k) = b2small(i,j,k)/(2.0*rho(i,j,k)) !ratio of magnetic field energy density over rest-mass energy density


           
           !===============COMPUTATION OF THE POYNTING FLUX======================================!
           !=====(see expression in appendix B of Komissarov & Barkov 2009, MNRAS 397, 1153)=====! 
           !=====================================================================================!
           ! S^i = -alp * T^i_0
           ! In Special Relativity it is equivalent to ExB = B^2 v^i - B^i (v^j B_j) [ideal-MHD limit]

!!$           !computing lower components of the velocity
!!$           if (conformal_state .eq. 0) then
!!$              velxlow = gxx(i,j,k)*velx(i,j,k) + gxy(i,j,k)*vely(i,j,k) + gxz(i,j,k)*velz(i,j,k)
!!$              velylow = gxy(i,j,k)*velx(i,j,k) + gyy(i,j,k)*vely(i,j,k) + gyz(i,j,k)*velz(i,j,k)
!!$              velzlow = gxz(i,j,k)*velx(i,j,k) + gyz(i,j,k)*vely(i,j,k) + gzz(i,j,k)*velz(i,j,k)
!!$           else
!!$              velxlow = psi(i,j,k)**4 * (gxx(i,j,k)*velx(i,j,k) + gxy(i,j,k)*vely(i,j,k) + gxz(i,j,k)*velz(i,j,k) )
!!$              velylow = psi(i,j,k)**4 * (gxy(i,j,k)*velx(i,j,k) + gyy(i,j,k)*vely(i,j,k) + gyz(i,j,k)*velz(i,j,k) )
!!$              velzlow = psi(i,j,k)**4 * (gxz(i,j,k)*velx(i,j,k) + gyz(i,j,k)*vely(i,j,k) + gzz(i,j,k)*velz(i,j,k) )
!!$           end if

           !computing lower components of the shift
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
              beta2 = 0.0d0
           end if

           ! Calculate u^i
           ux = (velx(i,j,k)-betax(i,j,k)/alp(i,j,k))*w_lorentz(i,j,k)
           uy = (vely(i,j,k)-betay(i,j,k)/alp(i,j,k))*w_lorentz(i,j,k)
           uz = (velz(i,j,k)-betaz(i,j,k)/alp(i,j,k))*w_lorentz(i,j,k)

           ! Calculate u_0
           utlow = (-alp(i,j,k) + velx(i,j,k)*betaxlow + vely(i,j,k)*betaylow + velz(i,j,k)*betazlow)*w_lorentz(i,j,k)

           ! Calculate the 4 vector b of the magnetic field
           bst=w_lorentz(i,j,k)*(Bx(i,j,k)*vlowx + By(i,j,k)*vlowy + Bz(i,j,k)*vlowz)/alp(i,j,k)
           bsx=(Bx(i,j,k)+alp(i,j,k)*bst*w_lorentz(i,j,k)*ux)/(w_lorentz(i,j,k))
           bsy=(By(i,j,k)+alp(i,j,k)*bst*w_lorentz(i,j,k)*uy)/(w_lorentz(i,j,k))
           bsz=(Bz(i,j,k)+alp(i,j,k)*bst*w_lorentz(i,j,k)*uz)/(w_lorentz(i,j,k))

           ! Calculate b_0
           bstlow = bst*(-alp(i,j,k)**2+beta2)+bsx*betaxlow + bsy*betaylow + bsz*betazlow

           !S^i = -alp*T^i_0 = -alp*(b^2 u^i u_0 + b^2/2 g^i_0 - b^i b_0) (being T only the MHD part of the energy-momentum tensor)
           poynting_flux_x(i,j,k) = -alp(i,j,k)*(b2small(i,j,k)*ux*utlow - bsx*bstlow)
           poynting_flux_y(i,j,k) = -alp(i,j,k)*(b2small(i,j,k)*uy*utlow - bsy*bstlow)
           poynting_flux_z(i,j,k) = -alp(i,j,k)*(b2small(i,j,k)*uz*utlow - bsz*bstlow)
           !S^2
           poynting_flux_norm(i,j,k) = gxx(i,j,k)*poynting_flux_x(i,j,k)*poynting_flux_x(i,j,k) + &
                2.0d0*gxy(i,j,k)*poynting_flux_x(i,j,k)*poynting_flux_y(i,j,k) + &
                2.0d0*gxz(i,j,k)*poynting_flux_x(i,j,k)*poynting_flux_z(i,j,k) + &
                gyy(i,j,k)*poynting_flux_y(i,j,k)*poynting_flux_y(i,j,k) + &
                2.0d0*gyz(i,j,k)*poynting_flux_y(i,j,k)*poynting_flux_z(i,j,k) + &
                gzz(i,j,k)*poynting_flux_z(i,j,k)*poynting_flux_z(i,j,k)
           
           
           !==================ALTERNATIVE COMPUTATION OF POYNTING FLUX=======================================!
           !! computation of the electric field E^i = - (vtilde_j x B_k)/sqrt(det)
           !! vtilde_i = alp*v_i-beta_i
           vtildelowx = alp(i,j,k)*vlowx - betaxlow
           vtildelowy = alp(i,j,k)*vlowy - betaylow
           vtildelowz = alp(i,j,k)*vlowz - betazlow
           Exlocal = -1.0d0*(vtildelowy*Blowz - vtildelowz*Blowy) !E^x*sqrt(det)
           Eylocal = -1.0d0*(vtildelowz*Blowx - vtildelowx*Blowz) !E^y*sqrt(det)
           Ezlocal = -1.0d0*(vtildelowx*Blowy - vtildelowy*Blowx) !E^z*sqrt(det)


!!$           call SpatialDeterminant(gxx(i,j,k),gxy(i,j,k),&
!!$                gxz(i,j,k),gyy(i,j,k), gyz(i,j,k),&
!!$                gzz(i,j,k),det)

           !! Poynting vector S_i (covariant) in a cartesian coordinate basis
           poynting_vector_x(i,j,k) = (Eylocal*Bz(i,j,k) - Ezlocal*By(i,j,k))!*sqrt(det) is contained in the def of E^i
           poynting_vector_y(i,j,k) = (Ezlocal*Bx(i,j,k) - Exlocal*Bz(i,j,k))!*sqrt(det) is contained in the def of E^i
           poynting_vector_z(i,j,k) = (Exlocal*By(i,j,k) - Eylocal*Bx(i,j,k))!*sqrt(det) is contained in the def of E^i
           
           !! Normal vector n^i of a sphere in a cartesian coordinate basis
           normalx = x(i,j,k)/r(i,j,k)
           normaly = y(i,j,k)/r(i,j,k)
           normalz = z(i,j,k)/r(i,j,k)
           
           !! Poynting scalar (S_i * n^i) in a cartesian coordinate basis
           poynting_scalar(i,j,k) = (poynting_vector_x(i,j,k) * normalx + &
                poynting_vector_y(i,j,k) * normaly + &
                poynting_vector_z(i,j,k) * normalz)
           



           !==================COMPUTATION OF DIVB AND DIVB2=======================================!
           !we do not compute divB and divB2 where we do not have enough points available for them
           !======================================================================================!

           if( (i+1>cctk_lsh(1)).or.(i-1<1).or.(j+1>cctk_lsh(2)).or.(j-1<1).or.(k+1>cctk_lsh(3)).or.(k-1<1)) cycle

           !when using the vector potential
           if( (i+2>cctk_lsh(1)).or.(i-2<1).or.(j+2>cctk_lsh(2)).or.(j-2<1).or.(k+2>cctk_lsh(3)).or.(k-2<1)) cycle
           


           divB(i,j,k) = 0.5d0*idx*(dBx(i+1,j,k)-dBx(i-1,j,k)) + &
                0.5d0*idy*(dBy(i,j+1,k)-dBy(i,j-1,k)) + &
                0.5d0*idz*(dBz(i,j,k+1)-dBz(i,j,k-1))

           divA(i,j,k) = 0.5d0*idx*(Avec_x(i+1,j,k)-Avec_x(i-1,j,k)) + &
                0.5d0*idy*(Avec_y(i,j+1,k)-Avec_y(i,j-1,k)) + &
                0.5d0*idz*(Avec_z(i,j,k+1)-Avec_z(i,j,k-1))


           !The following divB is conserved if divb is conserved,
           !where b is the magnetic field computed on the cell faces 
           !while B is computed at the cell center
           !note: Bx(i,j,k) = (bx(i+1/2,j,k)+bx(i-1/2,j,k))/2
           !divb(i,j,k) = (bx(i+1/2,j,k)-bx(i-1/2,j,k))/dx + (by(i,j+1/2,k)-bx(i,j-1/2,k))/dy +
           ! (bz(i,j,k+1/2)-bz(i,j,k-1/2))/dz

           divB2(i,j,k) = 0.25*idx*(dBx(i+1,j,k)-dBx(i,j,k) + dBx(i+1,j+1,k)-dBx(i,j+1,k) + &
                dBx(i+1,j,k+1)-dBx(i,j,k+1) + dBx(i+1,j+1,k+1)-dBx(i,j+1,k+1)) + &
                0.25*idy*(dBy(i,j+1,k)-dBy(i,j,k) + dBy(i+1,j+1,k)-dBy(i+1,j,k) + &
                dBy(i,j+1,k+1)-dBy(i,j,k+1) + dBy(i+1,j+1,k+1)-dBy(i+1,j,k+1)) + &
                0.25*idz*(dBz(i,j,k+1)-dBz(i,j,k) + dBz(i+1,j,k+1)-dBz(i+1,j,k) + &
                dBz(i,j+1,k+1)-dBz(i,j+1,k) + dBz(i+1,j+1,k+1)-dBz(i+1,j+1,k))



           if (B_norm(i,j,k)>0.0d0) then
              divB_normalized(i,j,k) = divB(i,j,k)/B_norm(i,j,k)
              divB2_normalized(i,j,k) = divB2(i,j,k)/B_norm(i,j,k)
           else
              divB_normalized(i,j,k) = 0.0d0
              divB2_normalized(i,j,k) = 0.0d0
           end if

           !Compute Lambda MRI (see McKinney, Tchekhovskoy, and Blandford 2012, http://arxiv.org/abs/1201.4163)
           !Calculates the angular velocity omega in Cartesian coordinates.
           omegax = alp(i,j,k)*velx(i,j,k) - betax(i,j,k)
           omegay = alp(i,j,k)*vely(i,j,k) - betay(i,j,k)
           !Transforms the angular velocity to cylindrical coordinates.
           if ((x(i,j,k)**2+y(i,j,k)**2)>0.0d0) then
              local_angular_velocity = (x(i,j,k)*omegay-y(i,j,k)*omegax)/(x(i,j,k)**2+y(i,j,k)**2)
           else
              local_angular_velocity = 0.0d0
           end if
           !Calculate b_i
           bsxlow = bst*betaxlow + bsx*gxx(i,j,k) + bsy*gxy(i,j,k) + bsz*gxz(i,j,k)
           bsylow = bst*betaylow + bsx*gxy(i,j,k) + bsy*gyy(i,j,k) + bsz*gyz(i,j,k)
           bszlow = bst*betazlow + bsx*gxz(i,j,k) + bsy*gyz(i,j,k) + bsz*gzz(i,j,k)
           !Calculate Va_x and lambda_x
           local_alfven_speed = sqrt(bsx*bsxlow/(rho(i,j,k)*(1.0d0+eps(i,j,k))+press(i,j,k)+b2small(i,j,k)))
           if ((local_alfven_speed>epsilon(local_alfven_speed)).and.(abs(local_angular_velocity)>1.0d-15)) &
                lambda_MRI_x(i,j,k) = 2.0d0*3.1416*abs(local_alfven_speed/local_angular_velocity)
           !Calculate Va_y and lambda_y
           local_alfven_speed = sqrt(bsy*bsylow/(rho(i,j,k)*(1.0d0+eps(i,j,k))+press(i,j,k)+b2small(i,j,k)))
           if ((local_alfven_speed>epsilon(local_alfven_speed)).and.(abs(local_angular_velocity)>1.0d-15)) &
                lambda_MRI_y(i,j,k) = 2.0d0*3.1416*abs(local_alfven_speed/local_angular_velocity)
           !Calculate Va_z and lambda_z
           local_alfven_speed = sqrt(bsz*bszlow/(rho(i,j,k)*(1.0d0+eps(i,j,k))+press(i,j,k)+b2small(i,j,k)))
           if ((local_alfven_speed>epsilon(local_alfven_speed)).and.(abs(local_angular_velocity)>1.0d-15)) &
                lambda_MRI_z(i,j,k) = 2.0d0*3.1416*abs(local_alfven_speed/local_angular_velocity)

        end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine Spritz_divB


subroutine Spritz_B_Spherical(CCTK_ARGUMENTS)
! Here we will decompose the magnetic field vector into spherical coordinates.
! Following this we will compute the norm of the poloidal component. The 
! toroidal component is merely the phi component of the vector.
! Convention: p = (rho, theta, phi)
!       rho = distance from origin to sphere containing point
!       theta = angle between z axis down great circle to point
!       phi = angle in xy plane from the x axis to the plane containing the
!               great circle
! The vectors Bmu correspond to the mu component of the magnetic field, while
!  the components Blowmu correspond to the mu component of the dual vector
!  field. Mu can be an element of (x,y,z,rho,theta,phi).
implicit none

  DECLARE_CCTK_ARGUMENTS
  integer :: i, j, k
  ! I am prefacing some local variables with my initials, since I expect
  !     it not to generate conflicts with other variables.
  ! Note the use of attempting camel casing for new variable names.
  !     workingVector is the junk vector for the new fields.
  CCTK_REAL :: at_r, at_rho, Blowx, Blowy, Blowz, Blowphi
  CCTK_REAL :: toroidalAxis(3), toroidalOrigin(3), toroidalRotation(3,3)
  CCTK_REAL :: newX, newY, newZ, Bupx, Bupy, Bupz, workingVector(3)
  CCTK_REAL :: det, B2big, B_norm_local

  !initialization
  Spritz_mean_B_Tor_numerator_temps = 0.0d0
  Spritz_mean_B_Pol_numerator_temps = 0.0d0
  Spritz_mean_B_Tor_denominator_temps = 0.0d0
  
  ! trivial initialisation of cartesian coordinate system transform
  ! TODO: Implement the following vectors so that they can in fact be used
  !     as parameters or variables.
  toroidalAxis = 0.0e0
  toroidalAxis(3) = 1.0e0
  toroidalOrigin = 0.0e0
  ! construct the rotation matrix
  call rotationMatrix(toroidalAxis,toroidalRotation)

  ! With the change of coordinates well defined we proceed to transform
  !     the vectors.

  !$OMP PARALLEL DO PRIVATE(k, j, i, workingVector, newX, newY, newZ, &
  !$OMP Blowx, Blowy, Blowz, Bupx, Bupy, Bupz, at_r, at_rho, Blowphi, det, &
  !$OMP B2big, B_norm_local)
  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
        ! First, we will transfer to the new coordinates.
        workingVector(1) = x(i,j,k)
        workingVector(2) = y(i,j,k)
        workingVector(3) = z(i,j,k)
        workingVector = workingVector - toroidalOrigin
        workingVector = matmul(toroidalRotation,workingVector)
        newX = workingVector(1)
        newY = workingVector(2)
        newZ = workingVector(3)

        ! Need the down components of the normal field. Pulled from subroutine 
        !       above by Bruno.
        workingVector(1) = gxx(i,j,k)*Bx(i,j,k) + gxy(i,j,k)*By(i,j,k) &
                + gxz(i,j,k)*Bz(i,j,k) ! this is Blowx old
        workingVector(2) = gxy(i,j,k)*Bx(i,j,k) + gyy(i,j,k)*By(i,j,k) &
                + gyz(i,j,k)*Bz(i,j,k) ! this is Blowy old
        workingVector(3) = gxz(i,j,k)*Bx(i,j,k) + gyz(i,j,k)*By(i,j,k) &
                + gzz(i,j,k)*Bz(i,j,k) ! this is Blowz old
        workingVector = matmul(toroidalRotation,workingVector)
        Blowx = workingVector(1)
        Blowy = workingVector(2)
        Blowz = workingVector(3)

        !B^j B_j
        B2big = Bx(i,j,k)*Blowx + By(i,j,k)*Blowy + Bz(i,j,k)*Blowz
        B_norm_local = sqrt(B2big)

        ! Now we make the up components of the new magnetic field.
        workingVector(1) = Bx(i,j,k)
        workingVector(2) = By(i,j,k)
        workingVector(3) = Bz(i,j,k)
        workingVector = matmul(toroidalRotation,workingVector)
        Bupx = workingVector(1)
        Bupy = workingVector(2)
        Bupz = workingVector(3)

        ! Some combinations of the coordinates that we will frequently need.
        at_r = sqrt(newX**2 + newY**2)
        at_rho = sqrt(at_r**2 + newZ**2)
        ! Need to check for coordinate singularities.
        if (abs(at_r)<=epsilon(at_r) .or. abs(at_rho)<=epsilon(at_rho)) then
                ! Singular behavior. Just set them to zero.
                Bphi(i,j,k) = 0.0e0
                Btheta(i,j,k) = 0.0e0
                Brho(i,j,k) = 0.0e0
                Blowphi = 0.0e0
                !Blowtheta = 0.0e0
                !Blowrho = 0.0e0
                B_Poloidal_norm(i,j,k) = 0.0e0
                B_Toroidal_norm(i,j,k) = 0.0e0
        else
        ! Phi component, up and down
        Bphi(i,j,k) = Bupy*newX - Bupx*newY
        Bphi(i,j,k) = Bphi(i,j,k)/at_r**2
        Blowphi = Blowy*newX - Blowx*newY
        ! further covariant quantities are shown commented for clarity if
        !       one wishes to extend analysis

        ! Rho component up and down
        Brho(i,j,k) = Bupx*newX +Bupy*newY +Bupz*newZ
        Brho(i,j,k) = Brho(i,j,k)/at_rho
        ! Blowrho = Blowx*newX + Blowy*newY
        ! Blowrho = (Blowrho(i,j,k) + Blowz*newZ)/at_rho

        ! Theta component, up and down.
        Btheta(i,j,k) = (Bupx*newX + Bupy*newY)*newZ
        Btheta(i,j,k) = (Btheta(i,j,k)/at_r - Bupz*at_r)/at_rho**2
        ! Blowtheta = Blowx*newX + Blowy*newY
        ! Blowtheta = Blowtheta(i,j,k)*newZ/at_r - Blowz*at_r

        ! Norm of the poloidal component of the magnetic field.
        B_Poloidal_norm(i,j,k)=sqrt(abs(B_norm_local**2-Blowphi*Bphi(i,j,k)))
        B_Toroidal_norm(i,j,k)=sqrt(abs(Blowphi*Bphi(i,j,k)))
        end if

        !mean magnetic field value
        !Int(rho*B*dV)/Int(rho*dV)
        call SpatialDeterminant(gxx(i,j,k),gxy(i,j,k),&
             gxz(i,j,k),gyy(i,j,k), gyz(i,j,k),&
             gzz(i,j,k),det)
        Spritz_mean_B_Tor_numerator_temps(i,j,k) = (rho(i,j,k)*B_Toroidal_norm(i,j,k)*w_lorentz(i,j,k)*sqrt(det))
        Spritz_mean_B_Pol_numerator_temps(i,j,k) = (rho(i,j,k)*B_Poloidal_norm(i,j,k)*w_lorentz(i,j,k)*sqrt(det))
        Spritz_mean_B_Tor_denominator_temps(i,j,k) = (rho(i,j,k)*w_lorentz(i,j,k)*sqrt(det))

        end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine Spritz_B_Spherical

! Construct a rotation matrix which takes a prescribed toroidalAxis and rotates
!       it onto the z axis.
subroutine rotationMatrix(toroidalAxis,toroidalRotation)
  implicit none
  CCTK_REAL :: toroidalAxis(3), toroidalRotation(3,3), rAxis, rhoAxis

  ! Initialise the rotation matrix to the identity
  toroidalRotation = 0.0e0
  toroidalRotation(1,1) = 1.0e0
  toroidalRotation(2,2) = 1.0e0
  toroidalRotation(3,3) = 1.0e0
  ! Set up some easier to use variables for constructing the rotation matrix
  rAxis = sqrt( toroidalAxis(1)**2 + toroidalAxis(2)**2 )
  rhoAxis = sqrt( rAxis**2 + toroidalAxis(3)**2 )
  ! Construct the rotation matrix, test for singular behavior.
  if ( rAxis .ne. 0.0e0 .and. rhoAxis .ne. 0.0e0) then
        ! first row
        toroidalRotation(1,1) = toroidalAxis(1)*toroidalAxis(3)/rAxis/rhoAxis
        toroidalRotation(1,2) = toroidalAxis(2)*toroidalAxis(3)/rAxis/rhoAxis
        toroidalRotation(1,3) = -rAxis/rhoAxis
        ! Second row
        toroidalRotation(2,1) = -toroidalAxis(2)/rAxis
        toroidalRotation(2,2) = toroidalAxis(1)/rAxis
        toroidalRotation(2,3) = 0.0e0
        ! Third row
        toroidalRotation(3,1) = toroidalAxis(1)/rhoAxis
        toroidalRotation(3,2) = toroidalAxis(2)/rhoAxis
        toroidalRotation(3,3) = toroidalAxis(3)/rhoAxis
  else if (rhoAxis .eq. 0.0e0) then
          ! somebody gave us a bad vector for the axis, default to identity
          write(0,*) 'toroidal field decomposition: error: vanishing axis'
  end if
end subroutine rotationMatrix


 /*@@
   @routine    Spritz_Magnetic_Energy_Global
   @date       
   @author     Bruno Giacomazzo
   @desc 

   Calculates the integral of the magnetic energy
   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Spritz_B_Spherical_Global(CCTK_ARGUMENTS)
      
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_INT :: ierr, Reduction_Handle, VarIndex

  call CCTK_ReductionHandle(Reduction_Handle, "sum")
 
  !Mean Toroidal Magnetic Field Value
  !Numerator=Int(rho*B*dV)
  call CCTK_VarIndex(VarIndex, "Spritz::Spritz_mean_B_Tor_numerator_temps")
  call CCTK_Reduce(ierr, cctkGH, -1, Reduction_Handle, &
       1, CCTK_VARIABLE_REAL, Spritz_mean_B_Tor_numerator, 1, VarIndex)
  if (ierr .ne. 0) then
     call CCTK_WARN(0, "Failed to integrate Spritz_mean_B_Tor_numerator_temps")
  end if
  Spritz_mean_B_Tor_numerator = Spritz_mean_B_Tor_numerator * grid_spacing_product_mhd

  !Mean Poloidal Magnetic Field Value
  !Numerator=Int(rho*B*dV)
  call CCTK_VarIndex(VarIndex, "Spritz::Spritz_mean_B_Pol_numerator_temps")
  call CCTK_Reduce(ierr, cctkGH, -1, Reduction_Handle, &
       1, CCTK_VARIABLE_REAL, Spritz_mean_B_Pol_numerator, 1, VarIndex)
  if (ierr .ne. 0) then
     call CCTK_WARN(0, "Failed to integrate Spritz_mean_B_Pol_numerator_temps")
  end if
  Spritz_mean_B_Pol_numerator = Spritz_mean_B_Pol_numerator * grid_spacing_product_mhd

  !Denominator=Int(rho*dV)
  call CCTK_VarIndex(VarIndex, "Spritz::Spritz_mean_B_Tor_denominator_temps")
  call CCTK_Reduce(ierr, cctkGH, -1, Reduction_Handle, &
       1, CCTK_VARIABLE_REAL, Spritz_mean_B_Tor_denominator, 1, VarIndex)
  if (ierr .ne. 0) then
     call CCTK_WARN(0, "Failed to integrate Spritz_mean_B_Tor_denominator_temps")
  end if
  Spritz_mean_B_Tor_denominator = Spritz_mean_B_Tor_denominator * grid_spacing_product_mhd

  !Mean Magnetic Field Values (poloidal and toroidal)
  Spritz_mean_B_Tor =  Spritz_mean_B_Tor_numerator/Spritz_mean_B_Tor_denominator
  Spritz_mean_B_Pol =  Spritz_mean_B_Pol_numerator/Spritz_mean_B_Tor_denominator


end subroutine Spritz_B_Spherical_Global




 /*@@
   @routine    Spritz_Magnetic_Energy_Local
   @date       
   @author     
   @desc 

   Compute the Magnetic Energy
   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Spritz_Magnetic_Energy_Local(CCTK_ARGUMENTS)
      
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: i, j, k, lst, rst1, rst2, rst3, nchars
  character(len=3) :: avoid_origin_parameter  
  CCTK_REAL :: bst,Blowx,Blowy,Blowz,vlowx,vlowy,vlowz,B2big,b2smalltemp,det,v2local
  CCTK_REAL :: curlv_x, curlv_y, curlv_z, curlv

  ! initialize all the temporary grid-variables 
  Spritz_EMenergy_temps = 0.d0
  Spritz_kinetic_energy_temps = 0.0d0
  Spritz_mean_B_numerator_temps = 0.0d0
  Spritz_mean_B_denominator_temps = 0.0d0
  Spritz_turbulent_energy_temps = 0.0d0
  Spritz_turbulent_energy_amplification_temps = 0.0d0
  
  !for short hand, right and left modified stencil values
  
  lst  = 1+Spritz_stencil
  rst1 = cctk_lsh(1) - Spritz_stencil
  rst2 = cctk_lsh(2) - Spritz_stencil
  rst3 = cctk_lsh(3) - Spritz_stencil

  !$OMP PARALLEL DO PRIVATE(k, j, i, det, vlowx, vlowy, vlowz, Blowx, Blowy, Blowz, &
  !$OMP B2big, bst, b2smalltemp, v2local, curlv_x, curlv_y, curlv_z, curlv)
  do k = lst, rst3
     do j = lst, rst2
        do i = lst, rst1



           call SpatialDeterminant(gxx(i,j,k),gxy(i,j,k),&
                gxz(i,j,k),gyy(i,j,k), gyz(i,j,k),&
                gzz(i,j,k),det)

              !VERY IMPORTANT NOTE:
              !B_norm and b2small are computed only when the output is requested for one of the
              !variables in Spritz_mhd_analysis and this happens only at ANALYSIS while for
              !Carpet we need to compute integral at POSTSTEP otherwise it would not work.
              !This is why I have to recompute quantities that are already computed in Spritz_divB :-(

              !v_j
              vlowx = gxx(i,j,k)*velx(i,j,k) + gxy(i,j,k)*vely(i,j,k) + gxz(i,j,k)*velz(i,j,k)
              vlowy = gxy(i,j,k)*velx(i,j,k) + gyy(i,j,k)*vely(i,j,k) + gyz(i,j,k)*velz(i,j,k)
              vlowz = gxz(i,j,k)*velx(i,j,k) + gyz(i,j,k)*vely(i,j,k) + gzz(i,j,k)*velz(i,j,k)
              
              !B_j
              Blowx = gxx(i,j,k)*Bx(i,j,k) + gxy(i,j,k)*By(i,j,k) + gxz(i,j,k)*Bz(i,j,k)
              Blowy = gxy(i,j,k)*Bx(i,j,k) + gyy(i,j,k)*By(i,j,k) + gyz(i,j,k)*Bz(i,j,k)
              Blowz = gxz(i,j,k)*Bx(i,j,k) + gyz(i,j,k)*By(i,j,k) + gzz(i,j,k)*Bz(i,j,k)
              
              !B^j B_j
              B2big = Bx(i,j,k)*Blowx + By(i,j,k)*Blowy + Bz(i,j,k)*Blowz

              !b^t
              bst = w_lorentz(i,j,k)*(Bx(i,j,k)*vlowx + By(i,j,k)*vlowy + Bz(i,j,k)*vlowz) !this is b^0 * alp

              !b^mu b_mu
              b2smalltemp = (B2big + bst**2)/w_lorentz(i,j,k)**2

              !This is n_{\mu} n_{\nu} T^{\mu \nu}_{EM}*sqrt(gamma)
              !Note that n_{\mu} n_{\nu} T^{\mu \nu}_{EM} is the EM energy measured in the Eulerian frame
              ! so one needs to divide by W to transform in the comoving frame.
              !Once in the comoving frame, the proper volume is given by sqrt(gamma)*W*dx^3
              !So the Lorentz factors cancel each others and I need to multiply only by sqrt(gamma)
              !see (astro-ph/0605331)
              Spritz_EMenergy_temps(i,j,k) = (B2big - 0.5d0*b2smalltemp)*sqrt(det)

              !Kinetic Energy (Newtonian expression)
              !0.5*rho*h*v^2*dV, with dV=w_lorentz*sqrt(det)*dx^3
              v2local = velx(i,j,k)*vlowx + vely(i,j,k)*vlowy + velz(i,j,k)*vlowz !v^i v_i
              Spritz_kinetic_energy_temps(i,j,k) = 0.5d0*(rho(i,j,k)+rho(i,j,k)*eps(i,j,k)+press(i,j,k))*v2local &
                   *w_lorentz(i,j,k)*sqrt(det)

              !mean magnetic field value
              !Int(rho*B*dV)/Int(rho*dV)
              Spritz_mean_B_numerator_temps(i,j,k) = (rho(i,j,k)*sqrt(B2big)*w_lorentz(i,j,k)*sqrt(det))
              Spritz_mean_B_denominator_temps(i,j,k) = (rho(i,j,k)*w_lorentz(i,j,k)*sqrt(det))


              !Turbulent Energy
              ! (see appendix of http://arxiv.org/abs/1302.7306)
              if(Spritz_use_zrake_macfadyen_amplification.ne.0) then
                 !NOTE: eps_KT(i,j,k) = zrake_c2*delta_w
                 Spritz_turbulent_energy_temps(i,j,k) = (eps_KT(i,j,k)/zrake_c2)*w_lorentz(i,j,k)*sqrt(det)

                 !this is the turbulent energy contained in the regions where we amplify the magnetic field
                 !compute curl(v)
                 curlv_x = DIFF_Y(velz,i,j,k)-DIFF_Z(vely,i,j,k)
                 curlv_y = DIFF_Z(velx,i,j,k)-DIFF_X(velz,i,j,k)
                 curlv_z = DIFF_X(vely,i,j,k)-DIFF_Y(velx,i,j,k)
                 curlv = sqrt(curlv_x**2 + curlv_y**2 + curlv_z**2)

                 if((eps_KT(i,j,k)>0.0d0).and.(curlv>0.5).and.(rho(i,j,k)>1.0d4*Spritz_rho_min)) then
                    Spritz_turbulent_energy_amplification_temps(i,j,k) = (eps_KT(i,j,k)/zrake_c2)*w_lorentz(i,j,k)*sqrt(det)
                 else
                    Spritz_turbulent_energy_amplification_temps(i,j,k) = 0.0d0
                 end if
              end if


!!$           end if
           
        end do
     end do
  end do
  !$OMP END PARALLEL DO
    
  ! Carpet does the following itself, but pugh does not.
  
  if (CCTK_IsThornActive("PUGH") /= 0) then

     call CCTK_WARN(0,"PUGH is not currently supported in Spritz_Analysis_MHD")
     
     ! Recompute the integrand on the symmetry and physical boundaries: devide by 2 on the boundary faces
     ! cctk_lbnd values start from zero, so must add one to get the total number of points.
     
     ! find the value of the "avoid_origin"  parameter
     call CCTK_ParameterValString(nchars,"avoid_origin","CartGrid3D",avoid_origin_parameter)
     
     ! if  avoid_origin = yes (default), then do not devide by two. 
     ! This gives the correct result for the trapezoidal integration rule on the symmetry boundaries, 
     ! while (even if it is wrong not to devide by two the points on the faces that are 
     ! physical outer boundaries) it should not spoil the result.
    
     if (avoid_origin_parameter == "no ") then

      
        if (cctk_lbnd(3) == 0) then
           
           do j = lst, rst2       
              do i = lst, rst1
                 
                 Spritz_EMenergy_temps(i,j,lst) = Spritz_EMenergy_temps(i,j,lst) * 0.5d0
                 
              end do
           end do
           
        end if
        
        
        if (cctk_ubnd(3)+1 == cctk_gsh(3)) then
           
           do j = lst, rst2
              do i = lst, rst1  
                 
                 Spritz_EMenergy_temps(i,j,rst3) = Spritz_EMenergy_temps(i,j,rst3) * 0.5d0
                 
              end do
           end do
           
        end if
        
        
        if (cctk_lbnd(2) == 0) then
           
           do k = lst, rst3
              do i = lst, rst1
                 
                 Spritz_EMenergy_temps(i,lst,k) = Spritz_EMenergy_temps(i,lst,k) * 0.5d0
                 
              end do
           end do
           
        end if
        
        
        if (cctk_ubnd(2)+1 == cctk_gsh(2)) then
           
           do k = lst, rst3
              do i = lst, rst1
                 
                 Spritz_EMenergy_temps(i,rst2,k) = Spritz_EMenergy_temps(i,rst2,k) * 0.5d0
                 
              end do
           end do
           
        end if
        
        
        if (cctk_lbnd(1) == 0) then
           
           do k = lst, rst3
              do j = lst, rst2
                 
                 Spritz_EMenergy_temps(lst,j,k) = Spritz_EMenergy_temps(lst,j,k) * 0.5d0
                 
              end do
           end do
           
        end if
        
        
        if (cctk_ubnd(1)+1 == cctk_gsh(1)) then
           
           do k = lst, rst3
              do j = lst, rst2
                 
                 Spritz_EMenergy_temps(rst1,j,k) = Spritz_EMenergy_temps(rst1,j,k) * 0.5d0
                 
              end do
           end do
          
        end if
        
     end if !end avoid_origin
      
  end if !end PUGH
    

  grid_spacing_product_mhd = product (cctk_delta_space) !must be computed in local mode and stored
  
end subroutine Spritz_Magnetic_Energy_Local


 /*@@
   @routine    Spritz_Magnetic_Energy_Global
   @date       
   @author     Bruno Giacomazzo
   @desc 

   Calculates the integral of the magnetic energy
   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Spritz_Magnetic_Energy_Global(CCTK_ARGUMENTS)
      
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_INT :: ierr, Reduction_Handle, VarIndex

  call CCTK_ReductionHandle(Reduction_Handle, "sum")
 
  call CCTK_VarIndex(VarIndex, "Spritz::Spritz_EMenergy_temps")
  call CCTK_Reduce(ierr, cctkGH, -1, Reduction_Handle, &
       1, CCTK_VARIABLE_REAL, Spritz_total_magnetic_energy, 1, VarIndex)
  if (ierr .ne. 0) then
     call CCTK_WARN(0, "Failed to integrate Spritz_EMenergy_temps")
  end if
  
  Spritz_total_magnetic_energy = Spritz_total_magnetic_energy * grid_spacing_product_mhd


  !Kinetic Energy
  call CCTK_VarIndex(VarIndex, "Spritz::Spritz_kinetic_energy_temps")
  call CCTK_Reduce(ierr, cctkGH, -1, Reduction_Handle, &
       1, CCTK_VARIABLE_REAL, Spritz_total_kinetic_energy, 1, VarIndex)
  if (ierr .ne. 0) then
     call CCTK_WARN(0, "Failed to integrate Spritz_kinetic_energy_temps")
  end if
  Spritz_total_kinetic_energy = Spritz_total_kinetic_energy * grid_spacing_product_mhd

  !Mean Magnetic Field Value
  !Numerator=Int(rho*B*dV)
  call CCTK_VarIndex(VarIndex, "Spritz::Spritz_mean_B_numerator_temps")
  call CCTK_Reduce(ierr, cctkGH, -1, Reduction_Handle, &
       1, CCTK_VARIABLE_REAL, Spritz_mean_B_numerator, 1, VarIndex)
  if (ierr .ne. 0) then
     call CCTK_WARN(0, "Failed to integrate Spritz_mean_B_numerator_temps")
  end if
  Spritz_mean_B_numerator = Spritz_mean_B_numerator * grid_spacing_product_mhd

  !Denominator=Int(rho*dV)
  call CCTK_VarIndex(VarIndex, "Spritz::Spritz_mean_B_denominator_temps")
  call CCTK_Reduce(ierr, cctkGH, -1, Reduction_Handle, &
       1, CCTK_VARIABLE_REAL, Spritz_mean_B_denominator, 1, VarIndex)
  if (ierr .ne. 0) then
     call CCTK_WARN(0, "Failed to integrate Spritz_mean_B_denominator_temps")
  end if
  Spritz_mean_B_denominator = Spritz_mean_B_denominator * grid_spacing_product_mhd

  !Mean Magnetic Field Value
  Spritz_mean_B =  Spritz_mean_B_numerator/Spritz_mean_B_denominator


  !Turbulent Energy
  call CCTK_VarIndex(VarIndex, "Spritz::Spritz_turbulent_energy_temps")
  call CCTK_Reduce(ierr, cctkGH, -1, Reduction_Handle, &
       1, CCTK_VARIABLE_REAL, Spritz_turbulent_energy, 1, VarIndex)
  if (ierr .ne. 0) then
     call CCTK_WARN(0, "Failed to integrate Spritz_turbulent_energy_temps")
  end if
  Spritz_turbulent_energy = Spritz_turbulent_energy * grid_spacing_product_mhd

  !Turbulent Energy
  call CCTK_VarIndex(VarIndex, "Spritz::Spritz_turbulent_energy_amplification_temps")
  call CCTK_Reduce(ierr, cctkGH, -1, Reduction_Handle, &
       1, CCTK_VARIABLE_REAL, Spritz_turbulent_energy_amplification, 1, VarIndex)
  if (ierr .ne. 0) then
     call CCTK_WARN(0, "Failed to integrate Spritz_turbulent_energy_amplification_temps")
  end if
  Spritz_turbulent_energy_amplification = Spritz_turbulent_energy_amplification * grid_spacing_product_mhd


end subroutine Spritz_Magnetic_Energy_Global



 /*@@
   @routine    Spritz_MHD_Initialize
   @date       
   @author     Bruno Giacomazzo
   @desc 

   Make sure all the MHD analysis quantities are initialized to zero
   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Spritz_MHD_Initialize(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  divB = 0.0d0
  divB2 = 0.0d0
  divA = 0.0d0
  B_norm = 0.0d0
  A_norm = 0.0d0
  divB_normalized = 0.0d0
  divB2_normalized = 0.0d0
  b2small = 0.0d0
  mhd_press_ratio = 0.0d0
  mhd_energy_ratio = 0.0d0
  
  poynting_flux_x = 0.0d0
  poynting_flux_y = 0.0d0
  poynting_flux_z = 0.0d0
  poynting_flux_norm = 0.0d0
  
  poynting_vector_x = 0.d0
  poynting_vector_y = 0.d0
  poynting_vector_z = 0.d0
  poynting_scalar   = 0.d0

  lambda_MRI_x = 0.0d0
  lambda_MRI_y = 0.0d0
  lambda_MRI_z = 0.0d0
  
  !divB = 0.0d0
  !divB2 = 0.0d0
  !divA = 0.0d0
  !B_norm = 0.0d0
  !divB_normalized = 0.0d0
  !divB2_normalized = 0.0d0
  !b2small = 0.0d0
  !mhd_press_ratio = 0.0d0
  !mhd_energy_ratio = 0.0d0

  mhd_c2p_mask = 0.0d0

  if (Spritz_check_BTor.ne.0) then
     Bphi = 0.0d0 
     Btheta = 0.0d0 
     Brho = 0.0d0
     B_Poloidal_norm = 0.0d0
     B_Toroidal_norm = 0.0d0
     Spritz_mean_B_Tor_numerator_temps = 0.0d0
     Spritz_mean_B_Pol_numerator_temps = 0.0d0
     Spritz_mean_B_Tor_denominator_temps = 0.0d0
     Spritz_mean_B_Tor_numerator = 0.0d0
     Spritz_mean_B_Pol_numerator = 0.0d0
     Spritz_mean_B_Tor_denominator = 0.0d0
     Spritz_mean_B_Tor = 0.0d0
     Spritz_mean_B_Pol = 0.0d0
  end if

  if (compute_EMenergy.ne.0) then
     Spritz_EMenergy_temps = 0.0d0
!!$     Spritz_EMenergy_temps_p = 0.0d0 !this causes segfault after initial data if used with more than one refinement level
!!$     Spritz_EMenergy_temps_p_p = 0.0d0 !this causes segfault after initial data if used with more than one refinement level
     Spritz_total_magnetic_energy = 0.0d0
     grid_spacing_product_mhd = 0.0d0
     Spritz_kinetic_energy_temps = 0.0d0
     Spritz_total_kinetic_energy = 0.0d0
     Spritz_mean_B_numerator_temps = 0.0d0
     Spritz_mean_B_denominator_temps = 0.0d0
     Spritz_mean_B_numerator = 0.0d0
     Spritz_mean_B_denominator = 0.0d0
     Spritz_mean_B = 0.0d0
     Spritz_turbulent_energy_temps = 0.0d0
     Spritz_turbulent_energy = 0.0d0
     Spritz_turbulent_energy_amplification_temps = 0.0d0
     Spritz_turbulent_energy_amplification = 0.0d0
  end if


end subroutine Spritz_MHD_Initialize
