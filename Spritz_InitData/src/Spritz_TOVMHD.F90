 /*@@
   @file      Spritz_TOVMHD.F90
   @date
   @author    Bruno Giacomazzo
   @desc
   Add a poloidal magnetic field to a TOV.
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


 /*@@
   @file      Spritz_TOVMHD.F90
   @date
   @author    Bruno Giacomazzo
   @desc
   Add a poloidal magnetic field to a TOV.
   @enddesc
@@*/

subroutine Spritz_TOVMHD(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT::nx,ny,nz,i,j,k
  CCTK_REAL::Pcut,Aphi_local,idx,idy,idz,cylrad2,x_local,y_local,z_local,det,maxima_x_l,maxima_y_l,maxima_z_l
  CCTK_INT :: k_lowerbound, j_lowerbound, i_lowerbound
  CCTK_INT :: iplus, iminus, jplus, jminus, kplus, kminus
  CCTK_INT :: iiplus, iiminus, jjplus, jjminus, kkplus, kkminus

  call CCTK_INFO("Adding a poloidal magnetic field to the TOV")

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  idx = 1.d0 / CCTK_DELTA_SPACE(1)
  idy = 1.d0 / CCTK_DELTA_SPACE(2)
  idz = 1.d0 / CCTK_DELTA_SPACE(3)

  !initialize the vector potential and magnetic field to zero
  Bx = 0.0
  Bx_p = 0.0
  dBx = 0.0
  dBx_p = 0.0
  By = 0.0
  By_p = 0.0
  dBy = 0.0
  dBy_p = 0.0
  Bz = 0.0
  Bz_p = 0.0
  dBz = 0.0
  dBz_p = 0.0

  Avec_x = 0.0
  Avec_x_p = 0.0
  Avec_y = 0.0
  Avec_y_p = 0.0
  Avec_z = 0.0
  Avec_z_p = 0.0

  psi_mhd = 0.0
  psi_mhd_p = 0.0

  if (tovmhd_timelevels==3) then
     Bx_p_p = 0.0
     dBx_p_p = 0.0
     By_p_p = 0.0
     dBy_p_p = 0.0
     Bz_p_p = 0.0
     dBz_p_p = 0.0
     Avec_x_p_p = 0.0
     Avec_y_p_p = 0.0
     Avec_z_p_p = 0.0
     psi_mhd_p_p = 0.0
  end if

  !compute the maximum of the pressure
  !call CCTK_ReductionArrayHandle (max_handle, 'maximum')
  !call CCTK_VarIndex(VarIndex, "Spritz::press")
  !call CCTK_Reduce(ierr, cctkGH, -1, max_handle,1, CCTK_VARIABLE_REAL, max_press, 1, VarIndex)
  !if (ierr .ne. 0) then
  !!  call CCTK_WARN(0, "Failed to compute the maximum of the pressure")
  !end if

  write(*,*) 'max_press=',max_press,"1/dx=",idx,"1/dy=",idy,"1/dz=",idz

  !cutoff on the pressure
  Pcut = press_cut * max_press

  !Done in order to avoid an intel compiler bug when using optimization
  !(it was causing segfault because of prefetching of
  ! maxima_x and maxima_y inside the following do loop)
  if(compute_DCandMaxima .ne. 0) then
     maxima_x_l = maxima_x_plus
     maxima_y_l = maxima_y_plus
     maxima_z_l = maxima_z_plus
  else
     maxima_x_l = 0.0d0
     maxima_y_l = 0.0d0
     maxima_z_l = 0.0d0
  end if

  do i=1,nx
     do j=1,ny
        do k=1,nz


           if(compute_DCandMaxima .ne. 0) then
              !cylradius is set to be the distance from the position of
              ! the maximum in the rest-mass density in the case the star
              ! is not centered in x=y=z=0
              cylrad2 = ((x(i,j,k)-maxima_x_l)**2+(y(i,j,k)-maxima_y_l)**2) !cylindrical radius squared
!!$              phi = atan2((y(i,j,k)-maxima_y_l),(x(i,j,k)-maxima_x_l))
              x_local = x(i,j,k)-maxima_x_l
              y_local = y(i,j,k)-maxima_y_l
              z_local = z(i,j,k)-maxima_z_l
           else
              cylrad2 = (x(i,j,k)**2+y(i,j,k)**2) !cylindrical radius squared
!!$              phi = atan2(y(i,j,k),x(i,j,k))
              x_local = x(i,j,k)
              y_local = y(i,j,k)
              z_local = z(i,j,k)
           end if

!!$           if (cylrad<=1.0d-10) phi=0.0d0 !to avoid NaN coming from atan

           if (CCTK_Equals(Avec_initial,"Dipole_Shibata").eq.1) then

              !Vector potential. NOTE: Aphi= A_\phi/cylrad2

              Aphi_local = (Ab1*cylrad0 /sqrt(2.0d0)) / (cylrad2 + z_local**2.0d0 + 0.5d0*cylrad0**2.0d0)**(3.0d0/2.0d0)
              Avec_x(i,j,k) = -y_local*Aphi_local
              Avec_y(i,j,k) = x_local*Aphi_local
              Avec_z(i,j,k) = 0.0d0

           elseif (CCTK_Equals(Avec_initial,"homogeneous").eq.1) then

              !Vector potential. NOTE: Aphi= A_\phi/cylrad2

              call SpatialDet(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det)

              Aphi_local = Ab1 * sqrt(det) / 3.0d0
              !Aphi_local = Ab1 * sqrt(cylrad2) / 3.0d0

              Avec_x(i,j,k) = -y_local*Aphi_local
              Avec_y(i,j,k) = x_local*Aphi_local
              Avec_z(i,j,k) = 0.0d0

           else

              if (rho(i,j,k)>Spritz_rho_min) then

                 !Vector potential. NOTE: this is A_\phi/cylrad2

                 if (CCTK_Equals(Avec_initial,"Bruno").eq.1) then

                    Aphi_local = Ab1*cylrad2*max((press(i,j,k)-Pcut),0.0d0)

                 elseif (CCTK_Equals(Avec_initial,"Standard").eq.1) then

                    Aphi_local = Ab1*max((press(i,j,k)-Pcut),0.0d0)

                 elseif (CCTK_Equals(Avec_initial,"Modified").eq.1) then

                    Aphi_local = Ab1*(max((press(i,j,k)-Pcut),0.0d0))**Avec_kappa

                 else
                    call CCTK_WARN(0,"Parameter Avec_initial is not valid!")
                 end if


!!$                 Avec_x(i,j,k) = -1.0d0*cylrad*Aphi_local*sin(phi)!*sign(1.0d0,x(i,j,k))
!!$                 Avec_y(i,j,k) = cylrad*Aphi_local*cos(phi)!*sign(1.0d0,x(i,j,k))
!!$                 Avec_z(i,j,k) = 0.0d0

                 Avec_x(i,j,k) = -y_local*Aphi_local
                 Avec_y(i,j,k) = x_local*Aphi_local
                 Avec_z(i,j,k) = 0.0d0

                 !if (abs(y(i,j,k))<=1.0d-10) Avec_x(i,j,k)=0.0d0
                 !if (abs(x(i,j,k))<=1.0d-10) Avec_y(i,j,k)=0.0d0

              else
                 !vector potential set to zero in the atmosphere
                 Avec_x(i,j,k) = 0.0d0
                 Avec_y(i,j,k) = 0.0d0
                 Avec_z(i,j,k) = 0.0d0
              end if

           end if

           Avec_x_p(i,j,k) = Avec_x(i,j,k)
           Avec_y_p(i,j,k) = Avec_y(i,j,k)
           Avec_z_p(i,j,k) = Avec_z(i,j,k)
 
           if (tovmhd_timelevels==3) then 
                Avec_x_p_p(i,j,k) = Avec_x(i,j,k)
                Avec_y_p_p(i,j,k) = Avec_y(i,j,k)
                Avec_z_p_p(i,j,k) = Avec_z(i,j,k)
           end if

        end do
     end do
  end do

  if(CCTK_Equals(curlA_method,"Staggered").eq.1) then
        k_lowerbound = 2
        j_lowerbound = 2
        i_lowerbound = 2
  else
        k_lowerbound = 1
        j_lowerbound = 1
        i_lowerbound = 1
  end if


  do i=i_lowerbound,nx
     do j=j_lowerbound,ny
        do k=k_lowerbound,nz

           call SpatialDet(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det)

           !B=curlA/sqrt(det)

           if (CCTK_Equals(curlA_method,"Centered").eq.1) then
 
              iplus  = 1
              iminus = 1
              jplus  = 1
              jminus = 1
              kplus  = 1
              kminus = 1

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

              Bx(i,j,k) = 0.5d0*(idy*(Avec_z(i,j+jplus,k)-Avec_z(i,j-jminus,k)) - idz*(Avec_y(i,j,k+kplus)-Avec_y(i,j,k-kminus)))
              By(i,j,k) = 0.5d0*(idz*(Avec_x(i,j,k+kplus)-Avec_x(i,j,k-kminus)) - idx*(Avec_z(i+iplus,j,k)-Avec_z(i-iminus,j,k)))
              Bz(i,j,k) = 0.5d0*(idx*(Avec_y(i+iplus,j,k)-Avec_y(i-iminus,j,k)) - idy*(Avec_x(i,j+jplus,k)-Avec_x(i,j-jminus,k)))

           elseif (CCTK_Equals(curlA_method,"SemiAnalytic").eq.1) then

              if (CCTK_Equals(Avec_initial,"Bruno").eq.0) then
                 call CCTK_WARN(0,"curlA_method=SemiAnalytic can be used only with Avec_initial=Bruno")
              end if


              !Semi-Analytic expressions
              Bx(i,j,k) = -(x(i,j,k)-maxima_x_l)*Ab1*((x(i,j,k)-maxima_x_l)**2+(y(i,j,k)-maxima_y_l)**2)* &
                   (max((press(i,j,k+1)-Pcut),0.0d0)-max((press(i,j,k-1)-Pcut),0.0d0))*idz*0.5d0

              By(i,j,k) = -(y(i,j,k)-maxima_y_l)*Ab1*((x(i,j,k)-maxima_x_l)**2+(y(i,j,k)-maxima_y_l)**2)* &
                   (max((press(i,j,k+1)-Pcut),0.0d0)-max((press(i,j,k-1)-Pcut),0.0d0))*idz*0.5d0

              Bz(i,j,k) = Ab1*((x(i,j,k)-maxima_x_l)**2+(y(i,j,k)-maxima_y_l)**2)* &
                   (4.0d0*max((press(i,j,k)-Pcut),0.0d0) + &
                   (y(i,j,k)-maxima_y_l)*(max((press(i,j+1,k)-Pcut),0.0d0)-max((press(i,j-1,k)-Pcut),0.0d0))*idy*0.5d0 + &
                   (x(i,j,k)-maxima_x_l)*(max((press(i+1,j,k)-Pcut),0.0d0)-max((press(i-1,j,k)-Pcut),0.0d0))*idx*0.5d0)

           elseif (CCTK_Equals(curlA_method,"Staggered").eq.1) then

              if (i .eq. 1) then
                 iiminus = 0
                 iiplus = 1
              else
                 iiminus = 1
                 iiplus = 0
              end if

              if (j .eq. 1) then
                 jjminus = 0
                 jjplus = 1
              else
                 jjminus = 1
                 jjplus = 0
              end if

              if (k .eq. 1) then
                 kkminus = 0
                 kkplus = 1
              else
                 kkminus = 1
                 kkplus = 0
              end if

              !Bx = dy Az - dz Ay
              !dBx is 0.5*(curl(A) at i+1/2 + curl(A) at i-1/2)
              Bx(i,j,k) = 0.5d0 * ((idy*(Avec_z(i+iiplus,j+jjplus,k+kkplus)-Avec_z(i+iiplus,j-jjminus,k+kkplus)) - &
                                    idz*(Avec_y(i+iiplus,j+jjplus,k+kkplus)-Avec_y(i+iiplus,j+jjplus,k-kkminus))) + &
                                   (idy*(Avec_z(i-iiminus,j+jjplus,k+kkplus)-Avec_z(i-iiminus,j-jjminus,k+kkplus)) - &
                                    idz*(Avec_y(i-iiminus,j+jjplus,k+kkplus)-Avec_y(i-iiminus,j+jjplus,k-kkminus))) )
              !By = -dx Az + dz Ax
              !dBy is 0.5*(curl(A) at j+1/2 + curl(A) at j-1/2)
              By(i,j,k) = 0.5d0 * ((idz*(Avec_x(i+iiplus,j+jjplus,k+kkplus)-Avec_x(i+iiplus,j+jjplus,k-kkminus)) - &
                                    idx*(Avec_z(i+iiplus,j+jjplus,k+kkplus)-Avec_z(i-iiminus,j+jjplus,k+kkplus))) + &
                                   (idz*(Avec_x(i+iiplus,j-jjminus,k+kkplus)-Avec_x(i+iiplus,j-jjminus,k-kkminus)) - &
                                    idx*(Avec_z(i+iiplus,j-jjminus,k+kkplus)-Avec_z(i-iiminus,j-jjminus,k+kkplus))) )
              !Bz= dx Ay - dy Ax
              !dBz is 0.5*(curl(A) at k+1/2 + curl(A) at k-1/2)
              Bz(i,j,k) = 0.5d0 * ((idx*(Avec_y(i+iiplus,j+jjplus,k+kkplus)-Avec_y(i-iiminus,j+jjplus,k+kkplus)) - &
                                    idy*(Avec_x(i+iiplus,j+jjplus,k+kkplus)-Avec_x(i+iiplus,j-jjminus,k+kkplus))) + &
                                   (idx*(Avec_y(i+iiplus,j+jjplus,k-kkminus)-Avec_y(i-iiminus,j+jjplus,k-kkminus)) - &
                                    idy*(Avec_x(i+iiplus,j+jjplus,k-kkminus)-Avec_x(i+iiplus,j-jjminus,k-kkminus))) )

           else
              call CCTK_WARN(0,"Not valid parameter for curlA_method. I don't know how to build the magnetic field!")
           end if

           Bx(i,j,k) = Bx(i,j,k)/sqrt(det)
           By(i,j,k) = By(i,j,k)/sqrt(det)
           Bz(i,j,k) = Bz(i,j,k)/sqrt(det)

           if (CCTK_Equals(Spritz_eos_type,"General").eq.1) then

              call Prim2ConGen_MHD(Spritz_eos_handle,gxx(i,j,k),gxy(i,j,k),gxz(i,j,k), &
                   gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det, &
                   dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),tau(i,j,k), &
                   rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k), &
                   eps(i,j,k),press(i,j,k),w_lorentz(i,j,k), &
                   Bx(i,j,k),By(i,j,k),Bz(i,j,k),dBx(i,j,k),dBy(i,j,k),dBz(i,j,k),alp(i,j,k), &
                   betax(i,j,k),betay(i,j,k),betaz(i,j,k))

           else if (CCTK_Equals(Spritz_eos_type,"Polytype").eq.1) then

              call Prim2ConPoly_MHD(Spritz_eos_handle,gxx(i,j,k),gxy(i,j,k),gxz(i,j,k), &
                   gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det, &
                   dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),tau(i,j,k), &
                   rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k), &
                   eps(i,j,k),press(i,j,k),w_lorentz(i,j,k), &
                   Bx(i,j,k),By(i,j,k),Bz(i,j,k),dBx(i,j,k),dBy(i,j,k),dBz(i,j,k),alp(i,j,k), &
                   betax(i,j,k),betay(i,j,k),betaz(i,j,k))

           else
              call CCTK_WARN(0, "Unrecognized EoS type");
           end if

           !Set previous timesteps
           dens_p(i,j,k) = dens(i,j,k)
           scon_p(i,j,k,1)   = sx(i,j,k)
           scon_p(i,j,k,2)   = sy(i,j,k)
           scon_p(i,j,k,3)   = sz(i,j,k)
           tau_p(i,j,k)  = tau(i,j,k)

           Bx_p(i,j,k)  = Bx(i,j,k)
           By_p(i,j,k)  = By(i,j,k)
           Bz_p(i,j,k)  = Bz(i,j,k)
           dBx_p(i,j,k) = dBx(i,j,k)
           dBy_p(i,j,k) = dBy(i,j,k)
           dBz_p(i,j,k) = dBz(i,j,k)


           if (tovmhd_timelevels==3) then
              dens_p_p(i,j,k) = dens(i,j,k)
              scon_p_p(i,j,k,1)   = sx(i,j,k)
              scon_p_p(i,j,k,2)   = sy(i,j,k)
              scon_p_p(i,j,k,3)   = sz(i,j,k)
              tau_p_p(i,j,k)  = tau(i,j,k)

              Bx_p_p(i,j,k)  = Bx(i,j,k)
              By_p_p(i,j,k)  = By(i,j,k)
              Bz_p_p(i,j,k)  = Bz(i,j,k)
              dBx_p_p(i,j,k) = dBx(i,j,k)
              dBy_p_p(i,j,k) = dBy(i,j,k)
              dBz_p_p(i,j,k) = dBz(i,j,k)
           end if


        end do
     end do
  end do

!!$  call CCTK_WARN(1,"Now should call Output used for debug")

end subroutine Spritz_TOVMHD


!!$subroutine Spritz_TOVMHD_output(CCTK_ARGUMENTS)
!!$
!!$  implicit none
!!$
!!$  DECLARE_CCTK_ARGUMENTS
!!$  DECLARE_CCTK_FUNCTIONS
!!$
!!$  CCTK_INT:: ierr
!!$
!!$  call CCTK_WARN(1,"Output used for debug")
!!$
!!$  call CCTK_OutputVarAsByMethod(ierr, cctkGH, "Spritz::Bz", "IOASCII_2D", "BzInit")
!!$
!!$  if (ierr .ne. 0) call CCTK_WARN(0,"Method not recognized!")
!!$
!!$  call CCTK_OutputVarAsByMethod(ierr, cctkGH, "Spritz::dBz", "IOASCII_2D", "dBzInit")
!!$
!!$  if (ierr .ne. 0) call CCTK_WARN(0,"Method not recognized!")
!!$
!!$  call CCTK_OutputVarAsByMethod(ierr, cctkGH, "Spritz::press", "IOASCII_2D", "pressInit")
!!$
!!$  if (ierr .ne. 0) call CCTK_WARN(0,"Method not recognized!")
!!$
!!$
!!$end subroutine Spritz_TOVMHD_output


 /*@@
   @file      Spritz_TOVMHD_maxima
   @date
   @author    Bruno Giacomazzo
   @desc
   Call the Spritz routine Spritz_FindMaxima to locate the maximum of the rest-mass density
   @enddesc
 @@*/


subroutine Spritz_TOVMHD_maxima(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  !DECLARE_CCTK_FUNCTIONS

  if(compute_DCandMaxima .ne. 0) then

    call CCTK_INFO("About to call FindMaxima")

    Spritz_maxima_iteration = -1
    maximum_density_plus = -1.d0
    maximum_density_minus = -1.d0

    CALL Spritz_FindMaxima_Guts(                                      &
         & cctkGH, cctk_lsh, cctk_iteration, Spritz_maxima_iteration, &
         & maximum_density_plus, maximum_density_minus,               &
         & maxima_x_plus,  maxima_y_plus,  maxima_z_plus,             &
         & maxima_i_plus,  maxima_j_plus,  maxima_k_plus,             &
         & maxima_x_minus, maxima_y_minus, maxima_z_minus,            &
         & maxima_i_minus, maxima_j_minus, maxima_k_minus,            &
         & rho, x, y, z )

    ! restore the initial values for use in Spritz
    Spritz_maxima_iteration = -1
    maximum_density_plus = -1.d0
    maximum_density_minus = -1.d0

    call CCTK_INFO("After call to FindMaxima")

  end if

end subroutine Spritz_TOVMHD_maxima
