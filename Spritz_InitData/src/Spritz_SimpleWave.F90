 /*@@
	@file		Spritz_SimpleWave.F90
	@date
	@author	Beatrice Giudici
	@desc
	Simple Wave Test for convergence
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
	@routine		Spritz_SimpleWave
	@date
	@author	Beatrice Giudici
	@desc
	Initial data for a simple wave test with sin-like function for the velocity
	(see Bernuzzi et al. 2016, and references therein)
	@enddesc
 @@*/

subroutine Spritz_SimpleWave(CCTK_ARGUMENTS)

	implicit none
	
	DECLARE_CCTK_ARGUMENTS
	DECLARE_CCTK_PARAMETERS
	DECLARE_CCTK_FUNCTIONS
	
	CCTK_INT :: nx, ny, nz, i, j, k
	CCTK_REAL :: pi, det, a, x_param, gam, c_0
	
! begin EOS Omni vars
	CCTK_INT  :: n = 1
	CCTK_INT  :: keytemp = 0
	CCTK_INT  :: anyerr = 0
	CCTK_INT  :: keyerr(1) = 0
	CCTK_REAL :: xpress(1) = 0.0d0
	CCTK_REAL :: xeps(1) = 0.0d0
	CCTK_REAL :: xtemp(1) = 0.0d0
	CCTK_REAL :: xye(1) = 0.0d0
	CCTK_REAL :: rf_precision = 1.0d-10
! end EOS Omni vars

	call CCTK_INFO("Setting initial data for a simple wave")
	
	pi = 4.0d0 * atan(1.0d0)
	
	nx = cctk_lsh(1)
	ny = cctk_lsh(2)
	nz = cctk_lsh(3)

	velx(:,:,:) = 0.0d0
   rho = 0.0d0
   eps = 0.0d0
   press = 0.0d0
   w_lorentz = 0.0d0
   simple_tmp = 0.0d0

! Setting the parameters as in Bernuzzi et al.
	a = 0.5d0
	x_param = 0.3d0
	gam = poly_gamma
	
	write(*,*) 'a = ', a, 'X = ', x_param, 'gamma = ', gam

! Compute the velocity	
	do i = 1, nx
		! Condition given by the Heaviside function
		if (abs(x(i,1,1)) < x_param) then
			velx(i,1,1) = a * (sin( (pi / 2.0d0) * ((x(i,1,1)/x_param) - 1.0d0) ))**6
		end if
	end do
	
	! Computing quantities at v=0
	simple_rho_0 = 1.0d0 ! a parameter

	! NOTE: The sound speed in the reference state is given by the exact solution, 
	! see Eq. 5.92 in A.M. Anile 'Relativistic fluids and magneto-fluids' (1990), 
	! where the pressure is given by a polytropic EOS

	c_0 = sqrt(gam * (poly_k*simple_rho_0**gam) * (gam - 1.0d0) / ((gam - 1.0d0) * &
			simple_rho_0 + gam*poly_k*simple_rho_0**gam))
	
	simple_eps_0 = (c_0**2) / (gam * (gam - 1.0d0 - c_0**2))
	
	write(*,*) 'rho_0 = ', simple_rho_0, 'eps_0 = ', simple_eps_0, 'c_0 = ', c_0 

	do j = 1, nx
	! Computing the sound speed according to the value of the Riemann invariant as in 
	! Eq. II.15 of E.P.T. Liang, Astrophys. J. 211, 361 (1977), or Eq. 5.93 in Anile (1990)
		simple_tmp(j,1,1) = ((sqrt(gam - 1.0d0) + c_0) / (sqrt(gam - 1.0d0) - c_0)) * &
	 		 ((1.0d0 + velx(j,1,1)) / (1.0d0 - velx(j,1,1))) ** (sqrt(gam - 1.0d0)/2.0d0)
	 
		c_s(j,1,1) = sqrt(gam - 1.0d0)*(simple_tmp(j,1,1) - 1.0d0) / (simple_tmp(j,1,1) + 1.0d0)
	
	! Computing hydro quantities (N.B. they are normalized over the corresponding quantities in 
	! the reference state) 
		eps(j,1,1) = (c_s(j,1,1)**2.0d0) / (gam * (gam - 1.0d0 - c_s(j,1,1)**2.0d0))
	 
	 	rho(j,1,1) = eps(j,1,1)**(1.0d0 / (gam - 1.0d0))
	 
	 	press(j,1,1) = eps(j,1,1)**(gam / (gam - 1.0d0))
	
	! Computing the Lorentz factor in a flat spacetime 
		w_lorentz(j,1,1) = 1.0d0 / sqrt(1.0d0 - velx(j,1,1)**2)
	end do

	gxx = 1.0d0
	gyy = 1.0d0
	gzz = 1.0d0
	gxy = 0.0d0
	gxz = 0.0d0
	gyz = 0.0d0
	alp = alp_value
	Bx = 0.0d0 
	By = 0.0d0 
	Bz = 0.0d0 
	dBx = Bx
	dBy = By
	dBz = Bz
	Avec_x = 0.0d0
	Avec_y = 0.0d0
	Avec_z = 0.0d0
	if (CCTK_EQUALS(initial_shift,"zero")) then
		betax = 0.0d0
		betay = 0.0d0
		betaz = 0.0d0
	end if
	
	do i = 1, nx
		! atmosphere
		
		if ( (rho(i,1,1) < Spritz_rho_min) .or. (velx(i,1,1) < 0.0d0) ) then
			rho(i,1,1) = rho_abs_min
			eps(i,1,1) = rho_abs_min**(gam - 1.0d0)
			velx(i,1,1) = 0.0d0
			w_lorentz(i,1,1) = 1.0d0
			
			xeps = 1.0d0
			
			call EOS_Omni_press(Spritz_polytrope_handle,keytemp,rf_precision,n,&
				  rho(i,1,1),xeps,xtemp,xye,press(i,1,1),keyerr,anyerr)
				  
		end if
		
		call SpatialDet(gxx(i,1,1),gxy(i,1,1),gxz(i,1,1),&
                gyy(i,1,1),gyz(i,1,1),gzz(i,1,1),det)

! always use polytype for these initial data 
! NOTE: there is no prim2conpolytype without magentic fields  
      call prim2conpolytype_mhd(Spritz_polytrope_handle,gxx(i,1,1),gxy(i,1,1),&
      	  gxz(i,1,1),gyy(i,1,1),gyz(i,1,1),gzz(i,1,1),det,dens(i,1,1),&
      	  sx(i,1,1),sy(i,1,1),sz(i,1,1),tau(i,1,1),rho(i,1,1),&
      	  velx(i,1,1),vely(i,1,1),velz(i,1,1),eps(i,1,1),press(i,1,1),&
      	  w_lorentz(i,1,1),Bx(i,1,1), By(i,1,1), Bz(i,1,1), &
           dBx(i,1,1), dBy(i,1,1), dBz(i,1,1),alp(i,1,1),&
           betax(i,1,1),betay(i,1,1),betaz(i,1,1))


	end do
	
! Planar symmetry
	do j = 1, ny
		do k = 1, nz
			velx(:,j,k) = velx(:,1,1)
			rho(:,j,k) = rho(:,1,1)
			eps(:,j,k) = eps(:,1,1)
			press(:,j,k) = press(:,1,1)
			w_lorentz(:,j,k) = w_lorentz(:,1,1)
			
			sx(:,j,k) = sx(:,1,1)
			dens(:,j,k) = dens(:,1,1)
			tau(:,j,k) = tau(:,1,1)
		end do
	end do
	
	
	Bx = 0.0d0 
	By = 0.0d0 
	Bz = 0.0d0 
	dBx = Bx
	dBy = By
	dBz = Bz
	Avec_x = 0.0d0
	Avec_y = 0.0d0
	Avec_z = 0.0d0
	
	vely(:,:,:) = 0.0d0
	velz(:,:,:) = 0.0d0
	
	densrhs = 0.0d0
	srhs = 0.0d0
	taurhs = 0.0d0
	
	Bxrhs = 0.0d0
   Byrhs = 0.0d0
   Bzrhs = 0.0d0
   
   Avec_x_rhs = 0.0d0
   Avec_y_rhs = 0.0d0
   Avec_z_rhs = 0.0d0
	
	simple_tmp = rho
	
	call CCTK_INFO("Finished initial data")
	
	return

end subroutine Spritz_SimpleWave


subroutine Spritz_SimpleWave_Analysis(CCTK_ARGUMENTS)

	implicit none
	
	DECLARE_CCTK_ARGUMENTS
	DECLARE_CCTK_PARAMETERS
	DECLARE_CCTK_FUNCTIONS
	
	if(velx(cctk_lsh(1),1,1) == 0.0d0) then
		simple_eps_0 = eps(cctk_lsh(1),1,1)
	else
		call cctk_warn(1,"The wave has reached the outer boundary: the computation of simple_eps is now wrong")
	end if
	
	simple_rho = (rho - simple_rho_0)/simple_rho_0
	simple_eps = (eps - simple_eps_0)/simple_eps_0
	
	return

end subroutine Spritz_SimpleWave_Analysis

