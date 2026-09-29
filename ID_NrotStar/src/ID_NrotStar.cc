/*@@
  @file      ID_NrotStar.cc
  @date
  @author    Lorenzo Sala, modifying WhiskyMHD_MagStar_InitData.C by Bruno Giacomazzo and Riccardo Ciolfi
  @desc
  Reads initial data for a rotating non-magnetized Neutron Star built with Lorene/nrotstar.
  Allows user to set an arbitrary EoS, given as path to CompOSE tables.
  @enddesc
  @@*/

#include <cassert>
#include <cstdio>
#include <vector>
#include <ios>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

// headers C
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fstream>

// headers Lorene
#include <star_rot.h>
#include <scalar.h>
#include <tensor.h>
#include <eos.h>
#include <nbr_spx.h>
#include <unites.h>
#include <metric.h>

// headers Cactus
#include "cctk_Functions.h"

using namespace std;
// define namespace here for old versions of Lorene that don't do so
namespace Lorene {}
using namespace Lorene;

#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velx_p (&vel_p[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely_p (&vel_p[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz_p (&vel_p[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velx_p_p (&vel_p_p[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely_p_p (&vel_p_p[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz_p_p (&vel_p_p[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
// #define Bx (&Bvec[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
// #define By (&Bvec[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
// #define Bz (&Bvec[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

extern "C"
void ID_NrotStar_initialise (CCTK_ARGUMENTS)
{
  using namespace Unites_mag;
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO ("Setting up LORENE ID_NrotStar initial data");
  FILE* nrotstar_data_file = fopen(filename, "r") ;

  int const local_grid_size = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
  double mm0;

  #define CUto10KM  0.147664
  #define reciprocCUto10KM  (1.0 / CUto10KM)

    try {
    Mg3d spectral_grid(nrotstar_data_file) ;
    Map_et mapping(spectral_grid, nrotstar_data_file) ;

    // read tabulated EoS from file eosname, provided by parfile
    // Eos* p_eos;
    // ifstream eosparfile(eosname);
    // p_eos = Eos::eos_from_file(eosparfile) ;
    // eosparfile.close() ;
    Eos* p_eos = Eos::eos_from_file(nrotstar_data_file) ;

    cout << *p_eos;
    CCTK_INFO("Calling Star_rot constructor... ");
    Star_rot star(mapping, *p_eos, nrotstar_data_file);
    CCTK_VINFO("Done constructor on reflevel %d.",int(log10(double(cctk_levfac[0]))/log10(2.0)));
    star.equation_of_state() ;
    star.update_metric() ;
    star.extrinsic_curvature() ;
    star.hydro_euler() ;

    int nx = cctk_lsh[0] ;
    int ny = cctk_lsh[1] ;
    int nz = cctk_lsh[2] ;

    Scalar sp_lapse = star.get_nn() ;//the lapse
    Vector sp_shift = star.get_beta();
    sp_shift.change_triad(mapping.get_bvect_cart()) ;
    Scalar sp_shift_x = sp_shift(1) ;
    Scalar sp_shift_y = sp_shift(2) ;
    Scalar sp_shift_z = sp_shift(3) ;

    Tensor sp_gamma = star.get_gamma().cov(); //covariant form
    sp_gamma.change_triad(mapping.get_bvect_cart()) ;
    Scalar sp_gam_xx = sp_gamma(1,1) ;
    Scalar sp_gam_xy = sp_gamma(1,2) ;
    Scalar sp_gam_xz = sp_gamma(1,3) ;
    Scalar sp_gam_yy = sp_gamma(2,2) ;
    Scalar sp_gam_yz = sp_gamma(2,3) ;
    Scalar sp_gam_zz = sp_gamma(3,3) ;

    //the extrinsic curvature tkij=1/B^2 kij
    Tensor sp_kij = star.get_bbb()*star.get_bbb()*star.get_tkij();
    sp_kij.change_triad(mapping.get_bvect_cart()) ;
    /*
    for (int i=1; i<=3; i++)
      for (int j=1; j<=i; j++) {
        sp_kij.set(i,j).dec_dzpuis(2) ;
    }*/

    Scalar sp_kij_xx = sp_kij(1,1) ; //covariant form...
    Scalar sp_kij_xy = sp_kij(1,2) ;
    Scalar sp_kij_xz = sp_kij(1,3) ;
    Scalar sp_kij_yy = sp_kij(2,2) ;
    Scalar sp_kij_yz = sp_kij(2,3) ;
    Scalar sp_kij_zz = sp_kij(3,3) ;

    //Cmp sp_density = poly_eos->get_m_0() * star.get_nbar() ; old presciption does not work with PP EOS
    mm0 = 1.000 ; // baryon mass set by default to 1, has already been converted.
    Scalar sp_density = star.get_nbar() * mm0 ;

    Vector sp_v = star.get_u_euler(); //contravariant representation!!
    sp_v.change_triad(mapping.get_bvect_cart()) ;
    Scalar sp_v_x = sp_v(1);
    Scalar sp_v_y = sp_v(2);
    Scalar sp_v_z = sp_v(3);

    Scalar sp_eps = star.get_ener() - sp_density ; // not epsilon!!!
    // but rho*epsilon ...


    //assign values to cactus variables on a cartesian grid
    for (int k=0; k<nz; k++) {
      for (int j=0; j<ny; j++) {
        for (int i=0; i<nx; i++) {
  	int i3D=CCTK_GFINDEX3D(cctkGH, i, j, k);

  	// convert position from cactus units to Lorene units (10Km)
  	double x_fd = x[i3D]*CUto10KM;
  	double y_fd = y[i3D]*CUto10KM;
  	double z_fd = z[i3D]*CUto10KM;

  	double r_sp = sqrt(x_fd*x_fd + y_fd*y_fd + z_fd*z_fd) ;
    double rcil = sqrt(x_fd*x_fd + y_fd*y_fd) ;
  	double theta_sp;
  	double phi_sp;

    if (z_fd>0.0) theta_sp = atan(rcil/z_fd);
    else if (z_fd<0.0) theta_sp = 2*acos(0.0)-atan(-rcil/z_fd);  // pi - atan
    else if (z_fd==0.0) theta_sp = acos(0.0);  //this goes from -pi/2 to +pi/2

    if (x_fd>0.0) phi_sp = atan(y_fd/x_fd);  //this goes from -pi/2 to +pi/2
    else if (x_fd<0.0) {
  		if (y_fd<0.0) phi_sp = atan(y_fd/x_fd)-2*acos(0.0); // phi-pi
  		else if (y_fd>=0.0) phi_sp = atan(y_fd/x_fd)+2*acos(0.0); // phi+pi
  	}
    else if (x_fd==0.0) {
        if (y_fd<0.0) phi_sp = -acos(0.0); // -pi/2
        else if (y_fd>0.0) phi_sp = acos(0.0); // pi/2
  	}

    if (phi_sp<0.0) phi_sp = phi_sp+4*acos(0.0); // add 2*pi=4*pi/2
  	if (r_sp<=1.0e-15) theta_sp=0;
  	if ((x_fd*x_fd+y_fd*y_fd)<=1.0e-15) phi_sp=0;

  	rho[i3D]   = sp_density.val_point(r_sp, theta_sp, phi_sp) ;
  	if(set_Lorene_eps) {
        // double totalenergy = sp_eps.val_point(r_sp, theta_sp, phi_sp)/rho[i3D];
        // if(isnan(totalenergy)) eps[i3D] = 0.0;
        eps[i3D]   = sp_eps.val_point(r_sp, theta_sp, phi_sp)/rho[i3D];
        // printf("setting %f\n", eps[i3D]);
    }
  	// v^j
    if ( (x_fd*x_fd+y_fd*y_fd) > 1.0e-15 ) {
        velx[i3D]  = sp_v_x.val_point(r_sp, theta_sp, phi_sp);
        vely[i3D]  = sp_v_y.val_point(r_sp, theta_sp, phi_sp);
        velz[i3D]  = sp_v_z.val_point(r_sp, theta_sp, phi_sp);
    }
    else {
        velx[i3D]  = 0.0;
        vely[i3D]  = 0.0;
        velz[i3D]  = 0.0;
    }

  	// g_ij
  	gxx[i3D]   = sp_gam_xx.val_point(r_sp, theta_sp, phi_sp);
  	gxy[i3D]   = sp_gam_xy.val_point(r_sp, theta_sp, phi_sp);
  	gxz[i3D]   = sp_gam_xz.val_point(r_sp, theta_sp, phi_sp);
  	gyy[i3D]   = sp_gam_yy.val_point(r_sp, theta_sp, phi_sp);
  	gyz[i3D]   = sp_gam_yz.val_point(r_sp, theta_sp, phi_sp);
  	gzz[i3D]   = sp_gam_zz.val_point(r_sp, theta_sp, phi_sp);
  	// the lapse
  	alp[i3D]   = sp_lapse.val_point(r_sp, theta_sp, phi_sp);
  	// and the shift beta^i
  	betax[i3D] = sp_shift_x.val_point(r_sp, theta_sp, phi_sp);
  	betay[i3D] = sp_shift_y.val_point(r_sp, theta_sp, phi_sp);
  	betaz[i3D] = sp_shift_z.val_point(r_sp, theta_sp, phi_sp);
  	// the extrinsic curvature K_ij
  	kxx[i3D]  = sp_kij_xx.val_point(r_sp, theta_sp, phi_sp);
  	kxy[i3D]  = sp_kij_xy.val_point(r_sp, theta_sp, phi_sp);
  	kxz[i3D]  = sp_kij_xz.val_point(r_sp, theta_sp, phi_sp);
  	kyy[i3D]  = sp_kij_yy.val_point(r_sp, theta_sp, phi_sp);
  	kyz[i3D]  = sp_kij_yz.val_point(r_sp, theta_sp, phi_sp);
  	kzz[i3D]  = sp_kij_zz.val_point(r_sp, theta_sp, phi_sp);

  	// Convert the various quantities to Cactus units
    rho[i3D]  *= 2.6873277924735397e-4; // LSMOD changed from 2.687182675e-4; //in Lorene rho_unit=1.66e17 (Kg/m^3), so this is 1.66E17/6.17714E20*10
    // if(rho[i3D] > 1.5E-3) cout <<  rho[i3D] << "\t" << eps[i3D] << endl;
  	// eps doesn't need to be converted, same units in Cactus and Lorene
  	// [K_ij] = 1/[L]
  	kxx[i3D]  *= CUto10KM;
  	kxy[i3D]  *= CUto10KM;
  	kxz[i3D]  *= CUto10KM;
  	kyy[i3D]  *= CUto10KM;
  	kyz[i3D]  *= CUto10KM;
  	kzz[i3D]  *= CUto10KM;
  	// v^j don't need to be converted, they are already in c unit
  	// (c is the speed of light)

        double debcheck = rho[i3D]*eps[i3D]*velx[i3D]*vely[i3D]*velz[i3D]* \
                          gxx[i3D]*gxy[i3D]*gxz[i3D]*gyy[i3D]*gyz[i3D]*gzz[i3D]* \
                          alp[i3D]*betax[i3D]*betay[i3D]*betaz[i3D] ;
        if  (isnan(debcheck)){
                             CCTK_WARN (CCTK_WARN_ABORT, "NaN(s) found reading ID.");
                             }

            }
        }
    }


    CCTK_INFO ("Setting time derivatives of lapse and shift");
    // These initial data assume stationarity
    if (CCTK_EQUALS (initial_dtlapse, "ID_NrotStar")) {
        #pragma omp parallel for
        for (int i=0; i<local_grid_size; ++i) {
            dtalp[i] = 0.0;
        }
    } else if (CCTK_EQUALS (initial_dtlapse, "none")) {
      // do nothing
    } else {
      CCTK_WARN (CCTK_WARN_ABORT, "internal error");
    }

    // CCTK_INFO ("Setting magnetic field to zero.");
    // #pragma omp parallel for
    // for (int i=0; i<local_grid_size; ++i) {
    //     Bx[i] = 0.;
    //   	By[i] = 0.;
    //   	Bz[i] = 0.;
    // }

    if (CCTK_EQUALS (initial_dtshift, "ID_NrotStar")) {
        #pragma omp parallel for
        for (int i=0; i<local_grid_size; ++i) {
            dtbetax[i] = 0.0;
            dtbetay[i] = 0.0;
            dtbetaz[i] = 0.0;
        }
    } else if (CCTK_EQUALS (initial_dtshift, "none")) {
      // do nothing
    } else {
      CCTK_WARN (CCTK_WARN_ABORT, "internal error");
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::lapse") > 1) {
        memcpy(alp_p, alp, local_grid_size*sizeof(CCTK_REAL));
        if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::lapse") > 2){
	        memcpy(alp_p_p, alp, local_grid_size*sizeof(CCTK_REAL));
	    }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::shift") > 1) {
        memcpy(betax_p, betax, local_grid_size*sizeof(CCTK_REAL));
        memcpy(betay_p, betay, local_grid_size*sizeof(CCTK_REAL));
        memcpy(betaz_p, betaz, local_grid_size*sizeof(CCTK_REAL));

        if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::shift") > 2) {
	        memcpy(betax_p_p, betax, local_grid_size*sizeof(CCTK_REAL));
	        memcpy(betay_p_p, betay, local_grid_size*sizeof(CCTK_REAL));
	        memcpy(betaz_p_p, betaz, local_grid_size*sizeof(CCTK_REAL));
	    }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::metric") > 1) {
        memcpy(gxx_p, gxx, local_grid_size*sizeof(CCTK_REAL));
        memcpy(gxy_p, gxy, local_grid_size*sizeof(CCTK_REAL));
        memcpy(gxz_p, gxz, local_grid_size*sizeof(CCTK_REAL));
        memcpy(gyy_p, gyy, local_grid_size*sizeof(CCTK_REAL));
        memcpy(gyz_p, gyz, local_grid_size*sizeof(CCTK_REAL));
        memcpy(gzz_p, gzz, local_grid_size*sizeof(CCTK_REAL));

        if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::metric") > 2) {
	        memcpy(gxx_p_p, gxx, local_grid_size*sizeof(CCTK_REAL));
	        memcpy(gxy_p_p, gxy, local_grid_size*sizeof(CCTK_REAL));
	        memcpy(gxz_p_p, gxz, local_grid_size*sizeof(CCTK_REAL));
	        memcpy(gyy_p_p, gyy, local_grid_size*sizeof(CCTK_REAL));
	        memcpy(gyz_p_p, gyz, local_grid_size*sizeof(CCTK_REAL));
	        memcpy(gzz_p_p, gzz, local_grid_size*sizeof(CCTK_REAL));
	    }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::curv") > 1) {
        memcpy(kxx_p, kxx, local_grid_size*sizeof(CCTK_REAL));
        memcpy(kxy_p, kxy, local_grid_size*sizeof(CCTK_REAL));
        memcpy(kxz_p, kxz, local_grid_size*sizeof(CCTK_REAL));
        memcpy(kyy_p, kyy, local_grid_size*sizeof(CCTK_REAL));
        memcpy(kyz_p, kyz, local_grid_size*sizeof(CCTK_REAL));
        memcpy(kzz_p, kzz, local_grid_size*sizeof(CCTK_REAL));

        if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::curv") > 2) {
        	memcpy(kxx_p_p, kxx, local_grid_size*sizeof(CCTK_REAL));
        	memcpy(kxy_p_p, kxy, local_grid_size*sizeof(CCTK_REAL));
        	memcpy(kxz_p_p, kxz, local_grid_size*sizeof(CCTK_REAL));
        	memcpy(kyy_p_p, kyy, local_grid_size*sizeof(CCTK_REAL));
        	memcpy(kyz_p_p, kyz, local_grid_size*sizeof(CCTK_REAL));
        	memcpy(kzz_p_p, kzz, local_grid_size*sizeof(CCTK_REAL));
	    }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::rho") > 1) {
        memcpy(rho_p, rho, local_grid_size*sizeof(CCTK_REAL));
        memcpy(velx_p, velx, local_grid_size*sizeof(CCTK_REAL));
        memcpy(vely_p, vely, local_grid_size*sizeof(CCTK_REAL));
        memcpy(velz_p, velz, local_grid_size*sizeof(CCTK_REAL));
        memcpy(press_p, press, local_grid_size*sizeof(CCTK_REAL));
        memcpy(w_lorentz_p, w_lorentz, local_grid_size*sizeof(CCTK_REAL));
        memcpy(eps_p, eps, local_grid_size*sizeof(CCTK_REAL));

        if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::rho") > 2) {
        	memcpy(rho_p_p, rho, local_grid_size*sizeof(CCTK_REAL));
        	memcpy(velx_p_p, velx, local_grid_size*sizeof(CCTK_REAL));
        	memcpy(vely_p_p, vely, local_grid_size*sizeof(CCTK_REAL));
        	memcpy(velz_p_p, velz, local_grid_size*sizeof(CCTK_REAL));
        	memcpy(press_p_p, press, local_grid_size*sizeof(CCTK_REAL));
        	memcpy(w_lorentz_p_p, w_lorentz, local_grid_size*sizeof(CCTK_REAL));
        	memcpy(eps_p_p, eps, local_grid_size*sizeof(CCTK_REAL));
	    }
    }

  CCTK_INFO ("Done.");
  } catch (ios::failure e) {
    CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Could not read initial data from file '%s': %s", filename, e.what());
  }
}
