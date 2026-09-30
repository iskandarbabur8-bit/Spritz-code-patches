/***********************************************************************************
    Copyright 2006 Scott C. Noble, Charles F. Gammie, Jonathan C. McKinney,
                   and Luca Del Zanna.

                        PVS_GRMHD

    This file was derived from PVS_GRMHD.  The authors of PVS_GRMHD include
    Scott C. Noble, Charles F. Gammie, Jonathan C. McKinney, and Luca Del Zanna.
    PVS_GRMHD is available under the GPL from:
    http://rainman.astro.uiuc.edu/codelib/

    You are morally obligated to cite the following  paper in his/her
    scientific literature that results from use of this file:

    [1] Noble, S. C., Gammie, C. F., McKinney, J. C., \& Del Zanna, L. \ 2006,
        Astrophysical Journal, 641, 626.

    PVS_GRMHD is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    PVS_GRMHD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PVS_GRMHD; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

    If the user has any questions, please direct them to Scott C. Noble at
    scn@astro.rit.edu  .

***********************************************************************************/
/***********************************************************************************

    This file has been adapted from GRHydro_Con2PrimM_pt_EOSOmni.c in order
    to work with Spritz (that will be probably called SpritzGRMHD).

***********************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <complex.h>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Arguments.h>
#include <cctk_Schedule.h>

/* Set this to be 1 if you want debug output */
#define DEBUG_CON2PRIMM (0)  /* CHANGED */


/* Adiabatic index used for the state equation */

#define MAX_NEWT_ITER   (30)       /* Max. # of Newton-Raphson iterations for find_root_2D(); */
#define NEWT_TOL        (1.0e-10)  /* Min. of tolerance allowed for Newton-Raphson iterations */
#define MIN_NEWT_TOL    (1.0e-10)  /* Max. of tolerance allowed for Newton-Raphson iterations */
#define EXTRA_NEWT_ITER (2)       //2
#define NEWT_RELAX      (1.0e-4)   /* Relaxation factor for Newton-Raphson */ // NOT USED NOW
#define FACTOR          (9.0e-1)   /* recurrence factor to find the Newton-Raphson relaxation factor */
#define MIN_RELAX       (1.0e-10)   /* minimum relaxation factor in order to avoid infinite loop or bad convergence */

#define NEWT_TOL2       (1.0e-15)  /* TOL of new 1D^*_{v^2} gnr2 method */
#define MIN_NEWT_TOL2   (1.0e-10)  /* TOL of new 1D^*_{v^2} gnr2 method */

#define W_TOO_BIG       (1.e20)    /* \gamma^2 (\rho_0 + u + p) is assumed
                                      to always be smaller than this.  This
                                      is used to detect solver failures */

#define  eps_dummy      (0.0)      /* Dummy eps variable for EOS_Omni functions */

#define FAIL_VAL        (1.e30)    /* Generic value to which we set variables when a problem arises */

/**************************************************
  The following functions assume a Gamma-law EOS:
***************************************************/

/* Local Globals */
struct LocGlob {
  CCTK_REAL Bsq,bsq,QdotBsq,VdotBsq,Qtsq,Qdotn,D,half_Bsq,tau ;      /*  - DON'T KNOW IF THEY ARE USEFUL*/
} ;

struct eosomnivars {
  CCTK_INT eoshandle,eoskeytemp;
  CCTK_REAL prec,eosprec,eos_y_e[1],eos_temp[1];
  CCTK_INT eoskeyerr[1],eosanyerr[1];
  CCTK_INT eosrl,eosit,eosi,eosj,eosk;         /*  */
  CCTK_REAL eosx,eosy,eosz;
  /* Matched-polytrope calibration constants for extrapolating cold_tabulated
   * EOS below its table rho_min (see Spritz_coldtab_calib in interface.ccl,
   * computed once by Spritz_Rho_Minima_Setup_Final in Spritz_Minima.cc).
   * Used by the atmosphere-reset branch below so it never has to call
   * EOS_Omni during a con2prim failure. */
  CCTK_REAL coldtab_rhomin,coldtab_press_lo,coldtab_eps_lo,coldtab_gamma_lo;
} ;

// Declarations:
static CCTK_REAL vsq_calc(CCTK_REAL W, struct LocGlob *lgp);

static CCTK_INT twod_newton_raphson( CCTK_REAL x[], CCTK_REAL gammaeos, struct LocGlob *lgp,
                                     void (*funcd) (CCTK_REAL [], CCTK_REAL [],
                                                    CCTK_REAL [], CCTK_REAL [][2],
                                                    CCTK_REAL *, CCTK_REAL *, CCTK_REAL, struct LocGlob *) );

static CCTK_INT threed_newton_raphson_omni( CCTK_REAL x[], CCTK_REAL *Spritzrhomin, CCTK_REAL rf_prec,
                                            struct eosomnivars *eosvars, struct LocGlob *lgp,
                                            void (*funcd) (CCTK_REAL [], CCTK_REAL [],
                                                           CCTK_REAL [], CCTK_REAL [][3],
                                                           CCTK_REAL *, CCTK_REAL *,
                                                           struct eosomnivars *eosvars, struct LocGlob *) );

static  void func_vsq( CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], CCTK_REAL [][2],
                       CCTK_REAL *f, CCTK_REAL *df, CCTK_REAL gammaeos, struct LocGlob *lgp);

static  void func_vsq_eosomni( CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], CCTK_REAL [][3],
                               CCTK_REAL *f, CCTK_REAL *df, struct eosomnivars *eosvars, struct LocGlob *lgp);

static CCTK_REAL x1_of_x0(CCTK_REAL x0, struct LocGlob *lgp ) ;


// EOS STUFF:
static CCTK_REAL eos_info(CCTK_REAL W, CCTK_REAL vsq, CCTK_REAL *dpdw, CCTK_REAL *dpdvsq, CCTK_REAL gammaeos, struct LocGlob *lgp);
/*static CCTK_REAL eos_info_eosomni(CCTK_REAL W, CCTK_REAL vsq, CCTK_REAL eps0, struct LocGlob *lgp);*/

/* pressure as a function of rho0 and u */
//static CCTK_REAL pressure_rho0_u(CCTK_REAL rho0, CCTK_REAL u, CCTK_REAL gammaeos)
//{
//  return((gammaeos - 1.)*u) ;
//}

static CCTK_REAL pressure_rho0_eps_eosomni(CCTK_REAL rho0, CCTK_REAL eps, CCTK_REAL* dpdrho, CCTK_REAL* dpdeps, struct eosomnivars *eosvars);

/* Pressure as a function of rho0 and w = rho0 + u + p */
static CCTK_REAL pressure_rho0_w(CCTK_REAL rho0, CCTK_REAL w,CCTK_REAL gammaeos)
{
  return((gammaeos-1.)*(w - rho0)/gammaeos) ;
}

void CCTK_FCALL CCTK_FNAME(Spritz_Con2Prim_Noble) (
    CCTK_INT  *handle, CCTK_INT *rl, CCTK_INT *ii, CCTK_INT *jj, CCTK_INT *kk, CCTK_REAL *x, CCTK_REAL *y, CCTK_REAL *z, CCTK_INT *keytemp, CCTK_REAL *eos_prec, CCTK_REAL *prec,   /*  */
    CCTK_REAL *gamma_eos,
    CCTK_REAL *dens_in,
    CCTK_REAL *sx_in, CCTK_REAL *sy_in, CCTK_REAL *sz_in,
    CCTK_REAL *tau_in, CCTK_REAL *Bconsx_in, CCTK_REAL *Bconsy_in, CCTK_REAL *Bconsz_in,
    CCTK_REAL *y_e_in, CCTK_REAL *temp_in,
    CCTK_REAL *rho,
    CCTK_REAL *velx, CCTK_REAL *vely, CCTK_REAL *velz,
    CCTK_REAL *epsilon, CCTK_REAL *pressure,
    CCTK_REAL *Bx, CCTK_REAL *By, CCTK_REAL *Bz, CCTK_REAL *bsq,
    CCTK_REAL *w_lorentz,
    CCTK_REAL *gxx, CCTK_REAL *gxy, CCTK_REAL *gxz,
    CCTK_REAL *gyy, CCTK_REAL *gyz, CCTK_REAL *gzz,
    CCTK_REAL *uxx, CCTK_REAL *uxy, CCTK_REAL *uxz,
    CCTK_REAL *uyy, CCTK_REAL *uyz, CCTK_REAL *uzz,
    CCTK_REAL *det,
    CCTK_INT  *epsnegative,
    CCTK_REAL *Spritz_rho_min,
    CCTK_REAL *coldtab_rhomin, CCTK_REAL *coldtab_press_lo, CCTK_REAL *coldtab_eps_lo, CCTK_REAL *coldtab_gamma_lo,
    CCTK_REAL *retval);

/**********************************************************************/
/**********************************************************************************

  Con2PrimM_pt():
  -----------------------------

   -- Attempts an inversion from GRMHD conserved variables to primitive variables assuming a guess.

   -- Uses the 2D method of Noble et al. (2006):
       -- Solves for two independent variables (W,v^2) via a 2D
          Newton-Raphson method
       -- Can be used (in principle) with a general equation of state.

   -- Minimizes two residual functions using a homemade Newton-Raphson routine.
       -- It is homemade so that it can catch exceptions and handle them correctly, plus it is
          optimized for this particular  problem.

  -- Note that the notation used herein is that of Noble et al. (2006) except for the argument
     list.


INPUT:  (using Spritz variable defintions)

   s[x,y,z]   =  scons[0,1,2]  = \alpha \sqrt(\gamma) T^0_i
   dens, tau  =  as defined in Spritz and are assumed to be densitized (i.e. with sqrt(\gamma))
   dens       =  D = \sqrt(\gamma) W \rho
   tau        =  \alpha^2 \sqrt(\gamma) T^{00} - D
   g[x,y,z][x,y,x] = spatial metric corresponding to \gamma
   u[x,y,z][x,y,z] = inverse of the spatial metric, g[x,y,z][x,y,z]
   det        =  \gamma (determinant of spatial metric)
   B[x,y,z]   =  Bvec[0,1,2]
   bsq        = b^\mu b_\mu

   epsnegative = (integer)
               = 0  if rho and epsilon are positive
              != 0  otherwise


  --  (users should set B[x,y,z] = 0  for hydrodynamic runs)


OUTPUT:  (using Spritz variable defintions)
   rho, eps   =  as defined in Spritz, primitive variables
   vel[x,y,z] =  as defined in Spritz, primitive variables


RETURN VALUE: of retval = (i*100 + j)  where
         i = 0 ->  Newton-Raphson solver either was not called (yet or not used)
                   or returned successfully;
             1 ->  Newton-Raphson solver did not converge to a solution with the
                   given tolerances;
             2 ->  Newton-Raphson procedure encountered a numerical divergence
                   (occurrence of "nan" or "+/-inf" ;

         j = 0 -> success
             1 -> failure: some sort of failure in Newton-Raphson;
             2 -> failure: unphysical vsq = v^2  value at initial guess;
             3 -> failure: W<0 or W>W_TOO_BIG
             4 -> failure: v^2 > 1
             ( used to be  5 -> failure: rho,uu <= 0   but now sets epsnegative to non-zero )

**********************************************************************************/
void CCTK_FCALL CCTK_FNAME(Spritz_Con2Prim_Noble)  (
    CCTK_INT  *handle, CCTK_INT *rl, CCTK_INT *ii, CCTK_INT *jj, CCTK_INT *kk, CCTK_REAL *x, CCTK_REAL *y, CCTK_REAL *z, CCTK_INT *keytemp, CCTK_REAL *eos_prec, CCTK_REAL *prec,    /* */
    CCTK_REAL *gamma_eos,
    CCTK_REAL *dens_in,
    CCTK_REAL *sx_in, CCTK_REAL *sy_in, CCTK_REAL *sz_in,
    CCTK_REAL *tau_in,
    CCTK_REAL *Bconsx_in, CCTK_REAL *Bconsy_in, CCTK_REAL *Bconsz_in,
    CCTK_REAL *y_e_in, CCTK_REAL* temp_in,
    CCTK_REAL *rho,
    CCTK_REAL *velx, CCTK_REAL *vely, CCTK_REAL *velz,
    CCTK_REAL *epsilon, CCTK_REAL *pressure,
    CCTK_REAL *Bx, CCTK_REAL *By, CCTK_REAL *Bz,
    CCTK_REAL *bsq,
    CCTK_REAL *w_lorentz,
    CCTK_REAL *gxx, CCTK_REAL *gxy, CCTK_REAL *gxz,
    CCTK_REAL *gyy, CCTK_REAL *gyz, CCTK_REAL *gzz,
    CCTK_REAL *uxx, CCTK_REAL *uxy, CCTK_REAL *uxz,
    CCTK_REAL *uyy, CCTK_REAL *uyz, CCTK_REAL *uzz,
    CCTK_REAL *det,
    CCTK_INT  *epsnegative,
    CCTK_REAL *Spritz_rho_min,
    CCTK_REAL *coldtab_rhomin, CCTK_REAL *coldtab_press_lo, CCTK_REAL *coldtab_eps_lo, CCTK_REAL *coldtab_gamma_lo,
    CCTK_REAL *retval)

{
  CCTK_REAL x_3d[3];
  CCTK_REAL sx, sy, sz;
  CCTK_REAL usx, usy, usz;
  CCTK_REAL tau, dens, gammaeos;
  CCTK_REAL QdotB,VdotB;
  CCTK_REAL rho0,u,p,w,gammasq,gamma,gtmp,W_last,W,vsq;
  CCTK_REAL g_o_WBsq, QdB_o_W;
  CCTK_REAL detg = (*det);
  CCTK_REAL sqrt_detg = sqrt(detg);
  CCTK_REAL inv_sqrt_detg = 1./sqrt_detg;
  CCTK_INT i,j, i_increase ;
  CCTK_INT nf, nfudgemax, failwarnmode, failinfomode;
  CCTK_REAL value_corr;
  CCTK_CHAR message[100];
  CCTK_REAL rho_atmosph;
  DECLARE_CCTK_PARAMETERS;

  struct LocGlob lg;
  struct eosomnivars eosvars;

  eosvars.eoshandle = *handle;
  //printf("handle = %i\n",*handle);
  eosvars.eoskeytemp = *keytemp;
  eosvars.eosprec = *eos_prec;
  eosvars.prec = *prec;
  eosvars.eos_y_e[0] = *y_e_in;
  eosvars.eos_temp[0] = *temp_in;
  eosvars.eoskeyerr[0] = 0;
  eosvars.eosanyerr[0] = 0;
  eosvars.eosrl = *rl;
  eosvars.eosit = 1;
  eosvars.eosi = *ii;
  eosvars.eosj = *jj;
  eosvars.eosk = *kk;
  eosvars.eosx = *x;
  eosvars.eosy = *y;
  eosvars.eosz = *z;
  eosvars.coldtab_rhomin   = *coldtab_rhomin;
  eosvars.coldtab_press_lo = *coldtab_press_lo;
  eosvars.coldtab_eps_lo   = *coldtab_eps_lo;
  eosvars.coldtab_gamma_lo = *coldtab_gamma_lo;

  gammaeos = *gamma_eos;

  nf =0;
  nfudgemax = 30;
  failinfomode = 1;
  failwarnmode = 0;

  /* Assume ok initially: */
  value_corr = 0.0;
  *retval = value_corr;
  *epsnegative = 0;

#if(DEBUG_CON2PRIMM)
  fprintf(stdout," !!! ENTERING DEBUG !!! \n");
  fprintf(stdout," *x1          = %26.16e \n", *x        );
  fprintf(stdout," *y1          = %26.16e \n", *y         );
  fprintf(stdout," *z1          = %26.16e \n", *z         );
  fprintf(stdout," *dens        = %26.16e \n", *dens_in         );
  fprintf(stdout," *sx          = %26.16e \n", *sx_in     );
  fprintf(stdout," *sy          = %26.16e \n", *sy_in     );
  fprintf(stdout," *sz          = %26.16e \n", *sz_in     );
  fprintf(stdout," *tau         = %26.16e \n", *tau_in    );
  fprintf(stdout," *Bconsx      = %26.16e \n", *Bconsx_in );
  fprintf(stdout," *Bconsy      = %26.16e \n", *Bconsy_in );
  fprintf(stdout," *Bconsz      = %26.16e \n", *Bconsz_in );
  fprintf(stdout," *rho         = %26.16e \n", *rho       );
  fprintf(stdout," *velx        = %26.16e \n", *velx      );
  fprintf(stdout," *vely        = %26.16e \n", *vely      );
  fprintf(stdout," *velz        = %26.16e \n", *velz      );
  fprintf(stdout," *epsilon     = %26.16e \n", *epsilon     );
  fprintf(stdout," *temp_in     = %26.16e \n", *temp_in     );
  fprintf(stdout," *y_e_in      = %26.16e \n", *y_e_in     );
  fprintf(stdout," *pressure    = %26.16e \n", *pressure    );
  fprintf(stdout," *Bx          = %26.16e \n", *Bx        );
  fprintf(stdout," *By          = %26.16e \n", *By        );
  fprintf(stdout," *Bz          = %26.16e \n", *Bz        );
  fprintf(stdout," *bsq         = %26.16e \n", *bsq       );
  fprintf(stdout," *w_lorentz   = %26.16e \n", *w_lorentz   );
  fprintf(stdout," *gxx         = %26.16e \n", *gxx       );
  fprintf(stdout," *gxy         = %26.16e \n", *gxy       );
  fprintf(stdout," *gxz         = %26.16e \n", *gxz       );
  fprintf(stdout," *gyy         = %26.16e \n", *gyy       );
  fprintf(stdout," *gyz         = %26.16e \n", *gyz       );
  fprintf(stdout," *gzz         = %26.16e \n", *gzz       );
  fprintf(stdout," *uxx         = %26.16e \n", *uxx       );
  fprintf(stdout," *uxy         = %26.16e \n", *uxy       );
  fprintf(stdout," *uxz         = %26.16e \n", *uxz       );
  fprintf(stdout," *uyy         = %26.16e \n", *uyy       );
  fprintf(stdout," *uyz         = %26.16e \n", *uyz       );
  fprintf(stdout," *uzz         = %26.16e \n", *uzz       );
  fprintf(stdout," *det         = %26.16e \n", *det       );
  fprintf(stdout," *epsnegative = %10d    \n", *epsnegative );
  fprintf(stdout," *retval      = %26.16e \n", *retval      );
  fprintf(stdout," !!!!!!!!!!!!!!!!!!!!!! \n");
  fflush(stdout);
#endif

  /* First undensitize all conserved variables : */
  sx   = (  *sx_in)   * inv_sqrt_detg;
  sy   = (  *sy_in)   * inv_sqrt_detg;
  sz   = (  *sz_in)   * inv_sqrt_detg;
  tau  = ( *tau_in)   * inv_sqrt_detg;
  dens = (*dens_in)   * inv_sqrt_detg;

  usx  =  (*uxx)*sx + (*uxy)*sy + (*uxz)*sz;
  usy  =  (*uxy)*sx + (*uyy)*sy + (*uyz)*sz;
  usz  =  (*uxz)*sx + (*uyz)*sy + (*uzz)*sz;

  *Bx = (*Bconsx_in) * inv_sqrt_detg;
  *By = (*Bconsy_in) * inv_sqrt_detg;
  *Bz = (*Bconsz_in) * inv_sqrt_detg;

  // Calculate various scalars (Q.B, Q^2, etc)  from the conserved variables:

  lg.Bsq =
    (*gxx) * (*Bx) * (*Bx) +
    (*gyy) * (*By) * (*By) +
    (*gzz) * (*Bz) * (*Bz) +
    2*(
       (*gxy) * (*Bx) * (*By) +
       (*gxz) * (*Bx) * (*Bz) +
       (*gyz) * (*By) * (*Bz) );

  QdotB = (sx * (*Bx) + sy * (*By) + sz * (*Bz)) ;

  VdotB = ((*gxx) * (*velx) + (*gxy) * (*vely) + (*gxz) * (*velz))* (*Bx) +
           ((*gxy) * (*velx) + (*gyy) * (*vely) + (*gyz) * (*velz))* (*By) +
           ((*gxz) * (*velx) + (*gyz) * (*vely) + (*gzz) * (*velz))* (*Bz);

  lg.VdotBsq = VdotB*VdotB;

  lg.QdotBsq = QdotB*QdotB ;

  lg.Qdotn = -(tau + dens) ;

  lg.Qtsq = (usx * sx  +  usy * sy  +  usz * sz) ;

  lg.D = dens;

  lg.half_Bsq = 0.5*lg.Bsq;

  lg.tau = tau;

  /* calculate W from last timestep and use for guess */
  vsq =
    (*gxx) * (*velx) * (*velx) +
    (*gyy) * (*vely) * (*vely) +
    (*gzz) * (*velz) * (*velz) +
    2*(
       (*gxy) * (*velx) * (*vely) +
       (*gxz) * (*velx) * (*velz) +
       (*gyz) * (*vely) * (*velz) );

  if( (vsq < 0.) && (fabs(vsq) < 1.0e-13) ) {
    vsq = fabs(vsq);
  }
  if(vsq < 0. || vsq > 1. ) {
    value_corr = 2.0;
    *retval = value_corr;
/*#if(DEBUG_CON2PRIMM)
    fprintf(stdout," *retval (inside loop)      = %26.16e \n", *retval      );
    fflush(stdout);
#endif*/
    return;
  }

  gammasq = 1. / (1. - vsq);
  gamma  = sqrt(gammasq);

  // Always calculate rho from D and gamma so that using D in EOS remains consistent
  //   i.e. you don't get positive values for dP/d(vsq) .
  rho0 = lg.D / gamma ;
  u = (*epsilon) * rho0;
  CCTK_REAL uold = u;

  CCTK_REAL dum1,dum2;
//  eosvars.eosanyerr[0]=5;
  p = pressure_rho0_eps_eosomni(rho0,*epsilon,&dum1,&dum2,&eosvars) ;  // EOSOMNI
  //  p = pressure_rho0_u(rho0,u,gammaeos) ;  // EOS
  // error handling
  if (eosvars.eosanyerr[0] != 0) {
    if(eosvars.eoskeyerr[0] == 668) {
      value_corr = 3.0;
      *retval = value_corr;
    }
  }

  w = rho0 + u + p;

  // we need this for the hot EOS failsafe tau reset later

  lg.bsq = lg.Bsq/(w*w) + lg.VdotBsq;

  W_last = w*gammasq ;

  lg.tau = tau;

  // Make sure that W is large enough so that v^2 < 1 :
  i_increase = 0;
  while( (( W_last*W_last*W_last * ( W_last + 2.*lg.Bsq )
	    - lg.QdotBsq*(2.*W_last + lg.Bsq) ) <= W_last*W_last*(lg.Qtsq-lg.Bsq*lg.Bsq))
	 && (i_increase < 10) ) {
    W_last *= 10.;
    i_increase++;
  }

  // Calculate W and vsq:
  x_3d[0] =  fabs( W_last );
  x_3d[1] = x1_of_x0( W_last, &lg ) ;

  //Use 2d NR for polytropes!
  if (*handle==1 || *handle==2) {
    // MODIFICATION //
    x_3d[2] = u;
    // END MODIFICATION //
    //rho_atmosph = rho_abs_min;
    value_corr = 1.0*threed_newton_raphson_omni( x_3d, Spritz_rho_min, Eos_rf_prec, &eosvars, &lg, func_vsq_eosomni ) ;   //twod_newton_raphson( x_3d, gammaeos, &lg, func_vsq )
    *retval = value_corr;

  } else {
    //USE 3d NR for non-polytropes!
    x_3d[2] = u;
    //rho_atmosph = rho_abs_min;
    value_corr = 1.0*threed_newton_raphson_omni( x_3d, Spritz_rho_min, Eos_rf_prec, &eosvars, &lg, func_vsq_eosomni ) ;
    *retval = value_corr;
  }

  W = x_3d[0];
  vsq = x_3d[1];

  /* Problem with solver, so return denoting error before doing anything further */
  if( ((*retval) != 0.) || (W == FAIL_VAL) ) {
    value_corr = value_corr*100.+1.;
    *retval = value_corr*100.+1.;
#if(DEBUG_CON2PRIMM)
    fprintf(stdout," !!! SECOND PRINT !!! \n" );
    fprintf(stdout," *retval      = %26.16e \n", *retval      );
    fprintf(stdout," *dens        = %26.16e \n", *dens_in         );
    fprintf(stdout," *sx          = %26.16e \n", *sx_in     );
    fprintf(stdout," *sy          = %26.16e \n", *sy_in     );
    fprintf(stdout," *sz          = %26.16e \n", *sz_in     );
    fprintf(stdout," *tau         = %26.16e \n", *tau_in    );
    fprintf(stdout," *Bconsx      = %26.16e \n", *Bconsx_in );
    fprintf(stdout," *Bconsy      = %26.16e \n", *Bconsy_in );
    fprintf(stdout," *Bconsz      = %26.16e \n", *Bconsz_in );
    fprintf(stdout," *rho         = %26.16e \n", *rho       );
    fprintf(stdout," *velx        = %26.16e \n", *velx      );
    fprintf(stdout," *vely        = %26.16e \n", *vely      );
    fprintf(stdout," *velz        = %26.16e \n", *velz      );
    fprintf(stdout," *epsilon     = %26.16e \n", *epsilon     );
    fprintf(stdout," *temp_in     = %26.16e \n", *temp_in     );
    fprintf(stdout," *y_e_in      = %26.16e \n", *y_e_in     );
    fprintf(stdout," *pressure    = %26.16e \n", *pressure    );
    fprintf(stdout," *Bx          = %26.16e \n", *Bx        );
    fprintf(stdout," *By          = %26.16e \n", *By        );
    fprintf(stdout," *Bz          = %26.16e \n", *Bz        );
    fprintf(stdout," *bsq         = %26.16e \n", *bsq       );
    fprintf(stdout," *w_lorentz   = %26.16e \n", *w_lorentz   );
    fprintf(stdout," *gxx         = %26.16e \n", *gxx       );
    fprintf(stdout," *gxy         = %26.16e \n", *gxy       );
    fprintf(stdout," *gxz         = %26.16e \n", *gxz       );
    fprintf(stdout," *gyy         = %26.16e \n", *gyy       );
    fprintf(stdout," *gyz         = %26.16e \n", *gyz       );
    fprintf(stdout," *gzz         = %26.16e \n", *gzz       );
    fprintf(stdout," *uxx         = %26.16e \n", *uxx       );
    fprintf(stdout," *uxy         = %26.16e \n", *uxy       );
    fprintf(stdout," *uxz         = %26.16e \n", *uxz       );
    fprintf(stdout," *uyy         = %26.16e \n", *uyy       );
    fprintf(stdout," *uyz         = %26.16e \n", *uyz       );
    fprintf(stdout," *uzz         = %26.16e \n", *uzz       );
    fprintf(stdout," *det         = %26.16e \n", *det       );
    fprintf(stdout," *epsnegative = %10d    \n", *epsnegative );
    fprintf(stdout," !!!!!!!!!!!!!!!!!!!!!! \n");
    fflush(stdout);
#endif
    return;
  }
  else{
    if(W <= 0. || W > W_TOO_BIG) {
      value_corr = 3.0;
      *retval = value_corr;
/*#if(DEBUG_CON2PRIMM)
      fprintf(stdout," *retval      = %26.16e \n", *retval      );
      fflush(stdout);
#endif*/
      return;
    }
  }

  // Calculate v^2:
  if( vsq >= 1. ) {
    value_corr = 4.0;
    *retval = value_corr;
/*#if(DEBUG_CON2PRIMM)
    fprintf(stdout," *retval      = %26.16e \n", *retval      );
    fflush(stdout);
#endif*/
    return;
  }

  // Recover the primitive variables from the scalars and conserved variables:
  gtmp = sqrt(1. - vsq);
  gamma = 1./gtmp ;
  rho0 = lg.D * gtmp;

  w = W * (1. - vsq) ;
  if (*handle==1 || *handle==2) {
    /*p = pressure_rho0_w(rho0,w,gammaeos) ;  // EOS
    u = w - (rho0 + p) ;
    *epsilon = u / rho0;
    // MODIFICATION //
    x_3d[2] = u;
    // END MODIFICATION // */              //  BIG MODIFICATION !!!!!!
    u=x_3d[2];
    *epsilon = u/rho0;
    tau = lg.tau;
    p = pressure_rho0_eps_eosomni(rho0,*epsilon,&dum1,&dum2,&eosvars) ;  // EOSOMNI
    *temp_in = eosvars.eos_temp[0];
  } else {
    u=x_3d[2];
    *epsilon = u/rho0;
    tau = lg.tau;
    p = pressure_rho0_eps_eosomni(rho0,*epsilon,&dum1,&dum2,&eosvars) ;  // EOSOMNI
    *temp_in = eosvars.eos_temp[0];
  }

  // update tau in case we have modified it in C2P

  *tau_in = tau*sqrt_detg;

  // User may want to handle this case differently, e.g. do NOT return upon
  // a negative rho/u, calculate v^i so that rho/u can be floored by other routine:
  if( (rho0 <= 0.) || (u < 0.) ) {   // =
    *epsnegative = 1;
    sprintf(message," PROBLEMS: ");
    CCTK_WARN(1,message);
    sprintf(message," Point is x = %12.10e y = %12.10e z= %12.10e", x,y,z);
    CCTK_WARN(1,message);
    sprintf(message," First variable is W = %12.10e and x_3d[0] = %12.10e", W, x_3d[0]);
    CCTK_WARN(1,message);
    sprintf(message," Third variable is vsq = %12.10e and x_3d[1] = %12.10e", vsq, x_3d[1]);
    CCTK_WARN(1,message);
    sprintf(message," Third variable is u = %12.10e and x_3d[2] = %12.10e", u, x_3d[2]);
    CCTK_WARN(1,message);
    sprintf(message," rho0 = %12.10e ", rho0);
    CCTK_WARN(1,message);
    sprintf(message," epsilon = %12.10e ", *epsilon);
    CCTK_WARN(1,message);
    sprintf(message," Handle is %d", (int)*handle);
    CCTK_WARN(1,message);
    CCTK_WARN(Spritz_c2p_warnlevel, "NEGATIVE ENERGY DENSITY!");
    return;
  }

  *rho = rho0;
  // *epsilon = u / rho0;
  *w_lorentz = gamma;
  *pressure = p ;

  g_o_WBsq = 1./(W+lg.Bsq);
  QdB_o_W  = QdotB / W;
  *bsq = lg.Bsq * (1.-vsq) + QdB_o_W*QdB_o_W;

  *velx = g_o_WBsq * ( usx + QdB_o_W*(*Bx) ) ;
  *vely = g_o_WBsq * ( usy + QdB_o_W*(*By) ) ;
  *velz = g_o_WBsq * ( usz + QdB_o_W*(*Bz) ) ;

  if (*rho <= rho_abs_min*(1.0+Spritz_atmo_tolerance) ) {
    *rho = rho_abs_min;
    *velx = 0.0;
    *vely = 0.0;
    *velz = 0.0;
    *w_lorentz = 1.0;
  }

#if(DEBUG_CON2PRIMM)
  fprintf(stdout," !!! OUTPUT !!! \n" );
  fprintf(stdout," *retvalo      = %26.16e \n", *retval      );
  fprintf(stdout," *x2        = %26.16e \n", *x        );
  fprintf(stdout," *y2        = %26.16e \n", *y         );
  fprintf(stdout," *z2        = %26.16e \n", *z         );
  fprintf(stdout," *denso        = %26.16e \n", *dens_in         );
  fprintf(stdout," *sxo                 = %26.16e \n", *sx_in     );
  fprintf(stdout," *syo         = %26.16e \n", *sy_in     );
  fprintf(stdout," *szo                 = %26.16e \n", *sz_in     );
  fprintf(stdout," *tauo        = %26.16e \n", *tau_in    );
  fprintf(stdout," *Bconsxo     = %26.16e \n", *Bconsx_in );
  fprintf(stdout," *Bconsyo     = %26.16e \n", *Bconsy_in );
  fprintf(stdout," *Bconszo      = %26.16e \n", *Bconsz_in );
  fprintf(stdout," *rhoo         = %26.16e \n", *rho      );
  fprintf(stdout," *velxo        = %26.16e \n", *velx     );
  fprintf(stdout," *velyo        = %26.16e \n", *vely     );
  fprintf(stdout," *velzo        = %26.16e \n", *velz     );
  fprintf(stdout," *epsilono     = %26.16e \n", *epsilon     );
  fprintf(stdout," *temp_ino     = %26.16e \n", *temp_in     );
  fprintf(stdout," *y_e_ino     = %26.16e \n", *y_e_in     );
  fprintf(stdout," *pressureo    = %26.16e \n", *pressure    );
  fprintf(stdout," *Bxo                 = %26.16e \n", *Bx        );
  fprintf(stdout," *Byo                 = %26.16e \n", *By        );
  fprintf(stdout," *Bzo                 = %26.16e \n", *Bz        );
  fprintf(stdout," *bsqo                = %26.16e \n", *bsq       );
  fprintf(stdout," *w_lorentzo   = %26.16e \n", *w_lorentz   );
  fprintf(stdout," *gxxo        = %26.16e \n", *gxx       );
  fprintf(stdout," *gxyo        = %26.16e \n", *gxy       );
  fprintf(stdout," *gxzo        = %26.16e \n", *gxz       );
  fprintf(stdout," *gyyo        = %26.16e \n", *gyy       );
  fprintf(stdout," *gyzo        = %26.16e \n", *gyz       );
  fprintf(stdout," *gzzo        = %26.16e \n", *gzz       );
  fprintf(stdout," *uxxo        = %26.16e \n", *uxx       );
  fprintf(stdout," *uxyo        = %26.16e \n", *uxy       );
  fprintf(stdout," *uxzo                = %26.16e \n", *uxz       );
  fprintf(stdout," *uyyo        = %26.16e \n", *uyy       );
  fprintf(stdout," *uyzo        = %26.16e \n", *uyz       );
  fprintf(stdout," *uzzo                = %26.16e \n", *uzz       );
  fprintf(stdout," *deto                = %26.16e \n", *det       );
  fprintf(stdout," *epsnegativeo = %10d    \n", *epsnegative );
  fprintf(stdout," !!!!!!!!!!!!!!! \n");
  fflush(stdout);
#endif

  /* done! */
  return;

}


/**********************************************************************/
/****************************************************************************
   vsq_calc():

      -- evaluate v^2 (spatial, normalized velocity) from
            W = \gamma^2 w

****************************************************************************/
static CCTK_REAL vsq_calc(CCTK_REAL W, struct LocGlob *lgp)
{
	CCTK_REAL Wsq,Xsq,Bsq_W;

	Wsq = W*W ;
	Bsq_W = (lgp->Bsq + W);
	Xsq = Bsq_W * Bsq_W;

	return(  ( Wsq * lgp->Qtsq  + lgp->QdotBsq * (Bsq_W + W)) / (Wsq*Xsq) ); // eq (28) Noble et al (2006)
}                                                                                // also eq (23) Del Zanna et al (2003)


/********************************************************************

  x1_of_x0():

    -- calculates v^2 from W  with some physical bounds checking;
    -- asumes x0 is already physical
    -- makes v^2 physical  if not;

*********************************************************************/

static CCTK_REAL x1_of_x0(CCTK_REAL x0, struct LocGlob *lgp )
{
  CCTK_REAL vsq;
  CCTK_REAL dv = 1.e-15;

  vsq = fabs(vsq_calc(x0,lgp)) ; // guaranteed to be positive

  return( ( vsq > 1. ) ? (1.0 - dv) : vsq   );

}

/********************************************************************

  validate_x():

    -- makes sure that x[0,1] have physical values, based upon
       their definitions:

*********************************************************************/

static void validate_x(CCTK_REAL x[2], CCTK_REAL x0[2] )
{

  const CCTK_REAL dv = 1.e-15;

  /* Always take the absolute value of x[0] and check to see if it's too big:  */
  x[0] = fabs(x[0]);
  x[0] = (x[0] > W_TOO_BIG) ?  x0[0] : x[0];

  x[1] = (x[1] < 0.) ?   0.       : x[1];  /* if it's too small */
  x[1] = (x[1] > 1.) ?  (1. - dv) : x[1];  /* if it's too big   */

  return;

}

/************************************************************

  twod_newton_raphson():

    -- performs Newton-Rapshon method on an 2d system for polytropes.

    -- inspired in part by Num. Rec.'s routine newt();

*****************************************************************/
static CCTK_INT twod_newton_raphson( CCTK_REAL x[], CCTK_REAL gammaeos, struct LocGlob *lgp,
			    void (*funcd) (CCTK_REAL [], CCTK_REAL [], CCTK_REAL [],
					   CCTK_REAL [][2], CCTK_REAL *,
					   CCTK_REAL *, CCTK_REAL, struct LocGlob *) )
{
  CCTK_REAL f, df, dx[2], x_old[2];
  CCTK_REAL resid[2], jac[2][2];
  CCTK_REAL errx, x_orig[2];
  CCTK_INT    n_iter, i_extra, doing_extra;
  CCTK_REAL vsq,vsq_old,W,W_old;
  const CCTK_REAL dv = (1.-1.e-15);

  CCTK_INT   keep_iterating;

  //fprintf(stdout," !!!! TWO D N-R USED !!!! \n");

  // Initialize various parameters and variables:
  errx = 1. ;
  df = f = 1.;
  i_extra = doing_extra = 0;
  x_old[0] = x_orig[0] = x[0] ;
  x_old[1] = x_orig[1] = x[1] ;

  vsq_old = vsq = W = W_old = 0.;
  n_iter = 0;

  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) {

    (*funcd) (x, dx, resid, jac, &f, &df, gammaeos, lgp);  /* returns with new dx, f, df */


    /* Save old values before calculating the new: */
    errx = 0.;
    x_old[0] = x[0] ;
    x_old[1] = x[1] ;

    /* Make the newton step: */
    x[0] += (NEWT_RELAX*dx[0])  ;    //
    x[1] += (NEWT_RELAX*dx[1])  ;    //

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);


    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    if( x[0] < 0. ) {  x[0] = fabs(x[0]);  }
    else {
     if(x[0] > W_TOO_BIG)  { x[0] = x_old[0] ; }
    }

    if( x[1] < 0. ) {  x[1] = 0.; }
    else {
      if( x[1] > 1. ) { x[1] = dv; }
    }

    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/

    if( (fabs(errx) <= MIN_NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {  // NEWT_TOL
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0))
        || (i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)


    /*  Check for bad untrapped divergences : */
  if( (!finite(f)) ||  (!finite(df)) ) {
    return(2);
  }


  if( fabs(errx) <= NEWT_TOL ){
    return(0);
  }
  else if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    return(0);
  }
  else {
    return(1);
  }

  return(0);

}


/************************************************************

  threed_newton_raphson_omni():

    -- performs Newton-Rapshon method on an 2d system for polytropes.

    -- inspired in part by Num. Rec.'s routine newt();

*****************************************************************/
static CCTK_INT threed_newton_raphson_omni( CCTK_REAL x[], CCTK_REAL *Spritzrhomin, CCTK_REAL rf_prec,
                                            struct eosomnivars *eosvars, struct LocGlob *lgp,
                                            void (*funcd) (CCTK_REAL [], CCTK_REAL [], CCTK_REAL [],
                                            CCTK_REAL [][3], CCTK_REAL *,
                                            CCTK_REAL *, struct eosomnivars *, struct LocGlob *) )
{
  DECLARE_CCTK_PARAMETERS;
  CCTK_REAL f, df, dx[3], x_old[3];
  CCTK_REAL resid[3], jac[3][3];
  CCTK_REAL errx, x_orig[3], x_entrance[3];
  CCTK_INT    n_iter, id, jd, i_extra, doing_extra;
  CCTK_REAL dW,dvsq,du,vsq_old,vsq,W,W_old,u,u_old;
  const CCTK_REAL dv = (1.-1.e-15);
  double relax_factor;
  CCTK_CHAR messnr[200];

  CCTK_INT   keep_iterating;
  CCTK_REAL rho_atm, press_atm, eps_atm;
  /*!!!!!!
  ! begin EOS Omni vars*/
  /*CCTK_REAL  Eos_rf_prec;*/
  CCTK_INT n,keytemp,anyerr,keyerr[1];  /* * */
  CCTK_REAL xrho[1],xpress[1],xeps[1];

  xpress[0] = 0.0;
  xeps[0] = eps_dummy;
  /*! end EOS Omni vars
  !!!!!!*/

  //fprintf(stdout," !!!! THREE D N-R USED !!!! \n");

  // Initialize various parameters and variables:
  errx = 1. ;
  df = f = 1.;
  i_extra = doing_extra = 0;
  x_old[0] = x_orig[0] = x[0] ;
  x_old[1] = x_orig[1] = x[1] ;
  x_old[2] = x_orig[2] = x[2] ;

  relax_factor = 1.0;
  //use_relax = 0;

  vsq_old = vsq = W = W_old = u = u_old = 0.;
  n_iter = 0;

  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) {

    /* Store entrance values */
    x_entrance[0] = x[0] ;
    x_entrance[1] = x[1] ;
    x_entrance[2] = x[2] ;

    (*funcd) (x, dx, resid, jac, &f, &df, eosvars, lgp);  /* returns with new dx, f, df */

    /* Save old values before calculating the new: */
    errx = 0.;
    x_old[0] = x[0] ;
    x_old[1] = x[1] ;
    x_old[2] = x[2] ;

    /* Make the newton step: */
    x[0] += (relax_factor*dx[0])  ;
    x[1] += (relax_factor*dx[1])  ;
    x[2] += (relax_factor*dx[2])  ;

    //printf("Updating vars: %g %g %g %g %g %g\n",x[0],dx[0],x[1],dx[1],x[2],dx[2]);

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);


    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    if( x[0] < 0. ) {
      x[0] = fabs(x[0]);

      //sprintf(messnr," ### W UNPHYS, TOO SMALL ");
      //CCTK_WARN(1,messnr);
      }
    else {
     if(x[0] > W_TOO_BIG)  {
       x[0] = x_old[0] ;

       //sprintf(messnr," ### W UNPHYS, TOO BIG ");
       //CCTK_WARN(1,messnr);
       }
    }

    if( x[1] < 0. ) {
      x[1] = 0.;

      //sprintf(messnr," ### vsq UNPHYS, TOO SLOW ");
      //CCTK_WARN(1,messnr);
      }
    else {
      if( x[1] > 1. ) {
        x[1] = dv;

        //sprintf(messnr," ### vsq UNPHYS, TOO FAST ");
        //CCTK_WARN(1,messnr);
      }
    }

    /* In case we obtain negative u, try to apply a relaxation factor to the Newton method... */
    if( x[2] < 0. ) {
      if(relax_factor > MIN_RELAX) {
        // Reset the initial values to entrance ones and set value for relax_factor...
        errx = 1. ;
        x[0] = x_entrance[0] ;
        x[1] = x_entrance[1] ;
        x[2] = x_entrance[2] ;
        keep_iterating = 1;
        relax_factor = relax_factor * FACTOR;
        n_iter = n_iter - 1; }   // Repeat the current iteration
        else {
        // Set u to value at atmosphere

        //x[2] = 0.;

        //sprintf(messnr," ### u UNPHYS ==> SETTING TO ATMOSPHERE VALUE IN N-R ");
        //CCTK_WARN(1,messnr);
        // Compute values of epsilon and pressure at atmosphere density.
        //
        // For cold_tabulated: use Spritz's own matched-polytrope
        // extrapolation below the table's rho_min (calibrated once, at
        // startup, in Spritz_Rho_Minima_Setup_Final / Spritz_Minima.cc, into
        // eosvars->coldtab_*) -- zero EOS_Omni calls here. This replaces the
        // EOS_Omni_press(keytemp=1) call that used to sit here: that relied
        // on EOS_Omni's own extrapolation, which hardcodes Gamma=2
        // regardless of the table's actual local slope (confirmed against
        // EOS_Omni_ColdEOSReadTable.F90), unlike IllinoisGRMHD/WhiskyTHC,
        // which both use the table's real local Gamma. The formula below is
        // the pure power-law ("multiplicative") form eps=eps_lo*ratio^
        // (Gamma_lo-1), P=P_lo*ratio^Gamma_lo, matching GRHayL's
        // NRPyEOS_hybrid_coldtab_compute_P_cold_and_eps_cold.c exactly --
        // verified numerically that the alternative additive form
        // (integrating deps=-P d(1/rho) from the matching point) goes
        // NEGATIVE at Spritz_rho_min for this table, consistent with that
        // file's own comment explaining why it was rejected there. The
        // multiplicative form is always positive and exactly continuous in
        // both P and eps at the matching point (ratio=1) by construction.
        //
        // For any other EOS table, fall back to the original EOS_Omni_press
        // call (only cold_tabulated has this extrapolation-Gamma bug).
        rho_atm = *Spritzrhomin;
        if (CCTK_EQUALS(Spritz_eos_table,"cold_tabulated")) {
          const CCTK_REAL Gamma_lo = eosvars->coldtab_gamma_lo;
          const CCTK_REAL ratio = rho_atm / eosvars->coldtab_rhomin;
          press_atm = eosvars->coldtab_press_lo * pow(ratio, Gamma_lo);
          eps_atm   = eosvars->coldtab_eps_lo * pow(ratio, Gamma_lo - 1.0);
        } else {
          xrho[0] = rho_atm;
          CCTK_INT keytemp_atm = 1;
          EOS_Omni_press(eosvars->eoshandle,keytemp_atm,eosvars->eosprec,1,xrho,xeps,eosvars->eos_temp,
                         eosvars->eos_y_e,xpress,eosvars->eoskeyerr,eosvars->eosanyerr);
          press_atm = xpress[0];
          eps_atm = xeps[0];
        }

        x[2] = eps_atm * rho_atm;
        keep_iterating = 1;
        // Set values also for other variables to the one of atmosphere...
        x[1] = 0.0;
        x[0] = (x[2] + rho_atm + press_atm)*(1.0/(1.0-x[1]));}

        doing_extra == 1;

        }

    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/

    if( (fabs(errx) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {  //eosvars->prec
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0))   //eosvars->prec
        || (i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)

  relax_factor = 1.0;

    /*  Check for bad untrapped divergences : */
  if( (!finite(f)) ||  (!finite(df)) ) {
    return(2);
  }


  if( fabs(errx) <= NEWT_TOL ){     // eosvars->prec
    return(0);
  }
  else if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){  //eosvars->prec   eosvars->prec
    return(0);
  }
  else {
    return(1);
  }

  return(0);

}

/**********************************************************************/
/*********************************************************************************
   func_vsq():

        -- calculates the residuals, and Newton step for general_newton_raphson();
        -- for this method, x=W,vsq here;

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
 *********************************************************************************/

static void func_vsq(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[],
		     CCTK_REAL jac[][2], CCTK_REAL *f, CCTK_REAL *df, CCTK_REAL gammaeos, struct LocGlob *lgp)
{


  CCTK_REAL  W, Wsq, vsq, p_tmp, dPdvsq, dPdW;
  CCTK_REAL t11, t16,t18,t2,t21,   t23,   t24,   t25,   t3,   t35,   t4,   t40;


  W = x[0];
  vsq = x[1];

  Wsq = W*W;

  p_tmp = eos_info(W, vsq, &dPdW, &dPdvsq, gammaeos, lgp);

  // These expressions were calculated using Mathematica, but made into efficient
  // code using Maple.  Since we know the analytic form of the equations, we can
  // explicitly calculate the Newton-Raphson step:

  t2 = -lgp->half_Bsq+dPdvsq;
  t3 = lgp->Bsq+W;
  t4 = t3*t3;
  t23 = 1/W;
  t16 = lgp->QdotBsq*t23*t23;
  t11 = lgp->Qtsq-vsq*t4+t16*(lgp->Bsq+W+W);
  t18 = -lgp->Qdotn-lgp->half_Bsq*(1.0+vsq)+0.5*t16-W+p_tmp;
  t24 = t16*t23;
  t25 = -1.0+dPdW-t24;
  t35 = t25*t3+(lgp->Bsq-2.0*dPdvsq)*(t16+vsq*W)*t23;
  //  t21 = 1/t3;
  //  t36 = 1/t35;
  t21 = 1/(t3*t35);
  dx[0] = -(t2*t11+t4*t18)*t21;
  t40 = -2*(vsq+t24)*t3;
  dx[1] = -(-t25*t11+t40*t18)*t21;
  //  detJ = t3*t35;
  jac[0][0] = t40;
  jac[0][1] = -t4;
  jac[1][0] = t25;
  jac[1][1] = t2;
  resid[0] = t11;
  resid[1] = t18;

  *df = -resid[0]*resid[0] - resid[1]*resid[1];

  *f = -0.5 * ( *df );

}

/**********************************************************************/
/*********************************************************************************
   func_vsq_eosomni():

        -- calculates the residuals, and Newton step for general_newton_raphson();
        -- for this method, x=W,vsq,u here;

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
 *********************************************************************************/

static void func_vsq_eosomni(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[],
                             CCTK_REAL jac[][3], CCTK_REAL *f, CCTK_REAL *df,
                             struct eosomnivars *eosvars, struct LocGlob *lgp)
{

  CCTK_REAL p_tmp;
  CCTK_REAL  W, vsq, u, epsilon,Wsq, dPdvsq, dPdW;
  CCTK_REAL res0, Lor, LorInv, rho0, Winv, QB2Winv2, drho0_dv, res1, j11, detJ,detJinv;
  CCTK_REAL QB2Winv3, j10,B2plusW,detJ_gcf,t36,B2plusW_sq, j00,j01,j12,j20,j21,j22;
  CCTK_REAL dpress_dv,dpress_du,dpdrho,dpdeps,c00,c01,c02,c10,c11,c12,c20,c21,c22;

  CCTK_INT nf, nfudgemax, failwarnmode, failinfomode;
  DECLARE_CCTK_PARAMETERS;

  nf =0;
  nfudgemax =30;
  failinfomode =2;
  failwarnmode =1;

  W = x[0];
  vsq = x[1];
  u = x[2];

  Wsq = W*W;
  Lor = 1.0/sqrt(1.0-vsq);
  LorInv = sqrt(1.0-vsq);
  rho0 = lgp->D * LorInv;
  epsilon = u/rho0;
  p_tmp = pressure_rho0_eps_eosomni(rho0,epsilon,&dpdrho,&dpdeps,eosvars);



  lgp->Qdotn = -(lgp->tau + lgp->D);
  B2plusW = lgp->Bsq+W;
  B2plusW_sq = B2plusW*B2plusW;
  Winv = 1/W;
  QB2Winv2 = lgp->QdotBsq*Winv*Winv;
  QB2Winv3=QB2Winv2*Winv;

  //Eq. 4: Qtsq-v^2(B^2+W)^2-QdotBsq(B^2+2W)/W^2=0 <-No u or p dependence

  resid[0] = lgp->Qtsq - vsq*B2plusW_sq+QB2Winv2*(lgp->Bsq+W+W);

  j00 = -2*(vsq+QB2Winv3)*B2plusW;
  j01 = -1.0*B2plusW_sq;


  //Eq. 5: -Qdotn - B^2(1+vsq)/2 + QdotBsq/2/W^2 - vsq W - rho_0(vsq) - u = 0
  //rho0 = D * sqrt(1-vsq)

  resid[1] = -lgp->Qdotn - lgp->half_Bsq*(1.0+vsq) + 0.5*QB2Winv2 - vsq*W - rho0 - u;

  drho0_dv = -0.5*lgp->D / LorInv;

  j10 = -vsq-QB2Winv3;
  j11 = -lgp->half_Bsq - W - drho0_dv;

  //Eq. 6: u+ p - W(1-vsq) + rho0 = 0  =>  p-W = -W vsq - rho0 - u

  resid[2] = u + p_tmp - W*(1.0-vsq) + rho0;

  //dp/dv = (dp/drho)u *  drho/dv
  //dp/drho_u = dp/drho_eps - eps/rho0 dpdeps
  dpress_dv = drho0_dv*dpdrho+
    u/2.0/lgp->D*pow(1.0-vsq,-1.5)*dpdeps;

  //dp/du = 1/rho dp deps, since rho0 is function of v only
  dpress_du=dpdeps/rho0;

  jac[0][0] = j00;
  jac[0][1] = j01;
  jac[0][2] = 0.0;
  jac[1][0] = j10;
  jac[1][1] = j11;
  jac[1][2] = -1.0;
  jac[2][0] = vsq-1.0;
  jac[2][1] = dpress_dv + W + drho0_dv;
  jac[2][2] = dpress_du + 1.0;

  c00 = j11*jac[2][2]+jac[2][1];
  c01 = -1*j01*jac[2][2];
  c02 = j01*jac[1][2];

  c10 = -jac[2][0]-j10*jac[2][2];
  c11 = j00*jac[2][2];
  c12 = j00;

  c20 = j10*jac[2][1]-j11*jac[2][0];
  c21 = j01*jac[2][0]-j00*jac[2][1];
  c22 = j00*j11-j01*j10;

  detJ=j00*c00+j01*c10;
  detJinv = 1/detJ;

  dx[0] = -(c00*resid[0]+c01*resid[1]+c02*resid[2])*detJinv;
  dx[1] = -(c10*resid[0]+c11*resid[1]+c12*resid[2])*detJinv;
  dx[2] = -(c20*resid[0]+c21*resid[1]+c22*resid[2])*detJinv;

  *df = -resid[0]*resid[0] - resid[1]*resid[1] - resid[2]*resid[2];

  *f = -0.5 * ( *df );

}

/**********************************************************************
 **********************************************************************

 The following routines specify the equation of state.  All routines
  above here should be indpendent of EOS.  If the user wishes
  to use another equation of state, the below functions must be replaced
  by equivalent routines based upon the new EOS.

 **********************************************************************
**********************************************************************/

/**********************************************************************/
/**********************************************************************
  eos_info():

      -- returns with all the EOS-related values needed;
 **********************************************************************/
static CCTK_REAL eos_info(CCTK_REAL W, CCTK_REAL vsq, CCTK_REAL *dpdw, CCTK_REAL *dpdvsq, CCTK_REAL gammaeos, struct LocGlob *lgp)
{
  register CCTK_REAL ftmp,gtmp;

  //fprintf(stdout," !!!!!!!!!!!!!!!!!!!!! \n" );
  //fprintf(stdout," !!! eos_info USED !!! \n" );
  //fprintf(stdout," !!!!!!!!!!!!!!!!!!!!! \n" );

  ftmp = 1. - vsq;
  gtmp = sqrt(ftmp);

  CCTK_REAL gam_m1_o_gam = ((gammaeos-1.)/gammaeos);

  *dpdw =  gam_m1_o_gam * ftmp ;
  *dpdvsq =  gam_m1_o_gam * ( 0.5 * lgp->D/gtmp  -  W ) ;

  return( gam_m1_o_gam * ( W * ftmp  -  lgp->D * gtmp )  );  // p

}

static CCTK_REAL pressure_rho0_eps_eosomni(CCTK_REAL rho,CCTK_REAL epsilon, CCTK_REAL* dpdrho, CCTK_REAL* dpdeps, struct eosomnivars *eosvars)
{

  CCTK_REAL rhopt[1],epspt[1],press[1];
  rhopt[0]=rho;
  epspt[0]=epsilon;

  EOS_Omni_press(eosvars->eoshandle,eosvars->eoskeytemp,eosvars->eosprec,1,
                 &rho,&epsilon,eosvars->eos_temp,
                 eosvars->eos_y_e,press,eosvars->eoskeyerr,eosvars->eosanyerr);

  EOS_Omni_DPressByDRho(eosvars->eoshandle,eosvars->eoskeytemp,eosvars->eosprec,1,
                 &rho,&epsilon,eosvars->eos_temp,
                 eosvars->eos_y_e,dpdrho,eosvars->eoskeyerr,eosvars->eosanyerr);

  EOS_Omni_DPressByDEps(eosvars->eoshandle,eosvars->eoskeytemp,eosvars->eosprec,1,
                 &rho,&epsilon,eosvars->eos_temp,
                 eosvars->eos_y_e,dpdeps,eosvars->eoskeyerr,eosvars->eosanyerr);

  return press[0];

}

/******************************************************************************
             END
 ******************************************************************************/


#undef DEBUG_CON2PRIMM
