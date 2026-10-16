/*@@
 * @file		Prim2ConInitial
 * @date		Wed Jun 30 14:55:48 CEST 2010
 * @author		David Radice
 * @desc
 * Setup the atmosphere and compute the conserved variables after the creation
 * of the initial data
 * @enddesc
  @@*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

/*#include "EOS_Base.h"*/

/* Spritz_coldtab_eos.c: centralized matched-polytrope replacement for
 * EOS_Omni_press below the cold_tabulated table's rho_min (real
 * locally-measured Gamma, not EOS_Omni's own hardcoded Gamma=2). Falls
 * through to the real EOS_Omni_press call unchanged otherwise. */
void Spritz_coldtab_press(
    CCTK_INT handle, CCTK_INT keytemp, CCTK_REAL rf_precision, CCTK_INT npoints,
    CCTK_REAL *rho, CCTK_REAL *eps, CCTK_REAL *temp, CCTK_REAL *ye,
    CCTK_REAL *press, CCTK_INT *keyerr, CCTK_INT *anyerr,
    CCTK_REAL coldtab_rhomin, CCTK_REAL coldtab_press_lo,
    CCTK_REAL coldtab_eps_lo, CCTK_REAL coldtab_gamma_lo,
    CCTK_REAL coldtab_gammath, CCTK_REAL coldtab_thfac);

#define MIN(X,Y) ((X)<(Y)?(X):(Y))

#define SQ(X) ((X)*(X))

#define WH_LOOP3(NAME,I,SI,EI,J,SJ,EJ,K,SK,EK)								\
	for(I = (SI); I < (EI); ++I)											\
		for(J = (SJ); J < (EJ); ++J)										\
			for(K = (SK); K < (EK); ++K)

#define WH_ENDLOOP3(NAME)

#define MHDNVARS 32
#define HDNVARS 19

/*@@
 * @routine _M_GetLevelVI
 * @desc
 * Get the pointer to the data of the wanted variable at the given time level
 * if the timelevel is not filled in then the data will be extrapolated back
 * in time with a first order extrapolation.
 * @enddesc
  @@*/
CCTK_REAL * _M_GetLevelVI(
		cGH const * restrict const cctkGH,
		int const level,
		int const varidx,
		int * restrict const extrapolation_flag
		) {
	DECLARE_CCTK_PARAMETERS

	assert(varidx >= 0);

	CCTK_REAL *  var;
	CCTK_REAL ** varh;
	int i, j, k, l;
	int ijk;

	CCTK_REAL * vijk;

	int maxlevel = CCTK_ActiveTimeLevelsVI(cctkGH, varidx);

	char const * varname;

	if(level < maxlevel || maxlevel < 3) {
		*extrapolation_flag = 0;
		return CCTK_VarDataPtrI(cctkGH, MIN(level,maxlevel), varidx);
	}
	else {
		varname = CCTK_VarName(varidx);
		CCTK_VInfo(CCTK_THORNSTRING, "Extrapolating back in time '%s'. "
				"If you don't want this to happen you should provide at least "
				"%d timelevels for '%s'.",varname,level+1,varname);

		*extrapolation_flag = 1;

		var = malloc(cctkGH->cctk_lsh[0]*cctkGH->cctk_lsh[1]*
				cctkGH->cctk_lsh[2]*sizeof(CCTK_REAL));
		if(!var) {
			CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
				   "Could not allocate temporary timelevel for '%s'",varname);
		}

		varh = malloc(maxlevel*sizeof(*varh));
		assert(varh);
		for(l = 0; l < maxlevel; ++l) {
			varh[l] = CCTK_VarDataPtrI(cctkGH, l, varidx);
		}
#pragma omp parallel private(i,j,k,l,ijk,vijk)
		{
			vijk = malloc(level*sizeof(CCTK_REAL));
			assert(vijk);

			WH_LOOP3(extrapolate, k,0,cctkGH->cctk_lsh[2],
					j,0,cctkGH->cctk_lsh[1], i,0,cctkGH->cctk_lsh[0]) {
				ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
				for(l = 0; l < maxlevel; ++l) {
					vijk[l] = varh[l][ijk];
				}
				for(l = maxlevel; l < level; ++l) {
					vijk[l] = 2*vijk[l-1] - vijk[l-2];
				}
				var[ijk] = vijk[level-1];
			} WH_ENDLOOP3(extrapolate);

			free(vijk);
		}
		free(varh);

		return var;
	}
}

/*@@
 * @routine _M_GetLevelVN
 * @desc
 * Get the pointer to the data of the wanted variable at the given time level
 * if the timelevel is not filled in then the data will be extrapolated back
 * in time with a first order extrapolation.
 * @enddesc
  @@*/
CCTK_REAL * _M_GetLevelVN(
		cGH const * restrict const cctkGH,
		int const level,
		char const * restrict const varname,
		int * restrict const extrapolation_flag
		) {
	return _M_GetLevelVI(cctkGH, level, CCTK_VarIndex(varname),
			extrapolation_flag);
}

/*@@
 * @routine _M_Prim2ConInitial
 * @desc
 * Setup the atmosphere and compute the conserved quantities on the given
 * timelevel.
 * @enddesc
  @@*/
void _M_Prim2ConInitial(
		cGH const * restrict const cctkGH,
		CCTK_REAL * restrict const _M_rho,			/* 0  */
	   	CCTK_REAL * restrict const _M_velx,			/* 1  */
	   	CCTK_REAL * restrict const _M_vely,			/* 2  */
		CCTK_REAL * restrict const _M_velz,			/* 3  */
		CCTK_REAL * restrict const _M_eps,			/* 4  */
		CCTK_REAL * restrict const _M_press,		/* 5  */
	   	CCTK_REAL * restrict const _M_dens,			/* 6  */
	   	CCTK_REAL * restrict const _M_sx,			/* 7  */
		CCTK_REAL * restrict const _M_sy,			/* 8  */
	   	CCTK_REAL * restrict const _M_sz,			/* 9  */
	   	CCTK_REAL * restrict const _M_tau,			/* 10 */
		CCTK_REAL * restrict const _M_w_lorentz,	/* 11 */
	   	CCTK_REAL * restrict const _M_gxx,			/* 12 */
	   	CCTK_REAL * restrict const _M_gxy,			/* 13 */
		CCTK_REAL * restrict const _M_gxz,			/* 14 */
	   	CCTK_REAL * restrict const _M_gyy,			/* 15 */
	   	CCTK_REAL * restrict const _M_gyz,			/* 16 */
		CCTK_REAL * restrict const _M_gzz,			/* 17 */
		CCTK_REAL * restrict const _M_psi,			/* 18 */
		CCTK_REAL * restrict const _M_alp,			/* 19 */
		CCTK_REAL * restrict const _M_betax,		/* 20 */
		CCTK_REAL * restrict const _M_betay,		/* 21 */
		CCTK_REAL * restrict const _M_betaz,		/* 22 */
		CCTK_REAL * restrict const _M_Bx,			/* 23 */
		CCTK_REAL * restrict const _M_By,			/* 24 */
		CCTK_REAL * restrict const _M_Bz,			/* 25 */
		CCTK_REAL * restrict const _M_dBx,			/* 26 */
		CCTK_REAL * restrict const _M_dBy,			/* 27 */
		CCTK_REAL * restrict const _M_dBz,			/* 28 */
		CCTK_REAL * restrict const _M_Y_e,    /* 29 */
		CCTK_REAL * restrict const _M_Y_e_con,    /* 30 */
		CCTK_REAL * restrict const _M_temp     /* 31 */
		) {
	DECLARE_CCTK_ARGUMENTS
	DECLARE_CCTK_PARAMETERS

	int i, j, k;
	int ijk;

	CCTK_REAL psi4pt;
	CCTK_REAL det;

        /*!!!!!!
        ! begin EOS Omni vars*/
        CCTK_INT  anyerr;
        /*CCTK_INT xanyerr;*/
        CCTK_INT nn,xkeytemp,keyerr[1];
        CCTK_REAL xrho[1],xpress[1],xeps1[1],xeps[1],xtemp[1],xye[1];
        nn = 1;
        xkeytemp = 0;
        anyerr = 0; /*xanyerr;*/
        keyerr[0] = 0;
        xeps[0] = 0.0;
        xtemp[0] = *_M_temp;
        xye[0] = *_M_Y_e;
        /*! end EOS Omni vars
        !!!!!!*/

#pragma omp parallel private(i,j,k,ijk,det,psi4pt,anyerr,xkeytemp,keyerr,xrho,xpress,xeps1,xeps,xtemp,xye)
	{
		WH_LOOP3(setupatm, k,0,cctkGH->cctk_lsh[2], j,0,cctkGH->cctk_lsh[1],
				i,0,cctkGH->cctk_lsh[0]) {
			ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

			/* check if we use a conformal metric (i.e. the value of the integer variable conformal_state is larger than 0) */
			if(*conformal_state>0 && _M_psi) {
			  /* if (i==0 && j==0 && k==0) printf("conformal_state=%i \n",*conformal_state); */ /* debug line */
			  psi4pt = SQ(SQ(_M_psi[ijk]));
			  SpatialDet(psi4pt*_M_gxx[ijk], psi4pt*_M_gxy[ijk],
				     psi4pt*_M_gxz[ijk], psi4pt*_M_gyy[ijk],
				     psi4pt*_M_gyz[ijk], psi4pt*_M_gzz[ijk],&det);
			}
			else {
			  /* if (i==0 && j==0 && k==0) printf("NOT USED, but conformal_state=%i \n",*conformal_state);*/ /* debug line */
			  SpatialDet(_M_gxx[ijk], _M_gxy[ijk], _M_gxz[ijk], _M_gyy[ijk],
				     _M_gyz[ijk], _M_gzz[ijk],&det);
			}

			if(_M_rho[ijk] < *Spritz_rho_min) {
				_M_rho[ijk]  = *Spritz_rho_min;
				_M_velx[ijk] = 0.0;
				_M_vely[ijk] = 0.0;
				_M_velz[ijk] = 0.0;
				_M_temp[ijk] = Spritz_hot_atmo_temp;

				/* Routed through Spritz_coldtab_press (Spritz_coldtab_eos.c),
				 * the same centralized matched-polytrope wrapper used by
				 * con2prim/prim2con/reconstruction: below the cold_tabulated
				 * table's rho_min it uses the table's real local Gamma
				 * (calibrated once by Spritz_ColdTabCalibration,
				 * Spritz_Minima.cc) instead of EOS_Omni's own hardcoded
				 * Gamma=2 -- so the initial-data atmosphere and any later
				 * con2prim-triggered atmosphere reset agree on eps/press at
				 * the same density. For any other EOS table, or any density
				 * at/above rho_min, falls through to the unchanged real
				 * EOS_Omni_press call. */
				xtemp[0] = _M_temp[ijk];
				xrho[0]  = _M_rho[ijk];
				xye[0]   = _M_Y_e[ijk];
				xkeytemp = 1; // in case we are in atmosphere, we already know the temperature. Therefore we want to keep it fixed.
				Spritz_coldtab_press(*Spritz_eos_handle,xkeytemp,Eos_rf_prec,nn,xrho,xeps,xtemp,xye,xpress,keyerr,&anyerr,
					*Spritz_coldtab_rhomin,*Spritz_coldtab_press_lo,
					*Spritz_coldtab_eps_lo,*Spritz_coldtab_gamma_lo,
					*Spritz_coldtab_gammath,*Spritz_coldtab_thfac);
				_M_eps[ijk]   = xeps[0];
				_M_press[ijk] = xpress[0];
			}
			else {
				if(CCTK_Equals(Spritz_eos_type,"General")) {
					if(CCTK_Equals(Spritz_eos_table,"nuc_eos")) {
						Prim2ConGen_Hot_MHD(*Spritz_eos_handle, *Constrain_Temp_now, *Spritz_reflevel,
							  0.0, 0.0, 0.0, _M_gxx[ijk], _M_gxy[ijk],
							  _M_gxz[ijk], _M_gyy[ijk], _M_gyz[ijk], _M_gzz[ijk],
							  det,&_M_dens[ijk],&_M_sx[ijk],&_M_sy[ijk],&_M_sz[ijk],
							  &_M_tau[ijk], _M_rho[ijk], _M_velx[ijk], _M_vely[ijk],
							  _M_velz[ijk], &_M_eps[ijk],&_M_press[ijk],
							  &_M_w_lorentz[ijk], _M_Bx[ijk], _M_By[ijk], _M_Bz[ijk],
							  &_M_dBx[ijk],&_M_dBy[ijk],&_M_dBz[ijk], _M_alp[ijk],
							  _M_betax[ijk], _M_betay[ijk], _M_betaz[ijk],
								&_M_temp[ijk], _M_Y_e[ijk]);
					}
					else {
					  Prim2ConGen_MHD(*Spritz_eos_handle, _M_gxx[ijk], _M_gxy[ijk],
							  _M_gxz[ijk], _M_gyy[ijk], _M_gyz[ijk], _M_gzz[ijk],
							  det,&_M_dens[ijk],&_M_sx[ijk],&_M_sy[ijk],&_M_sz[ijk],
							  &_M_tau[ijk], _M_rho[ijk], _M_velx[ijk], _M_vely[ijk],
							  _M_velz[ijk], _M_eps[ijk],&_M_press[ijk],
							  &_M_w_lorentz[ijk], _M_Bx[ijk], _M_By[ijk], _M_Bz[ijk],
							  &_M_dBx[ijk],&_M_dBy[ijk],&_M_dBz[ijk], _M_alp[ijk],
							  _M_betax[ijk], _M_betay[ijk], _M_betaz[ijk],
							  *Spritz_coldtab_rhomin, *Spritz_coldtab_press_lo,
							  *Spritz_coldtab_eps_lo, *Spritz_coldtab_gamma_lo,
							  *Spritz_coldtab_gammath, *Spritz_coldtab_thfac);
					}
				}
				else {
				  Prim2ConPoly_MHD(*Spritz_eos_handle, _M_gxx[ijk], _M_gxy[ijk],
						   _M_gxz[ijk], _M_gyy[ijk], _M_gyz[ijk], _M_gzz[ijk],
						   det,&_M_dens[ijk],&_M_sx[ijk],&_M_sy[ijk],&_M_sz[ijk],
						   &_M_tau[ijk], _M_rho[ijk], _M_velx[ijk], _M_vely[ijk],
						   _M_velz[ijk], &_M_eps[ijk],&_M_press[ijk],
						   &_M_w_lorentz[ijk], _M_Bx[ijk], _M_By[ijk], _M_Bz[ijk],
						   &_M_dBx[ijk],&_M_dBy[ijk],&_M_dBz[ijk], _M_alp[ijk],
						   _M_betax[ijk], _M_betay[ijk], _M_betaz[ijk]);
				}
			}
			if(*evolve_Y_e==1) _M_Y_e_con[ijk] = _M_Y_e[ijk] * _M_dens[ijk];
		} WH_ENDLOOP3(setupatm);
	}
}

/*@@
 * @routine Spritz_Prim2ConInitial
 * @desc
 * Setup the atmosphere and compute the conserved quantities on all the
 * timelevels.
 * @enddesc
  @@*/
void Spritz_Prim2ConInitial(CCTK_ARGUMENTS) {
	DECLARE_CCTK_ARGUMENTS
	DECLARE_CCTK_PARAMETERS

	// CCTK_INFO("Initializing the atmosphere and the conserved variables");

	int extrapolation_flag[MHDNVARS];
	int vardb[MHDNVARS] = {
			CCTK_VarIndex("HydroBase::rho"),		/* 0  */
			CCTK_VarIndex("HydroBase::vel[0]"),		/* 1  */
			CCTK_VarIndex("HydroBase::vel[1]"),		/* 2  */
			CCTK_VarIndex("HydroBase::vel[2]"),		/* 3  */
			CCTK_VarIndex("HydroBase::eps"),		/* 4  */
			CCTK_VarIndex("HydroBase::press"),		/* 5  */
			CCTK_VarIndex("Spritz::dens"),			/* 6  */
			CCTK_VarIndex("Spritz::scon[0]"),		/* 7  */
			CCTK_VarIndex("Spritz::scon[1]"),		/* 8  */
			CCTK_VarIndex("Spritz::scon[2]"),		/* 9  */
			CCTK_VarIndex("Spritz::tau"),			/* 10 */
			CCTK_VarIndex("HydroBase::w_lorentz"),		/* 11 */
			CCTK_VarIndex("ADMBase::gxx"),			/* 12 */
			CCTK_VarIndex("ADMBase::gxy"),			/* 13 */
			CCTK_VarIndex("ADMBase::gxz"),			/* 14 */
			CCTK_VarIndex("ADMBase::gyy"),			/* 15 */
			CCTK_VarIndex("ADMBase::gyz"),			/* 16 */
			CCTK_VarIndex("ADMBase::gzz"),			/* 17 */
			CCTK_VarIndex("StaticConformal::psi"),	/* 18 */
			CCTK_VarIndex("ADMBase::alp"),			/* 19 */
			CCTK_VarIndex("ADMBase::betax"),		/* 20 */
			CCTK_VarIndex("ADMBase::betay"),		/* 21 */
			CCTK_VarIndex("ADMBase::betaz"),		/* 22 */
			CCTK_VarIndex("Spritz::Bx"),			/* 23 */
			CCTK_VarIndex("Spritz::By"),			/* 24 */
			CCTK_VarIndex("Spritz::Bz"),			/* 25 */
			CCTK_VarIndex("Spritz::dBx"),			/* 26 */
			CCTK_VarIndex("Spritz::dBy"),			/* 27 */
			CCTK_VarIndex("Spritz::dBz"),			/* 28 */
			CCTK_VarIndex("Hydrobase::Y_e"),			/* 29 */
			CCTK_VarIndex("Spritz::Y_e_con"),			/* 30 */
			CCTK_VarIndex("Hydrobase::temperature") /* 31 */
		};
	CCTK_REAL * varptr[MHDNVARS];
	int idx;
	int nvars = MHDNVARS;

	int levl;
	int nlevels;
	if(CCTK_Equals(initial_data_setup_method,"init_all_levels")) {
		nlevels = timelevels;
	}
	else if(CCTK_Equals(initial_data_setup_method,"init_two_levels")) {
		nlevels = MIN(2,timelevels);
	}
	else {
		nlevels = 1;
	}

	/*
	 * We need to locally compute Spritz_rho_min as Spritz_rho_central may not
	 * be the maximum value of the density
	 */
	if(initial_rho_abs_min > 0 || rho_abs_min > 0) {
		*Spritz_rho_min = (initial_rho_abs_min > 0) ?
			initial_rho_abs_min : rho_abs_min;
	}
	else {
		CCTK_REAL max_rho;
		int const reduction_handle = CCTK_ReductionHandle("maximum");
		int const local_vardb = {CCTK_VarIndex("HydroBase::rho")};
		int const ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
				CCTK_VARIABLE_REAL, &max_rho, 1, local_vardb);
		assert(!ierr);

		*Spritz_rho_min = ((initial_rho_rel_min > 0) ?
			initial_rho_rel_min : rho_rel_min) * max_rho;
		if(initial_atmosphere_factor > 0) {
			*Spritz_rho_min *= initial_atmosphere_factor;
		}
	}

	for(levl = 0; levl < nlevels; ++levl) {
		for(idx = 0; idx < nvars; ++idx) {
			varptr[idx] = _M_GetLevelVI(cctkGH, levl,
					vardb[idx], &extrapolation_flag[idx]);
			/* psi storage could be disabled */
			assert(idx == 18 || varptr[idx]);
		}
		for( ; idx < MHDNVARS; ++idx) {
			extrapolation_flag[idx] = 0;
			varptr[idx] = NULL;
		}

		_M_Prim2ConInitial(cctkGH, varptr[0], varptr[1], varptr[2], varptr[3],
			   varptr[4], varptr[5], varptr[6], varptr[7], varptr[8],
			   varptr[9], varptr[10], varptr[11], varptr[12], varptr[13],
			   varptr[14], varptr[15], varptr[16], varptr[17], varptr[18],
			   varptr[19], varptr[20], varptr[21], varptr[22], varptr[23],
			   varptr[24], varptr[25], varptr[26], varptr[27], varptr[28],
			   varptr[29], varptr[30], varptr[31]);

		for(idx = 0; idx < MHDNVARS; ++idx) {
		if(extrapolation_flag[idx]) {
				free(varptr[idx]);
			}
		}
	}
}
