 /*@@
   @file      Spritz_Minima.cc
   @date      Tue Aug 29 18:52:10 2006
   @author
   @desc
        Sets up the scalars used for the atmosphere, after initial data.
   @enddesc
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

//#include "Carpet/Carpet/src/carpet.hh"
#include "carpet.hh"

#ifdef HAVE_CARPET
using namespace Carpet;
#endif

#ifdef __cplusplus
  extern "C" {
#endif

/*#include "EinsteinBase/EOS_Base/src/EOS_Base.h"*/

    /* Scheduled functions */
    void Spritz_Rho_Minima_Setup_Final(CCTK_ARGUMENTS);

    void Spritz_Rho_Minima_Setup_Final_PUGH(CCTK_ARGUMENTS);

#ifdef __cplusplus
  } /* extern "C" */
#endif

 /*@@
   @routine    Spritz_Rho_Minima_Setup_Final
   @date       Tue Aug 29 18:54:05 2006
   @author     Luca Baiotti
   @desc
        After initial data, set up the scalar Spritz_rho_min used for checking the atmosphere.
        This is computed taking into account the actual maximum rest-mass density on the grid.
	Version for Carpet.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory
@@*/

void Spritz_Rho_Minima_Setup_Final(CCTK_ARGUMENTS)
{

#ifdef HAVE_CARPET

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL  max_rho;

  static int flag = true;

  if (flag) // Actually run this reduction only the first time the routine is called.
    {
      if (rho_abs_min > 0.0)
	{
	  *Spritz_rho_min = rho_abs_min;
	}
      else
	{
	  // Go to global mode

	  BEGIN_GLOBAL_MODE (cctkGH) {

            // Find the global maximum of rho

            const int Reduction_Handle = CCTK_ReductionHandle ("maximum");

            const CCTK_INT input_array_variable_indices= {CCTK_VarIndex("HydroBase::rho")};

            const int ierr = CCTK_Reduce(cctkGH,
					 -1, // target processors; -1 -> all
					 Reduction_Handle,
					 1,  // number of output variables
					 CCTK_VARIABLE_REAL,
					 &max_rho,
					 1,  // number of variables to be reduced
					 input_array_variable_indices);

            if (ierr != 0)
	      {
		CCTK_WARN(0, "Failed to compute the global maximum of rho");
	      }

            *Spritz_rho_min = max_rho * rho_rel_min;

            // Go back to local mode
	  } END_GLOBAL_MODE;

	}
      // After this has run once, set the flag so that this does not run again
      flag = false;
    } // end if (flag)
  else
    {
      return;
    }

#define  eps_dummy   0.0
  CCTK_REAL rho_atm, press_atm, eps_atm, dpdrho_atm, dpdeps_atm;

  /*!!!!!!
  ! begin EOS Omni vars*/
  /*CCTK_REAL  Eos_rf_prec;*/
  CCTK_INT n,keytemp,anyerr,keyerr[1];  /* * */
  CCTK_REAL xrho[1],xpress[1],xpressinit[1],xeps[1],xtemp[1],xye[1];
  CCTK_REAL xdpdeps[1],xdpdrho[1],xcs2[1];

  n = 1;
  keytemp = 0;
  keyerr[0] = 0;
  xpress[0] = 0.0;
  xpressinit[0] = eps_dummy;
  xeps[0] = eps_dummy;
  xtemp[0] = 0.0;
  xye[0] = 0.0;
  xcs2[0] = 0.0;
  /*! end EOS Omni vars
  !!!!!!*/

  rho_atm      = *Spritz_rho_min;
  /*press_atm    = EOS_Pressure         ( *Spritz_polytrope_handle, rho_atm, eps_dummy );
  eps_atm      = EOS_SpecificIntEnergy( *Spritz_polytrope_handle, rho_atm, eps_dummy );
  dpdrho_atm   = EOS_DPressByDRho     ( *Spritz_polytrope_handle, rho_atm, eps_dummy );
  dpdeps_atm   = EOS_DPressByDEps     ( *Spritz_polytrope_handle, rho_atm, eps_dummy );*/

  /*!!!!!!*/
  xrho[0] = rho_atm;
  xtemp[0]= Spritz_hot_atmo_temp;
  xye[0]  = Spritz_Y_e_max;
  keytemp = 1;

  if(CCTK_EQUALS(Spritz_eos_table,"nuc_eos")) {
      EOS_Omni_press(*Spritz_eos_handle,keytemp,Eos_rf_prec,n,xrho,xeps,xtemp,xye,xpress,keyerr,&anyerr);
      press_atm = xpress[0];
      eps_atm = xeps[0];
      EOS_Omni_cs2(*Spritz_eos_handle,keytemp,Eos_rf_prec,n,xrho,xeps,xtemp,xye,xcs2,keyerr,&anyerr);
      *Spritz_Cs2_min = xcs2[0];
  }
  else {
      EOS_Omni_press(*Spritz_eos_handle,keytemp,Eos_rf_prec,n,xrho,xeps,xtemp,xye,xpress,keyerr,&anyerr);
      press_atm = xpress[0];
      /* NOTE: EOS_Omni_EpsFromPress used to be called here to set eps_atm,
       * but with keytemp=1 the EOS_Omni_press call above already writes the
       * cold-curve eps directly into xeps -- the EpsFromPress call was pure
       * redundancy (same handle, same keytemp, no new inputs). It also has
       * no case(5) (cold_tabulated) and hard-aborts the run if this EOS is
       * in use, so it's dropped entirely; xeps[0] already holds the right
       * value from the call above. */
      eps_atm = xeps[0];
      EOS_Omni_DPressByDEps(*Spritz_eos_handle,keytemp,Eos_rf_prec,n,xrho,xeps,xtemp,xye,xdpdeps,keyerr,&anyerr);
      dpdeps_atm = xdpdeps[0];
      EOS_Omni_DPressByDRho(*Spritz_eos_handle,keytemp,Eos_rf_prec,n,xrho,xeps,xtemp,xye,xdpdrho,keyerr,&anyerr);
      dpdrho_atm = xdpdrho[0];
      /*!!!!!!*/
      *Spritz_Cs2_min = ( dpdrho_atm + press_atm * dpdeps_atm / ( rho_atm*rho_atm ) )
        / ( 1.0 + eps_atm + press_atm / rho_atm );
  }
  if ( *Spritz_Cs2_min < 0.0 ) CCTK_WARN( CCTK_WARN_ABORT, "Spritz_Cs2_min is negative.");

#undef  eps_dummy

#endif

  //Debug stuff
  // char warnline;

  //CCTK_VInfo (CCTK_THORNSTRING, "STEP 2: compute rho max; rho min: %13.12e \n",*Spritz_rho_min);

 //  printf(warnline,"STEP 2: compute rho max; rho min: %13.12e \n",*Spritz_rho_min);
 //CCTK_WARN(1,warnline);

  //  printf(warnline,"STEP 3: recompute ID with new atmosphere; rho min: ', Spritz_rho_min, 'reflev: ',%i), Spritz_rho_min, Spritz_reflevel);



}


 /*@@
   @routine    Spritz_Rho_Minima_Setup_Final_PUGH
   @date       Tue Aug 29 18:57:42 2006
   @author     Luca Baiotti
   @desc
         As above, but for PUGH.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory
@@*/

void Spritz_Rho_Minima_Setup_Final_PUGH(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL  max_rho;

  if (rho_abs_min > 0.0)
    {
      *Spritz_rho_min = rho_abs_min;
    }
  else
    {
      // Find the global maximum of rho

      const int Reduction_Handle = CCTK_ReductionHandle("maximum");

      const CCTK_INT input_array_variable_indices={CCTK_VarIndex("HydroBase::rho")};

      const int ierr = CCTK_Reduce(cctkGH,
			 -1, // target processors; -1 -> all
			 Reduction_Handle,
			 1,  // number of output variables
			 CCTK_VARIABLE_REAL,
			 &max_rho,
			 1,  // number of variables to be reduced
			 input_array_variable_indices);

      if (ierr != 0)
	{
	  CCTK_WARN(0, "Failed to compute the global maximum of rho");
	}

      *Spritz_rho_min = max_rho * rho_rel_min;
    }

#define  eps_dummy   0.0
  CCTK_REAL rho_atm, press_atm, eps_atm, dpdrho_atm, dpdeps_atm;

  /*!!!!!!
  ! begin EOS Omni vars*/
  /*CCTK_REAL  Eos_rf_prec;*/
  CCTK_INT n,keytemp,anyerr,keyerr[1]; /* * */
  CCTK_REAL xrho[1],xpress[1],xpressinit[1],xeps[1],xtemp[1],xye[1];
  CCTK_REAL xdpdeps[1],xdpdrho[1],xcs2[1];

  n = 1;
  keytemp = 0;
  keyerr[0] = 0;
  xpress[0] = 0.0;
  xpressinit[0] = eps_dummy;
  xeps[0] = eps_dummy;
  xtemp[0] = 0.0;
  xye[0] = 0.0;
  /*! end EOS Omni vars
  !!!!!!*/

  rho_atm      = *Spritz_rho_min;
  /*press_atm    = EOS_Pressure         ( *Spritz_polytrope_handle, rho_atm, eps_dummy );
  eps_atm      = EOS_SpecificIntEnergy( *Spritz_polytrope_handle, rho_atm, eps_dummy );
  dpdrho_atm   = EOS_DPressByDRho     ( *Spritz_polytrope_handle, rho_atm, eps_dummy );
  dpdeps_atm   = EOS_DPressByDEps     ( *Spritz_polytrope_handle, rho_atm, eps_dummy );*/

  /*!!!!!!*/
  xrho[0] = rho_atm;
  xtemp[0]= Spritz_hot_atmo_temp;
  xye[0]  = Spritz_Y_e_max;
  keytemp = 1;
  EOS_Omni_press(*Spritz_eos_handle,keytemp,Eos_rf_prec,n,xrho,xeps,xtemp,xye,xpress,keyerr,&anyerr);
  press_atm = xpress[0];
  /* NOTE: see the identical case above -- keytemp=1 already makes the
   * EOS_Omni_press call write the cold-curve eps into xeps directly, so the
   * EpsFromPress call that used to sit here was pure redundancy (and has no
   * case(5)/cold_tabulated support, hard-aborting the run). Dropped. */
  eps_atm = xeps[0];
  EOS_Omni_DPressByDEps(*Spritz_eos_handle,keytemp,Eos_rf_prec,n,xrho,xeps,xtemp,xye,xdpdeps,keyerr,&anyerr);
  dpdeps_atm = xdpdeps[0];
  EOS_Omni_DPressByDRho(*Spritz_eos_handle,keytemp,Eos_rf_prec,n,xrho,xeps,xtemp,xye,xdpdrho,keyerr,&anyerr);
  dpdrho_atm = xdpdrho[0];
  EOS_Omni_DPressByDRho(*Spritz_eos_handle,keytemp,Eos_rf_prec,n,xrho,xeps,xtemp,xye,xcs2,keyerr,&anyerr);
  *Spritz_Cs2_min = xcs2[0];

  if ( *Spritz_Cs2_min < 0.0 )
    { CCTK_WARN( CCTK_WARN_ABORT, "Spritz_Cs2_min is negative."); }
#undef  eps_dummy

}
