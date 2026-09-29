 /*@@
   @file      Spritz_InitZeroVar.F90
   @date      Tue Apr 23 15:44 2019
   @author    Federico Cipolletta
   @desc
   Initialize GRMHD vars to 0.0d0
   @enddesc
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

 /*@@
   @routine    Spritz_Initialize
   @date      Tue Apr 23 15:44 2019
   @author    Federico Cipolletta
   @desc
   Initialize GRMHD vars to 0.0d0
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/

subroutine Spritz_Initialize(CCTK_ARGUMENTS)

  USE Spritz_Scalars

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  call CCTK_INFO('Initializing GRMHD variables to zero...')

  dens = 0.0d0
  rho = 0.0d0
  tau = 0.0d0
  eps = 0.0d0
  press = 0.0d0
  scon = 0.0d0
  vel = 0.0d0
  w_lorentz = 0.0d0
  Bx = 0.0d0
  By = 0.0d0
  Bz = 0.0d0
  dBx = 0.0d0
  dBy = 0.0d0
  dBz = 0.0d0
  Avec_x = 0.0d0
  Avec_y = 0.0d0
  Avec_z = 0.0d0
  psi_mhd = 0.0d0
  Ex = 0.0d0
  Ey = 0.0d0
  Ez = 0.0d0

  densrhs = 0.0d0
  taurhs = 0.0d0
  srhs = 0.0d0
  Bxrhs = 0.0d0
  Byrhs = 0.0d0
  Bzrhs = 0.0d0
  Avec_x_rhs = 0.0d0
  Avec_y_rhs = 0.0d0
  Avec_z_rhs = 0.0d0
  psi_mhdrhs = 0.0d0

  Y_e = 0.0d0
  if(evolve_Y_e.ne.0) then
    Y_e_con = 0.0d0
    Y_e_con_rhs = 0.0d0
  end if
  if(GetTemp_FromThermo.ne.0) then
    temperature = 0.0d0
    entropy = 0.0d0
  end if

  if(compute_DCandMaxima.ne.0)then
    maxrho_global = 0.0d0
    maxrho_global_plus = 0.0d0
    maxrho_global_minus = 0.0d0
  end if

end subroutine Spritz_Initialize
