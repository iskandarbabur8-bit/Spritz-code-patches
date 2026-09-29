// Spritz_RegisterVars.cc
//
// converted from F90 to improve readability and maintainability
//
// Frank Loeffler
//
// MHD variables added by Bruno Giacomazzo
//

#include <cstdio>
#include <string>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

using namespace std;

// Utility functions to register variables with MoL
// Note: We could check for the return value here, but MoL issues a
//       level 0 warning in that case anyway. If that changes in the
//       future, a check can simply be inserted here.
static void register_evolved(string v1, string v2)
{
  MoLRegisterEvolvedGroup(CCTK_GroupIndex(v1.c_str()), CCTK_GroupIndex(v2.c_str()));
}
static void register_constrained(string v1)
{
  MoLRegisterConstrainedGroup(CCTK_GroupIndex(v1.c_str()));
}
static void register_saveandrestore(string v1)
{
  MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex(v1.c_str()));
}

// Main function called by Cactus to register variables with MoL

extern "C" void Spritz_Register(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // We need some aliased functions, so we first check if they are available
  string needed_funs[5] = {"MoLRegisterEvolvedGroup",
                           "MoLRegisterConstrainedGroup",
                           "MoLRegisterSaveAndRestoreGroup",
                           "MoLRegisterEvolved",
                           "MoLRegisterConstrained"};
  for (int i = 0; i < 5; i++)
    if (!CCTK_IsFunctionAliased(needed_funs[i].c_str()))
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "The function \"%s\" has not been aliased!",
                 needed_funs[i].c_str());

  // Now we can set which variables have to be registered as which type with MoL
  register_constrained("HydroBase::rho");
  register_constrained("HydroBase::press");
  register_constrained("HydroBase::eps");
  register_constrained("HydroBase::vel");
  register_constrained("HydroBase::w_lorentz");

  if (CCTK_EQUALS(evolution_method, "Spritz"))
  {
    // dens and scon
    register_evolved("Spritz::dens", "Spritz::densrhs");
    register_evolved("Spritz::scon", "Spritz::srhs");

    // tau
    if (CCTK_EQUALS(Spritz_eos_type, "General"))
      register_evolved("Spritz::tau" , "Spritz::taurhs");
    else if (CCTK_EQUALS(Spritz_eos_type, "Polytype"))
      register_constrained("Spritz::tau");
    else
      CCTK_WARN(0, "Don't recognize the type of EOS!");

    //register_evolved("Spritz::Spritz_mhd_Avec_var", "Spritz::Spritz_mhd_Avec_rhs_var");
    register_evolved("Spritz::Avec_x", "Spritz::Avec_x_rhs");
    register_evolved("Spritz::Avec_y", "Spritz::Avec_y_rhs");
    register_evolved("Spritz::Avec_z", "Spritz::Avec_z_rhs");
    register_constrained("Spritz::Spritz_mhd_vector_var");
    register_constrained("Spritz::Spritz_mhd_evolved_var");
    register_evolved("Spritz::psi_mhd", "Spritz::psi_mhdrhs");

    // lapse, metric, curv (Note that MoL does not register SaveAndRestore variables if they are registered otherwise by other thorns, so the below registration actually happens only when using fixed spacetime.)
    register_saveandrestore("admbase::lapse");
    register_saveandrestore("admbase::metric");
    register_saveandrestore("admbase::curv");

    if (*stress_energy_state)
      {
	register_constrained("TmunuBase::stress_energy_scalar");
	register_constrained("TmunuBase::stress_energy_vector");
	register_constrained("TmunuBase::stress_energy_tensor");
      }

    // shift
    if (!CCTK_EQUALS(initial_shift, "none"))
    {
      register_saveandrestore("admbase::shift");
    }

    // note that we set in pararamcheck
    // evolve_y_e and evolve_temper, but MoLRegister
    // happens before paramcheck...
    if (CCTK_EQUALS(Y_e_evolution_method,"Spritz")) {
      register_constrained("HydroBase::Y_e");
      register_evolved("Spritz::Y_e_con",
                       "Spritz::Y_e_con_rhs");
      }

    if (CCTK_EQUALS(temperature_evolution_method,"Spritz")) {
      register_constrained("HydroBase::temperature");
      register_constrained("HydroBase::entropy");
    }

  }
  else if (CCTK_EQUALS(evolution_method, "none"))
  {
    register_constrained("Spritz::dens");
    register_constrained("Spritz::scon");
    register_constrained("Spritz::tau");
  }
}
