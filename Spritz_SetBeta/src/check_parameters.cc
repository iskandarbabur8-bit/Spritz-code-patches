#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <fstream>
#include <string.h>

extern "C"

bool checkfilepresence(std::string filename) {
    std::ifstream file(filename, std::ios_base::in);
    if(file.is_open()) {
        file.close();
        return true;
    }
    else return false;
}

void Spritz_SetBeta_check_parameters (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS (initial_temperature, "Spritz_SetBeta") or
      CCTK_EQUALS (initial_Y_e,         "Spritz_SetBeta"))
  {
      if (not ( CCTK_EQUALS (initial_temperature, "Spritz_SetBeta") or
                CCTK_EQUALS (initial_Y_e,         "Spritz_SetBeta") or
                not CCTK_EQUALS (beta_datafile,   "")))
      {
        CCTK_PARAMWARN("The parameters HydroBase::initial_temperature and HydroBase::initial_Y_e must be set both to the value \"Spritz_SetBeta\" or \"none\".Parameter Spritz_SetBeta::Beta_DataFile must be set to the beta-equilibrated file path.");
      }

      if (CCTK_EQUALS (beta_df_format,"Compose")) {
          //check existence of EoS file that are to be used.
          std::string filepath = beta_datafile;
          if(!checkfilepresence(filepath+".t")) CCTK_ERROR("Error in reading EoS files. eos.t may not exist, reexecute the python script.");
          if(!checkfilepresence(filepath+".nb")) CCTK_ERROR("Error in reading EoS files. eos.nb may not exist, reexecute the python script.");
          if(!checkfilepresence(filepath+".yqb")) CCTK_ERROR("Error in reading EoS files. eos.yqb may not exist, reexecute the python script.");
          if(!checkfilepresence(filepath+".thermo")) CCTK_ERROR("Error in reading EoS files. eos.thermo may not exist, reexecute the python script.");
      }
      else if (CCTK_EQUALS (beta_df_format,"Lorene")) {
          //check existence of EoS file that are to be used.
          std::string filepath = beta_datafile;
          if(!checkfilepresence(filepath+"LoreneInput0.01.txt")) CCTK_ERROR("Error in reading EoS files. LoreneInput*.txt may not exist.");
          if(!checkfilepresence(filepath+"ye.txt")) CCTK_ERROR("Error in reading EoS files. ye.txt may not exist.");
      }
      else {
		  CCTK_ERROR("Beta Data file's format not recognized.");
	  }
  }
}
