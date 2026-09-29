/*@@
  @file      SetBetaEquilibrium.cc
  @date
  @author    Lorenzo Sala
  @desc
  Sets up initial temperature and Y_e data reading from file.
  @enddesc
  @@*/

#include <cassert>
#include <iostream>
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
#include <sstream>

// headers Cactus
#include "cctk_Functions.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

using namespace std;

extern "C"

double linearinterp(vector<double> xarray, int xindex, double thisx, vector<double> yarray) {
    // 1st order interpolation
    double x1, x2, y1, y2, dx;
    y1 = yarray[xindex];
    if(xindex != xarray.size() - 1) {
        x1 = xarray[xindex];
        x2 = xarray[xindex+1];
        y2 = yarray[xindex+1];
    } else {
        x1 = xarray[xindex-1];
        x2 = xarray[xindex];
        y2 = yarray[xindex];
    }
    dx = thisx - x1;
    return y1 + dx * (y2 - y1)/(x2 - x1);
}

//lagrangian interpolation routine
double lagrinterp(double thisx, double xpoints[], double ypoints[], int n) {
	double result = 0; // Initialize result
	for (int i=0; i<n; i++) {
		double term = ypoints[i];
		for (int j=0;j<n;j++) {
			if (j!=i)
				term = term*(thisx - xpoints[j]) / (xpoints[i] - xpoints[j]);
		}
		result += term;
	}
	return result;
}

// linear interpolation at borders, lagrangian 5-points inside borders.
double interp1D(vector<double> xarray, int xindex, double thisx, vector<double> yarray) {
    double result = 0.;
    if(xindex < 4 || xindex >= xarray.size()-4) {
        result = linearinterp(xarray, xindex, thisx, yarray);
    }//linear interpolation between the first and last 3 points
    else { //lagrange intepolation
        double xpoints[5] = {xarray[xindex-2],xarray[xindex-1],xarray[xindex],xarray[xindex+1],xarray[xindex+2]};
        double ypoints[5] = {yarray[xindex-2],yarray[xindex-1],yarray[xindex],yarray[xindex+1],yarray[xindex+2]};
        result = lagrinterp(thisx,xpoints,ypoints,5);
        // if(xindex > 213) { //DEBUG STUFF
        //     cout << xindex << " " << xpoints[0] << " " << xpoints[1] << " " << xpoints[2] << " " << xpoints[3] << " " << xpoints[4] << endl;
        //     cout << xindex << " " << ypoints[0] << " " << ypoints[1] << " " << ypoints[2] << " " << ypoints[3] << " " << ypoints[4] << endl;
        //     cout << thisx << " " << result << endl;
        // }
    }
    return result;
}

void Spritz_SetBetaEquilibrium (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO ("Setting up initial Beta Equilibrated data for Compose format.");

  string readme_filename = beta_datafile;
  readme_filename = readme_filename.erase(readme_filename.length()-3) + "README";
  ifstream readme_file(readme_filename.c_str());
  bool ifreadme = readme_file.is_open();
  if(ifreadme){
    readme_file.ignore(1000,'\n');
    string infoline; getline(readme_file,infoline);
    CCTK_VINFO("Betaslice info: %s\n", infoline.c_str());
  }
  readme_file.close();


  string temp_filename = beta_datafile;
  temp_filename += ".t"; //reading temperature (if it is an entropyslice, this won't be used)
  ifstream temp_file(temp_filename.c_str());
    temp_file.ignore(1000,'\n'); //get rid of first two lines
    temp_file.ignore(1000,'\n');
    double tempslice; temp_file >> tempslice;
  temp_file.close();

  // read electron fraction data, as function of nb.
    string ye_filename = beta_datafile;
    ye_filename += ".yqb"; //reading electron fraction at beta equilibrium
    ifstream ye_file(ye_filename.c_str());
    vector<double> ye_beta;
    while(1) { //fill array yearray
        double tmp_ye;
        ye_file >> tmp_ye;
        if(ye_file.eof()) break;
        ye_beta.push_back(tmp_ye);
    }
    ye_file.close();

    // read thermo quantities before nb, 'cause we need m_neutron_MeV
    double m_neutron_MeV;
    const double c_siinv2 = 1.1126500560536185E-17; //2.99792458E+8 ;	 ///< Velocity of light [m/s] //same as in Lorene
    const double mev_si = 1.602176565E-13 ;   ///< One MeV [J] //same as in Lorene

    string thermo_filename = beta_datafile;
    thermo_filename += ".thermo"; //reading thermo qties in beta-equilibrium
    ifstream betaeq_data_file(thermo_filename.c_str());
    vector<double> betaentropy, energies;
    double tmpentropy, tmpenergy, dummy;
    int tempidx, rhoidx, yeidx; //not used
    betaeq_data_file >> m_neutron_MeV >> dummy; //get neutron mass
    betaeq_data_file.ignore(1000,'\n'); //get rid of remaining part of first line
    while(1) { //fill array ye_index
      betaeq_data_file >> tempidx >> rhoidx >> yeidx >> dummy >> tmpentropy >> dummy >> dummy >> dummy >> dummy >> tmpenergy;
      if(betaeq_data_file.eof()) break;
      betaentropy.push_back(tmpentropy);
      energies.push_back(tmpenergy);
    }
    betaeq_data_file.close();

  // set up rho, as function of nb and electron fraction previously found
    double geo_to_si_rho = 6.17714E17;
    double nbtorho_convfac = log10(m_neutron_MeV * mev_si * c_siinv2 / geo_to_si_rho) + 42.0;

    string rho_filename = beta_datafile;
    rho_filename += ".nb"; //reading baryon density
    ifstream rho_file(rho_filename.c_str());
      rho_file.ignore(1000,'\n'); //get rid of first two lines
      rho_file.ignore(1000,'\n');
      vector<double> cactuslogrho;
      while(1) { //fill array cactuslogrho
        double nb;
        rho_file >> nb;
        if(rho_file.eof()) break;
        cactuslogrho.push_back( log10(nb) + nbtorho_convfac);
      }
      rho_file.close();
      if(ye_beta.size() != cactuslogrho.size()) CCTK_WARN(0,"Error in EoS tables: ye size not eq to rho size.");
      if(energies.size() != cactuslogrho.size()) CCTK_WARN(0,"Error in EoS tables: energy size not eq to rho size.");
      if(betaentropy.size() != cactuslogrho.size()) CCTK_WARN(0,"Error in EoS tables: entropy size not eq to rho size.");

  int const local_grid_size = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

  vector<double> temperatures;
  string cs_filename = beta_datafile;
  cs_filename += ".t_s"; //if eos is an entropyslice, read corresponding temperatures (depending on density)
  ifstream constentropy_temp_file(cs_filename.c_str());
  // eos is entropyslice iff eos.t_s exists.
  bool entropyslice = constentropy_temp_file.is_open();
  if(entropyslice){
      CCTK_INFO ("Temperature will be set at constant entropy.");
      while(1) {
          double tmp;
          constentropy_temp_file >> tmp;
          if(constentropy_temp_file.eof()) break;
          temperatures.push_back(tmp);
      }
    if(temperatures.size() != cactuslogrho.size()) CCTK_WARN(0,"Error in EoS tables: temp size not eq to rho size.");
  }
  else CCTK_VINFO("Temperature[MeV] will be set to %f", tempslice);
  constentropy_temp_file.close();

  // set up Ye_hot_atmo value, choosing atmo density
    int rhoindex = 0; //first of all, get closest index to the left
    for (int i = 0; i < cactuslogrho.size() && cactuslogrho[i] < log10(rho_abs_min); ++i) rhoindex = i;
    double thislogrho = log10(rho_abs_min);
    double hot_atmo_Y_e = interp1D(cactuslogrho,rhoindex,thislogrho,ye_beta);
    if(CCTK_IsThornActive("Spritz")) {
        int chk = 1;
        chk = CCTK_ParameterSet ("Spritz_hot_atmo_Y_e","Spritz",to_string(hot_atmo_Y_e).c_str());
        if(chk == 0) CCTK_VINFO("Spritz_hot_atmo_Y_e has been set to %s", to_string(hot_atmo_Y_e).c_str());
        else         CCTK_INFO("Error in setting Spritz_hot_atmo_Y_e");
        std::ostringstream streamObj;
        streamObj << pow(10.0,cactuslogrho[0]) * 1.01; //just to be sure that it is greater than the minimum in table
        chk = CCTK_ParameterSet ("rho_abs_eos_min",    "Spritz",streamObj.str().c_str());
        if(chk == 0) CCTK_VINFO("rho_abs_eos_min has been set to %s", streamObj.str().c_str());
        else         CCTK_INFO("Error in setting rho_abs_eos_min");
        if(rho_abs_min < rho_abs_eos_min) {
            CCTK_VINFO("rho_abs_min is %f, whereas the minimum in table is %f.", rho_abs_min, rho_abs_eos_min);
            CCTK_PARAMWARN("rho_abs_min must be greater than the minimum in table.");
        }
    }
    else if(CCTK_IsThornActive("GRHydro")) {
        CCTK_ParameterSet ("GRHydro_hot_atmo_Y_e","GRHydro",to_string(hot_atmo_Y_e).c_str());
        CCTK_VINFO("GRHydro_hot_atmo_Y_e has been set to %s", to_string(hot_atmo_Y_e).c_str());
    }
    else CCTK_WARN(2,"No Spritz nor GRHydro?");

  CCTK_INFO ("Setting up initial Temperature, Y_e and entropy.");
  try {
    for (int k=0; k < cctk_lsh[2]; k++) {
      for (int j=0; j < cctk_lsh[1]; j++) {
        for (int i=0; i < cctk_lsh[0]; i++) {
    	  int i3D = CCTK_GFINDEX3D(cctkGH, i, j, k);
          double thislogrho = log10(MAX(rho[i3D], rho_abs_min));

          //now set temperature, entropy and Y_e
          int rhoindex = 0; //first of all, get closest index to the left
          for (int i = 0; i < cactuslogrho.size() && cactuslogrho[i] < thislogrho; ++i) rhoindex = i;

          //then, use it to set electron fraction
          Y_e[i3D] = interp1D(cactuslogrho,rhoindex,thislogrho,ye_beta);
          Y_e[i3D] = MAX(ID_Ye_min,MIN(Y_e[i3D],ID_Ye_max)); //Y_e is constrained into bounds

          if(entropyslice) temperature[i3D] = interp1D(cactuslogrho,rhoindex,thislogrho,temperatures);
          else temperature[i3D] = tempslice;
          temperature[i3D] = MAX(ID_T_min,MIN(temperature[i3D],ID_T_max));

          entropy[i3D] = interp1D(cactuslogrho,rhoindex,thislogrho,betaentropy);
          if(get_EoS_eps) {
              eps[i3D] = interp1D(cactuslogrho,rhoindex,thislogrho,energies);
          }
          //during Prim2ConInitial, temperature and entropy will be overwritten using rho, eps, ye.
          //pressure will be set accordingly using EoS.
          // if(rhoindex > 213){ //DEBUG
          //     cout << thislogrho << " " << temperature[i3D] << " " << Y_e[i3D] << " " << eps[i3D] << endl;
          // }
        }
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::temperature") > 1) {
      memcpy(temperature_p, temperature, local_grid_size*sizeof(CCTK_REAL));
      memcpy(Y_e_p, Y_e, local_grid_size*sizeof(CCTK_REAL));
      if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::temperature") > 2) {
        memcpy(temperature_p_p, temperature, local_grid_size*sizeof(CCTK_REAL));
        memcpy(Y_e_p_p, Y_e, local_grid_size*sizeof(CCTK_REAL));
	  }
    }
    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::eps") > 1) {
      memcpy(eps_p, eps, local_grid_size*sizeof(CCTK_REAL));
      if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::eps") > 2) {
        memcpy(eps_p_p, eps, local_grid_size*sizeof(CCTK_REAL));
	  }
    }
  CCTK_VINFO ("Beta set on on reflevel %d.",int(log10(double(cctk_levfac[0]))/log10(2.0)));
  } catch (ios::failure e) {
    CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Could not read initial beta_equilibrated data from file '%s': %s", beta_datafile, e.what());
  }
//  return;
}

void Spritz_SetBetaEquilibrium_Lorene (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO ("Setting up initial Beta Equilibrated data for Lorene format.");

/*string readme_filename = beta_datafile;
  readme_filename = readme_filename.erase(readme_filename.length()-3) + "README";
  ifstream readme_file(readme_filename.c_str());
  bool ifreadme = readme_file.is_open();
  if(ifreadme){
    readme_file.ignore(1000,'\n');
    string infoline; getline(readme_file,infoline);
    CCTK_VINFO("Betaslice info: %s\n", infoline.c_str());
  }
  readme_file.close();*/

  // read electron fraction data, as function of nb.
    string ye_filename = beta_datafile;
    ye_filename += "ye.txt"; //reading electron fraction at beta equilibrium
    ifstream ye_file(ye_filename.c_str());
    vector<double> ye_beta;
    while(1) { //fill array yearray
        double tmp_ye;
        ye_file >> tmp_ye;
        if(ye_file.eof()) break;
        ye_beta.push_back(tmp_ye);
    }
    ye_file.close();

    /* Federico: This hard-coded value comes from considering 
     * neutron mass as 1.66e-24 g, as in Lorene */
    const double m_neutron_MeV = 931.1916672930752; 
    const double c_siinv2 = 1.1126500560536185E-17; //2.99792458E+8 ;	 ///< Velocity of light [m/s] //same as in Lorene
    const double mev_si = 1.602176565E-13 ;   ///< One MeV [J] //same as in Lorene

    // set up rho, as function of nb and electron fraction previously found
    double geo_to_si_rho = 6.17714E17;
    double nbtorho_convfac = log10(m_neutron_MeV * mev_si * c_siinv2 / geo_to_si_rho) + 42.0;

    /* Federico: read nb, total internal energy and P from the other file */
    string thermo_filename = beta_datafile;
    thermo_filename += "LoreneInput0.01.txt";
    ifstream betaeq_data_file(thermo_filename.c_str());

    /* Federico: Temperature seems to be contained in the first line of the
     * LoreneInput*.txt file, preceeded by # char. 
     * We may want to recover such variable */
    string temperature_line;                                      // Lorenzo
//  betaeq_data_file >> temperature_line;
    /* Federico: strip first character if it is a '#' */
//  if (temperature_line[0] == '#') temperature_line++;
    getline(betaeq_data_file, temperature_line);                  // Lorenzo
    if (temperature_line[0] == '#') temperature_line.erase(0, 1); // Lorenzo

    /* Federico: Convert read temperature string to a double */
    // Lorenzo:  eptr is declared, but not initialised
//  char *eptr;
    double tempslice;
    tempslice = strtod(temperature_line.c_str(), NULL/*&eptr*/);  // Lorenzo
    vector<double> temperatures;
    CCTK_VINFO("Temperature[MeV] will be set to %f", tempslice);
    
    /* Federico: Hard-coded value to ignore the next number of lines */
    int linesign = 8;
    vector<double> cactuslogrho, totenergies, pressures, energies;
    /* Federico: Ignoring first lines */
    for (int i=0;i<linesign;i++) {
        betaeq_data_file.ignore(1000,'\n');
    }

    /* Federico: Reading other lines */
    while(1) {
      double tmpidx, nb, tmpenergy, pr;
      betaeq_data_file >> tmpidx >> nb >> tmpenergy >> pr;
      if(betaeq_data_file.eof()) break;
      cactuslogrho.push_back( log10(nb) + nbtorho_convfac);
      totenergies.push_back(tmpenergy);
      pressures.push_back(pr);
    }
    betaeq_data_file.close();
     
    if(ye_beta.size() != cactuslogrho.size()) CCTK_WARN(0,"Error in EoS tables: ye size not eq to rho size.");
    if(totenergies.size() != cactuslogrho.size()) CCTK_WARN(0,"Error in EoS tables: energy size not eq to rho size.");
    if(pressures.size() != cactuslogrho.size()) CCTK_WARN(0,"Error in EoS tables: pressure size not eq to rho size.");

    /* Federico: Not sure of the following, because the energy 
     * that we read is not the eps variable 
     * It seems that the second column is
     * 
     * totenergy = (1+eps)*rho
     * 
     * therefore
     * 
     * eps = (totenergy/rho) - 1*/
    for (int i=0;i<totenergies.size();i++) {
       energies.push_back( ( totenergies[i] / (pow(10.0,cactuslogrho[i])) ) - 1.0 );
    }

    int const local_grid_size = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
    
    // set up Ye_hot_atmo value, choosing atmo density
    int rhoindex = 0; //first of all, get closest index to the left
    for (int i = 0; i < cactuslogrho.size() && cactuslogrho[i] < log10(rho_abs_min); ++i) rhoindex = i;
    double thislogrho = log10(rho_abs_min);
    double hot_atmo_Y_e = interp1D(cactuslogrho,rhoindex,thislogrho,ye_beta);
    if(CCTK_IsThornActive("Spritz")) {
        int chk = 1;
        chk = CCTK_ParameterSet ("Spritz_hot_atmo_Y_e","Spritz",to_string(hot_atmo_Y_e).c_str());
        if(chk == 0) CCTK_VINFO("Spritz_hot_atmo_Y_e has been set to %s", to_string(hot_atmo_Y_e).c_str());
        else         CCTK_INFO("Error in setting Spritz_hot_atmo_Y_e");
        std::ostringstream streamObj;
        streamObj << pow(10.0,cactuslogrho[0]) * 1.01; //just to be sure that it is greater than the minimum in table
        chk = CCTK_ParameterSet ("rho_abs_eos_min",    "Spritz",streamObj.str().c_str());
        if(chk == 0) CCTK_VINFO("rho_abs_eos_min has been set to %s", streamObj.str().c_str());
        else         CCTK_INFO("Error in setting rho_abs_eos_min");
        if(rho_abs_min < rho_abs_eos_min) {
            CCTK_VINFO("rho_abs_min is %f, whereas the minimum in table is %f.", rho_abs_min, rho_abs_eos_min);
            CCTK_PARAMWARN("rho_abs_min must be greater than the minimum in table.");
        }
    }
    else if(CCTK_IsThornActive("GRHydro")) {
        CCTK_ParameterSet ("GRHydro_hot_atmo_Y_e","GRHydro",to_string(hot_atmo_Y_e).c_str());
        CCTK_VINFO("GRHydro_hot_atmo_Y_e has been set to %s", to_string(hot_atmo_Y_e).c_str());
    }
    else CCTK_WARN(2,"No Spritz nor GRHydro?");

    CCTK_INFO ("Setting up initial Temperature, Y_e and entropy.");
    try {
      for (int k=0; k < cctk_lsh[2]; k++) {
        for (int j=0; j < cctk_lsh[1]; j++) {
          for (int i=0; i < cctk_lsh[0]; i++) {
    	    int i3D = CCTK_GFINDEX3D(cctkGH, i, j, k);
            double thislogrho = log10(MAX(rho[i3D], rho_abs_min));

            //now set temperature, entropy and Y_e
            int rhoindex = 0; //first of all, get closest index to the left
            for (int i = 0; i < cactuslogrho.size() && cactuslogrho[i] < thislogrho; ++i) rhoindex = i;

            //then, use it to set electron fraction
            Y_e[i3D] = interp1D(cactuslogrho,rhoindex,thislogrho,ye_beta);
            Y_e[i3D] = MAX(ID_Ye_min,MIN(Y_e[i3D],ID_Ye_max)); //Y_e is constrained into bounds

            temperature[i3D] = tempslice;
            temperature[i3D] = MAX(ID_T_min,MIN(temperature[i3D],ID_T_max));

            /* Federico: the entropy info seems not to be available in this case... */
            //entropy[i3D] = interp1D(cactuslogrho,rhoindex,thislogrho,betaentropy);

            eps[i3D] = interp1D(cactuslogrho,rhoindex,thislogrho,energies);
            //during Prim2ConInitial, temperature and entropy will be overwritten using rho, eps, ye.
            //pressure will be set accordingly using EoS.
            // if(rhoindex > 213){ //DEBUG
            //     cout << thislogrho << " " << temperature[i3D] << " " << Y_e[i3D] << " " << eps[i3D] << endl;
            // }
          }
        }
      }

      if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::temperature") > 1) {
        memcpy(temperature_p, temperature, local_grid_size*sizeof(CCTK_REAL));
        memcpy(Y_e_p, Y_e, local_grid_size*sizeof(CCTK_REAL));
        if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::temperature") > 2) {
          memcpy(temperature_p_p, temperature, local_grid_size*sizeof(CCTK_REAL));
          memcpy(Y_e_p_p, Y_e, local_grid_size*sizeof(CCTK_REAL));
  	  }
      }
      if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::eps") > 1) {
        memcpy(eps_p, eps, local_grid_size*sizeof(CCTK_REAL));
        if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::eps") > 2) {
          memcpy(eps_p_p, eps, local_grid_size*sizeof(CCTK_REAL));
	    }
      }
      CCTK_VINFO ("Beta set on on reflevel %d.",int(log10(double(cctk_levfac[0]))/log10(2.0)));
      } catch (ios::failure e) {
        CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Could not read initial beta_equilibrated data from file '%s': %s", beta_datafile, e.what());
    }
//  return;
}
