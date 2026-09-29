#!/usr/bin/env python
import h5py
import numpy as np
import linecache
import math
import sys
import os
from S2C_utilities import *
from beta_calc import *

def temperatureslicer():
    if (len(sys.argv) < 2):
        sys.exit("ERROR: Invalid number of arguments: you need to pass as argument the path to h5file.")
    tablepath = sys.argv[1]

    print("Reading from .h5 file ",tablepath)
    with h5py.File(tablepath, "r") as hdf5file:
        nb_n            = np.array(hdf5file['pointsrho'])[0]
        temp_n          = np.array(hdf5file['pointstemp'])[0]
        ye_n            = np.array(hdf5file['pointsye'])[0]
        logrho          = np.array(hdf5file['logrho'])
        logtemp         = np.array(hdf5file['logtemp'])
        ye              = np.array(hdf5file['ye'])
        energy_shift    = np.array(hdf5file['energy_shift'])[0]
        logenergy       = np.array(hdf5file['logenergy'])
        logpress        = np.array(hdf5file['logpress'])
        entropy        = np.array(hdf5file['entropy'])
        mu_e            = np.array(hdf5file['mu_e'])
        muhat           = np.array(hdf5file['muhat'])
        munu            = np.array(hdf5file['munu'])
    print("Data have been read successfully.")

    print("Slicer: choose a temperature value and you'll get a slice. Choose in range [MeV]",pow(10,logtemp[0]),"-->",pow(10,logtemp[-1]),end='... ')
    chosentemp = max(pow(10,logtemp[0]),min(float(input()),pow(10,logtemp[-1])))
    chosenlogtemp = math.log10(chosentemp)
    tempindex = np.argmax(logtemp > chosenlogtemp) - 1
    print("Slice will be taken at temperature index",tempindex)
    slicedfilepath = "../" + os.path.basename(tablepath)[:-3] + "_slice" + str(tempindex)

    if not (os.path.isdir(slicedfilepath)):
        print("Creating directory ",slicedfilepath)
        os.mkdir(slicedfilepath)

    # mu_l = munu + muhat
    mu_l = mu_e - muhat - m_e_mev
    ye_beta, mu_l_beta, ye_beta_index = beta_computation(nb_n,temp_n,mu_l,ye)
    BE_fullplot(ye_beta,logtemp,logrho,slicedfilepath+'/BE_plots.png')

    ye_beta_index_slice = ye_beta_index[int(tempindex),:] #start slicing ye_beta_index
    ye_beta_slice       = ye_beta[int(tempindex),:]
    pressure_beta_slice = np.empty([len(logrho)])
    energy_beta_slice   = np.empty([len(logrho)])
    dummy = str(0)+" "

    with open(slicedfilepath + "/eos" + ".thermo", 'w') as outfile:
        print("Slicing at Temp index", tempindex)
        print("Writing results on",outfile.name,"...")
        thermoheader = "{:.8f}   {:.8f}   1\n".format(m_n_mev,m_p_mev)
        outfile.write(thermoheader)  #.thermo header
        for rho_index in range(0,nb_n): #for every rho, at fixed temperature, interpolate over ye to get quantities
            ye_index = int(ye_beta_index_slice[rho_index])
            thisye = ye_beta_slice[rho_index]
            thislogpress  = getinterpqty(ye, thisye, ye_index, logpress[:,tempindex,rho_index])
            thisentropy   = getinterpqty(ye, thisye, ye_index, entropy[:,tempindex,rho_index])
            thislogenergy = getinterpqty(ye, thisye, ye_index, logenergy[:,tempindex,rho_index])

            qty1 = math.pow(10, thislogpress - logrho[rho_index] - pressconv)
            qty2 = thisentropy
            qty7 = 10**(thislogenergy - logkappa) - 10**(math.log10(energy_shift) - logkappa)
            outfile.write(str(1) + " " + str(rho_index) + " " + str(ye_index) + " " + str(qty1) + " " + str(qty2) + " " + dummy + dummy + dummy + dummy + str(qty7) + "\n")
            pressure_beta_slice[rho_index] = thislogpress;
            energy_beta_slice[rho_index]   = thislogenergy;

    with open(slicedfilepath + "/eos.t",'w') as outfile:
        print("Writing results on",outfile.name,"...")
        outfile.write(str(1) + "\n")
        outfile.write(str(1) + "\n")
        outfile.write(str(math.pow(10,logtemp[tempindex])))

    writestandardtxt(slicedfilepath, nb_n, ye_n, logrho, ye_beta_slice, ye)

    writeonreadme(slicedfilepath + "/README", 'temperature of '+str(math.pow(10,logtemp[tempindex]))+' MeV.')
    printplots(logrho,ye_beta_slice,pressure_beta_slice,energy_beta_slice,'t',str(math.pow(10,logtemp[tempindex])),slicedfilepath+'/plots.png')
    return slicedfilepath

if __name__ == "__main__":
    temperatureslicer()
