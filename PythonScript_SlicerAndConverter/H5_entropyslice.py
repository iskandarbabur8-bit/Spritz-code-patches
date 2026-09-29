#!/usr/bin/env python
import h5py
import numpy as np
import linecache
import math
import sys
import os
from scipy import interpolate
from S2C_utilities import *
from beta_calc import *


def entropyslicer():
    if (len(sys.argv) < 2):
        sys.exit("ERROR: Invalid number of arguments: you need to pass as argument the EOS folder, which should be the EOS name.")
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

    print("Slicer: choose an entropy value and you'll get a slice. Choose in range [kb units], entropy per baryon ",np.amin(entropy),"-->",np.amax(entropy),end='... ')
    chosenentropy = float(input())
    slicedfilepath = "../" + os.path.basename(tablepath)[:-3] + "_entropyslice" + str(int(chosenentropy))

    if not (os.path.isdir(slicedfilepath)):
        print("Creating directory ",slicedfilepath)
        os.mkdir(slicedfilepath)

    mu_l = mu_e - muhat - m_e_mev
    ye_beta, mu_l_beta, ye_beta_index = beta_computation(nb_n,temp_n,mu_l,ye)
    BE_fullplot(ye_beta,logtemp,logrho,slicedfilepath+'/BE_plots.png')

    entropybeta = np.zeros([nb_n,temp_n])
    for rho_index, _ in enumerate(logrho):
        for temp_index, _ in enumerate(logtemp):
            thisye    = ye_beta[temp_index,rho_index]
            thisyeidx = int(ye_beta_index[temp_index,rho_index])
            entropybeta[rho_index,temp_index] = getinterpqty(ye,thisye,thisyeidx,entropy[:, temp_index, rho_index]) #interpolate in ye to get beta-eq entropy.
    entropyplot(entropybeta,logtemp,logrho,chosenentropy,slicedfilepath+'/Entropy_plots.png')

    # #now get contourplot at selected entropy and get values from it.
    get_contour_with_pyplot = True #set True if you want to use matplotlib, False if you prefer the other routine.
    if(get_contour_with_pyplot):
        print("Slicing entropy with matplotlib.")
        X, Y = np.meshgrid(logtemp,logrho)
        cs = plt.contour(Y,X,entropybeta, [chosenentropy], cmap='viridis');
        contourlevel = np.asarray(cs.collections[0].get_paths()[0].vertices.T); #get contour values
        entropy_interp = interpolate.interp1d(contourlevel[0],contourlevel[1],kind='cubic',fill_value=logtemp[0],bounds_error=False); #intepolate contour values
        temp_entropyslice = np.full([len(logrho)],logtemp[0])
        for rhoindex,logrhoval in enumerate(logrho):
            temp_entropyslice[rhoindex] = max(min(entropy_interp(logrhoval),logtemp[-1]),logtemp[0])
    else:
        print("Slicing entropy with parabolic interpolation.")
        temp_entropyslice = np.full([len(logrho)],logtemp[0])
        for rhoindex,_ in enumerate(logrho):
            Sminindex = np.argmin(np.absolute(entropybeta[rhoindex,:] - chosenentropy))
            temp_tmp, _ = lookforintersection(Sminindex, logtemp, entropybeta[rhoindex,:] - chosenentropy)
            temp_entropyslice[rhoindex] = logtemp[0] if temp_tmp is None else temp_tmp
    #now we have temp_entropyslice[rhoindex]

    ye_beta_slice       = np.empty([len(logrho)])
    pressure_beta_slice = np.empty([len(logrho)])
    energy_beta_slice   = np.empty([len(logrho)])
    nb_array            = np.empty([len(logrho)])
    dummy = str(0)+" "

    with open(slicedfilepath + "/eos" + ".thermo", 'w') as outfile:
        print("Writing results on",outfile.name,"...")
        thermoheader = "{:.8f}   {:.8f}   1\n".format(m_n_mev,m_p_mev)
        outfile.write(thermoheader)  #.thermo header
        for rho_index in range(0,nb_n):
            logtempval = temp_entropyslice[rho_index]
            temp_index = np.argmax(logtemp > logtempval) - 1
            yeval = getinterpqty(logtemp,logtempval,temp_index,ye_beta[:,rho_index])
            ye_index = int(ye_beta_index[temp_index,rho_index])

            thislogpress = get2Dinterpqty(ye,yeval,ye_index,
                                          logtemp,logtempval,temp_index,
                                          logpress[:,:,rho_index])
            qty1         = 10**(thislogpress - logrho[rho_index] - pressconv)

            thislogenergy = get2Dinterpqty(ye,yeval,ye_index,
                                          logtemp,logtempval,temp_index,
                                          logenergy[:,:,rho_index])
            qty7          = 10**(thislogenergy - logkappa) - 10**(math.log10(energy_shift) - logkappa)

            qty2 = chosenentropy

            outfile.write(str(1) + " " + str(rho_index) + " " + str(ye_index) + " " + str(qty1) + " " + str(qty2) + " " + dummy + dummy + dummy + dummy + str(qty7) + "\n")
            pressure_beta_slice[rho_index] = thislogpress;
            energy_beta_slice[rho_index] = thislogenergy;
            ye_beta_slice[rho_index] = yeval


    with open(slicedfilepath + "/eos.t_s", 'w') as outfile:
        print("Writing results on",outfile.name,"...")
        for rho_index in range(0,nb_n):
            tmptemp = pow(10,temp_entropyslice[rho_index])
            outfile.write(str(tmptemp) + "\n")

    with open(slicedfilepath + "/eos.t", 'w') as outfile:
        print("Writing results on",outfile.name,"...")
        outfile.write(str(1) + "\n")
        outfile.write(str(1) + "\n")
        outfile.write(str(0))

    writestandardtxt(slicedfilepath, nb_n, ye_n, logrho, ye_beta_slice, ye)

    writeonreadme(slicedfilepath + "/README", 'specific entropy of '+str(chosenentropy)+' kb/baryon.')
    printplots(logrho,ye_beta_slice,pressure_beta_slice,temp_entropyslice,'s',str(int(chosenentropy)),slicedfilepath+'/plots.png')

    return slicedfilepath


if __name__ == "__main__":
    entropyslicer()
