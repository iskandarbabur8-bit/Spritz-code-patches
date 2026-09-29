#!/usr/bin/env python
import h5py
import numpy as np
import linecache
import math
import sys
import os
from beta_calc import *
from S2C_utilities import *


def betatable():
    if (len(sys.argv) < 2):
        sys.exit("ERROR: Invalid number of arguments: you need to pass as argument the path to h5file.")
    tablepath = sys.argv[1]

    print("Reading from .h5 file ",tablepath)
    hdf5file        = h5py.File(tablepath, "r")

    nb_n            = np.array(hdf5file['pointsrho'])[0]
    temp_n          = np.array(hdf5file['pointstemp'])[0]
    ye_n            = np.array(hdf5file['pointsye'])[0]
    logrho          = np.array(hdf5file['logrho'])
    logtemp         = np.array(hdf5file['logtemp'])
    ye              = np.array(hdf5file['ye'])
    mu_e            = np.array(hdf5file['mu_e'])
    # mu_p            = np.array(hdf5file['mu_p'])
    # mu_n            = np.array(hdf5file['mu_n'])
    muhat           = np.array(hdf5file['muhat'])
    munu            = np.array(hdf5file['munu'])

    hdf5file.close()
    print("Data have been read successfully.")

    slicedfilepath = "../" + os.path.basename(tablepath)[:-3] + "_betah5file"
    if not (os.path.isdir(slicedfilepath)):
        print("Creating directory ",slicedfilepath)
        os.mkdir(slicedfilepath)

    # mu_l = munu + muhat
    mu_l = mu_e - muhat - 0.510998946
    ye_beta, mu_l_beta, ye_beta_index = beta_computation(nb_n,temp_n,mu_l,ye)
    BE_fullplot(ye_beta,logtemp,logrho,slicedfilepath+'/BE_plots.png')

    print("Writing on .h5 file...")
    hdf5path = slicedfilepath+'/betaslice.h5'
    hdf5file = h5py.File(hdf5path, "w")
    hdf5file.create_dataset("logrho",       data = logrho,       dtype='<f8')
    hdf5file.create_dataset("logtemp",      data = logtemp,      dtype='<f8')
    hdf5file.create_dataset("ye_beta",      data = ye_beta,      dtype='<f8')
    hdf5file.create_dataset("mu_l_beta",    data = mu_l_beta,    dtype='<f8')
    hdf5file.close()
    print("HDF5 file saved to", hdf5path)

if __name__ == "__main__":
    betatable()
