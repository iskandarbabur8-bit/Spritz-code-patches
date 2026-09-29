#!/usr/bin/env python
import h5py
import numpy as np
import linecache
import math
import sys
import os
import re
import datetime
import scipy.constants as const
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from shutil import copy2
from scipy import interpolate

m_n = 1.6749274125E-24 #g
mev_si = 1.602176565E-13
m_n_mev = 939.565413 # c_light = 299792458
m_p_mev = 938.272081 # c_light = 299792458
m_e_mev = 0.510998946

logkappa = math.log10(const.c**2 * 1.0E4) # #useful constant for energy conversion
pressconv = 17.980714384531204 # math.log10(m_n * 1.0E-3 * 1.0E-4 / mev_si)  #useful constant for pressure conversion


def getinterpqty(xarray,thisx,xidx,yarray):
    if(xidx < 4 or xidx == len(xarray)-5):
        return linearinterpo(xarray,thisx,xidx,yarray)
    else:
        try: #standard interpolation with cubic splines
            f = interpolate.interp1d(xarray[xidx-1:xidx+3],yarray[xidx-1:xidx+3],kind='cubic')
            res = f(thisx)
        except:
            # spline is more accurate than lagrange. We use the latter if too small domain.
            poly = interpolate.lagrange(xarray[xidx-1:xidx+3],yarray[xidx-1:xidx+3])
            res = poly(thisx)
        return res

def linearinterpo(xarray,thisx,xidx,yarray):
    if(xidx == len(xarray) - 1): xidx -= 1 #if xidx is at border, decrease it for extrapolation
    x1 = xarray[xidx]
    x2 = xarray[xidx+1]
    dx = thisx - x1
    y1 = yarray[xidx]
    y2 = yarray[xidx+1]
    return y1 + (dx)*(y2-y1)/(x2-x1)

def get2Dinterpqty(xarray,thisx,xidx,yarray,thisy,yidx,zarray):
    #cubic splines on both x and y.
    if(thisx>xarray[-1] or (thisx<xarray[0] and (xarray[0]-thisx)/xarray[0]>0.1)): print("Warning on x",thisx,xarray[0])
    if(thisy>yarray[-1] or thisy<yarray[0]): print("Warning on y",thisy,yarray[0])
    thisx = max(thisx,xarray[0])
    xstencil = range(max(0,xidx-2),min(xidx+3,xarray.size-1))
    ystencil = range(max(0,yidx-2),min(yidx+3,yarray.size-1))

    tempoy = np.empty([len(xstencil)])
    for i,ix in enumerate(xstencil):
        try:    tempoy[i] = interpolate.interp1d(yarray[ystencil],zarray[ix,ystencil],kind='cubic')(thisy)
        except: tempoy[i] = interpolate.lagrange(yarray[ystencil],zarray[ix,ystencil])(thisy) #lagrange if not enough points to compute derivatives
    #tempy contains y-interpolated values in thisy
    try:    res = interpolate.interp1d(xarray[xstencil],tempoy,kind='cubic')(thisx)
    except: res = interpolate.lagrange(xarray[xstencil],tempoy)(thisx) #lagrange if not enough points to compute derivatives
    return res

def writeonreadme(path, detailstring):
    print("Writing README on",path,"... ")
    outfile = open(path,'w')
    outfile.write("Slice taken at " + detailstring + "\n")
    outfile.write("Slice has been taken on " + str(datetime.datetime.now()) + ", using " + sys.argv[0] + " on " + os.path.basename(sys.argv[1]) + "\n")
    outfile.write("Output includes:" + "\n")
    outfile.write("\teos.nb \t- 2-line header, then baryon density gridpoints [fm-3]\n")
    outfile.write("\teos.t \t- 2-line header, then temperature gridpoints [MeV], only one point\n")
    if 'entropy' in detailstring: outfile.write("\teos.t_s \t- Temperature [MeV] as function of density at constant entropy\n")
    outfile.write("\teos.yq \t- 2-line header, then electron fraction gridpoints\n")
    outfile.write("\teos.yqb \t- electron fraction at beta-equilibrium\n")
    outfile.write("\teos.thermo \t- 1-line header, then (CompOSE units):\n")
    outfile.write("\t\t tempindex | rhoindex | yeindex(beta-equilibrated) | pressure | entropy | 4 dummy numbers | eps(c^2 units)\n")
    outfile.close()
    return

def printplots(logrho,yebetaarray,pressarray,temparray,mode,sliceval,savepath):
    print("Saving figures on",savepath)

    sizeratios = [2,1,1] if mode=='s' else [3,1]
    slicestr = 'S = ' + sliceval + ' kb' if (mode=='s') else 'T = ' + sliceval + ' MeV'

    f, canvas = plt.subplots(len(sizeratios),1,sharex=True,gridspec_kw = {'height_ratios':sizeratios})
    f.suptitle('Beta-equilibrated slice at ' + slicestr)

    canvas[0].grid()
    canvas[0].plot(10**logrho,10**pressarray)
    canvas[0].set_yscale('log')
    canvas[0].set_xscale('log')
    canvas[0].set_ylabel(r'Pressure [dyne/cm${}^2$]')

    canvas[1].plot(10**logrho,yebetaarray,'C3')
    canvas[1].set_xlabel(r'Rest Mass Density [g/cm${}^3$]')
    canvas[1].set_ylabel(r'$Y_{e,\beta eq}$')

    if(mode=='s'):
        canvas[2].grid()
        canvas[2].plot(10**logrho,10**temparray,'C5')
        canvas[2].set_yscale('log')
        canvas[2].set_xlabel(r'Rest Mass Density [g/cm${}^3$]')
        canvas[2].set_ylabel('T [MeV]')

    f.subplots_adjust(hspace=0)
    plt.savefig(savepath,dpi=500)
    return

def BE_fullplot(ye_beta, logtemp, logrho,savepath):
    plt.figure()
    X, Y = np.meshgrid(logrho, logtemp)
    Z = ye_beta
    contours = plt.contour(X, Y, Z, [0.05,0.10,0.20,0.30,0.40,0.50], colors='black',linewidths=0.8)
    plt.clabel(contours, inline=True, fontsize=8)
    plt.contourf(X, Y, Z, 50, cmap='viridis',alpha=0.75)
    cbar = plt.colorbar();
    plt.xlabel(r'Rest Mass Density [g/cm${}^3$, logscale]')
    plt.ylabel(r'Logarithmic Temperature [MeV, logscale]')
    cbar.set_label('Electron Fraction', rotation=90)
    print("Saving figures on",savepath)
    plt.savefig(savepath,dpi=1000)
    return

def entropyplot(entropybeta,logtemp,logrho,chosenentropy,savepath): #prints plot of entropy in beta-equilibrium, emphasizing the chosen one
    plt.figure()
    X, Y = np.meshgrid(logrho,logtemp)
    Z = entropybeta.transpose()
    contours = plt.contour(X, Y, Z, chosenentropy, colors='white',linewidths=1.2)
    plt.clabel(contours, inline=True, fontsize=10, color='white')
    levels = np.logspace(math.log10(max(1.0E-3,Z.min())), math.log10(Z.max()), 50)
    plt.contourf(X, Y, Z, levels=levels, cmap='viridis',alpha=0.85,norm=colors.LogNorm())
    plt.xlabel(r'Rest Mass Density [g/cm${}^3$, logscale]')
    plt.ylabel(r'Logarithmic Temperature [MeV, logscale]')
    cbar = plt.colorbar();
    cbar.set_label(r'Entropy [$k_b$/baryon]', rotation=90)
    cbar.set_ticks(np.logspace(-3, 12,6))
    print("Saving figures on",savepath)
    plt.savefig(savepath,dpi=1000)
    return

def writestandardtxt(slicedfilepath, nb_n, ye_n, logrho, ye_beta_slice, ye):

    with open(slicedfilepath + "/eos.nb", 'w') as outfile:
        convfac = math.log10(m_n) + 39.0 #fm-3 to cm-3
        print("Writing results on",outfile.name,"...")
        outfile.write(str(1)+ "\n")
        outfile.write(str(nb_n) + "\n")
        for rho_index in range(0,nb_n):
            nb = math.pow(10,logrho[rho_index] - convfac)
            outfile.write(str(nb) + "\n")

    with open(slicedfilepath + "/eos.yqb", 'w') as outfile:
        print("Writing results on",outfile.name,"...")
        for rho_index in range(0,nb_n):
            outfile.write(str(ye_beta_slice[rho_index]) + "\n")

    with open(slicedfilepath + "/eos.yq", 'w') as outfile:
        print("Writing results on",outfile.name,"...")
        outfile.write(str(1)+ "\n")
        outfile.write(str(ye_n) + "\n")
        for ye_index in range(0,ye_n):
            outfile.write(str(ye[ye_index]) + "\n")
