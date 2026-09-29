#!/usr/bin/env python
import h5py
import numpy as np
import linecache
import math
import sys
import os

def interpparab1D(x,y): #return intersection with 0
    if(x.size!=3 or y.size!=3):
        print("Size is different from three!!")
        print(x,y)
    #get quadratic coefficients and delta
    a = (x[1]*y[0] - x[2]*y[0] - x[0]*y[1] + x[2]*y[1] + x[0]*y[2] - x[1]*y[2]) / ((x[0] - x[1])*(x[0] - x[2])*(x[1] - x[2]))
    b = (-(x[1]**2*y[0]) + x[2]**2*y[0] + x[0]**2*y[1] - x[2]**2*y[1] - x[0]**2*y[2] + x[1]**2*y[2]) / ((x[0] - x[1])*(x[0] - x[2])*(x[1] - x[2]))
    c = (x[1]**2*x[2]*y[0] - x[1]*x[2]**2*y[0] - x[0]**2*x[2]*y[1] + x[0]*x[2]**2*y[1] + x[0]**2*x[1]*y[2] - x[0]*x[1]**2*y[2])/((-x[0] + x[1])*(x[1] - x[2])*(-x[0] + x[2]))
    delta = b**2 - 4*a*c
    if(delta < 0): return None
    if(a==0):
        if(b==0): return None
        else:     return -c/b
    res1 = (-b+math.sqrt(delta))/(2*a)
    res2 = (-b-math.sqrt(delta))/(2*a)
    check1 = x[0] <= res1 <= x[2] #if 1st solution is in range
    check2 = x[0] <= res2 <= x[2] #if 2nd solution is in range
    if( check1 and check2 and delta>0 ): return None #if two different solutions in range
    elif(check1): return res1
    elif(check2): return res2
    else: return None


def lookforintersection(index, xarray, yarray):
    neighbourhood  = yarray[max(0,index-3):index+3]
    xneighbourhood = xarray[max(0,index-3):index+3]
    if (np.all(neighbourhood < 0) or np.all(neighbourhood > 0)): #if there's no intersection in range, do nothing.
        return None,None
    else: #else, get the switching point
        for i,_ in enumerate(neighbourhood[:-1]):
            if (neighbourhood[i+1]*neighbourhood[i] <= 0):
                #switching point has been found
                swp = i
                pass
        if(swp==0):
            startpoint = 0
        elif(swp>=neighbourhood.size-3):
            startpoint = neighbourhood.size-3
        else:
            startpoint = swp if(abs(neighbourhood[swp+1])<abs(neighbourhood[swp])) else swp-1
    threex = xneighbourhood[startpoint:startpoint+3]
    threey =  neighbourhood[startpoint:startpoint+3]

    xval      = interpparab1D(threex,threey)
    xidx_prev = None if xval is None else (np.argmax(xarray > xval) - 1) #previous index on x array
    return xval, xidx_prev

def beta_computation(nb_n,temp_n,mu_l,ye):
    #mu_p = mu_q + mu_b + mp
    #mu_n = mu_b + mn
    #mu_e = mu_l - mu_q + me
    #mu_l = mu_e + mu_p - mu_n + mn - mp -me = mu_e - muhat - me

    print("Computing Ye in beta-equilibrium ... ",end='')
    ye_beta = np.empty([temp_n,nb_n])
    mu_l_beta = np.empty([temp_n,nb_n])
    ye_beta_index = np.empty([temp_n,nb_n])

    for rhoindex in range(0,nb_n):
        for tempindex in range(0,temp_n):
            yeminindex = np.argmin(np.absolute(mu_l[:,tempindex,rhoindex])) #this is the index of the nearest ye to zero
            ye_beta_index[tempindex,rhoindex] = yeminindex
            # An interpolation gives better precision: the following line looks for the intersection in the neighbourhood of yeminindex
            ye_beta_tmp, ye_beta_index_tmp    = lookforintersection(yeminindex, ye, mu_l[:,tempindex,rhoindex])
            ye_beta[tempindex,rhoindex]       = ye[yeminindex] if ye_beta_tmp is None else ye_beta_tmp
            ye_beta_index[tempindex,rhoindex] = yeminindex if ye_beta_index_tmp is None else ye_beta_index_tmp
            mu_l_beta[tempindex,rhoindex]     = mu_l[yeminindex,tempindex,rhoindex] if ye_beta_tmp is None else 0.
    print("Done.")
    return ye_beta, mu_l_beta, ye_beta_index
