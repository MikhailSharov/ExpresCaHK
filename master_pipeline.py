#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: mikhailsharov

Purpose: The purpose of this pipeline is that given some echelle spectrum from
EXPRES, extracted via Ryan Petersburg extraction method, the output of the pipeline
will be a S-value, analogous to that of the MWO.

Inputs: An array of .fits spectrum files, which are present in the working directory of the
pipeline.

Outputs: An array of S-values, which are linearly calibrated to the MWO values.

"""

import numpy as np
from astropy import units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Box1DKernel
from scipy import *
import scipy.optimize as opt
import pandas as pd

#INSERT FILE NAMES HERE AS AN ARRAY:
expresfiles = []

nsofile = 'NSOspec.csv'                                         #NSO SOLAR .CSV

#RV Correction-----------------------------------------------------------------

nsoflux = np.asarray((pd.read_csv(nsofile)['Flux']))
nsowave = np.asarray((pd.read_csv(nsofile)['Wavelength']))

#The following function finds the index of the most similar element in another array
def near_index(arr, num):
    arr = np.asarray(arr)
    index = (np.abs(arr - num)).argmin()
    return index
#The following function finds the chi squared of a particular region
def chis(x):
    chi = []
    expressx = expresx.copy()
    for k in range(len(expresx)):
        expressx[k] = expressx[k]/x
    for k in range(len(expressx)):
        if expressx[k] > 3961 and expressx[k] < 3963:
            nsoequi = near_index(nsox, expressx[k])
            chi.append(((expresy[k] - nsoy[nsoequi])**2)/(nsoy[nsoequi]))
    return sum(chi)

for i in range(len(expresfiles)):
    expresflux = []
    expreswave = []
    expresx = []
    expresy = []
    xtemp = []
    ytemp = []
    nsox = []
    nsoy = []
    expresfile = fits.open(expresfiles[i])
    expresdata = expresfile[1].data.copy()
    
    for j in range(len(expresdata['bary_wavelength'][6])):
        expresflux.append(expresdata['spectrum'][6][j]/expresdata['continuum'][6][j])
        expreswave.append(expresdata['bary_wavelength'][6][j])
    for j in range(len(expreswave)):
        if expreswave[j] > 3960 and expreswave[j] < 3965:
            expresx.append(expreswave[j]/1)
            expresy.append(expresflux[j])
    for j in range(len(nsowave)):
        if nsowave[j] > 3959 and nsowave[j] < 3966:
            nsox.append(nsowave[j])
            nsoy.append(nsoflux[j])
    nsox = nsox[::-1]
    nsoy = nsoy[::-1]
    
    for j in np.arange(0.999,1.001,0.000001):
        xtemp.append(j)
        ytemp.append(chis(j))

