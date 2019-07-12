# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 11:11:00 2019

@author: msharov
"""
#IMPORTING PACKAGES-------------------------------------------------------------------------------------
import pandas as pd
import numpy as np
from astropy import units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Box1DKernel
from scipy import *
from scipy.optimize import leastsq
from math import *
import scipy.optimize as opt
from PyAstronomy import pyasl

#IMPORTING FILES----------------------------------------------------------------------------------------
expresfiles = ['165341_190506.1112.fits']
nsofile = 'NSOspec.csv'

returns = []

nsoflux = np.asarray((pd.read_csv(nsofile)['Flux']))
nsowave = np.asarray((pd.read_csv(nsofile)['Wavelength']))
#DEFINING A REQUIRED FUNCTION ----------------------------------------------------------------------------
def near_index(arr, num):
    arr = np.asarray(arr)
    index = (np.abs(arr - num)).argmin()
    return index

#MAIN------------------------------------------------------------------------------------------------

def fun(x):
    temp = []
    expressx = expresx.copy()
    for k in range(len(expresx)):
        expressx[k] = expressx[k]/x
    for k in range(len(expressx)):
        if expressx[k] > 3961 and expressx[k] < 3963:
            nsoequi = near_index(nsox, expressx[k])
            temp.append(((expresy[k] - nsoy[nsoequi])**2)/(nsoy[nsoequi]))
    return sum(temp)

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
        ytemp.append(fun(j))

plt.clf()
#plt.plot(expresx,expresy,label = 'expres RV Corrected')
#plt.plot(nsox, nsoy, label = 'NSO')
#plt.legend()
plt.plot(xtemp,ytemp)

print(len(nsox))
print(len(expresx))