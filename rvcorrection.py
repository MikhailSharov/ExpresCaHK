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

#IMPORTING FILES----------------------------------------------------------------------------------------
expresfiles = ['101501_190503.1087.fits','103095_190503.1088.fits','141004_180602.1050.fits','86728_190210.1119.fits','95128_190210.1126.fits','109358_190210.1158.fits','144579_190427.1085.fits','146233_190427.1090.fits','157214_180625.1150.fits','165341_190506.1112.fits']
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
adjust = []
def fun(x):
    temp = []
    expressx = expresx.copy()
    for k in range(len(expresx)):
        expressx[k] = expressx[k]/x
    for k in range(len(expressx)):
        if expressx[k] > 5169 and expressx[k] < 5201:
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
    
    for j in range(len(expresdata['bary_wavelength'][42])):
        expresflux.append(expresdata['spectrum'][42][j]/expresdata['continuum'][42][j])
        expreswave.append(expresdata['bary_wavelength'][42][j])
    for j in range(len(expreswave)):
        if expreswave[j] > 5170 and expreswave[j] < 5175:
            expresx.append(expreswave[j]/1)
            expresy.append(expresflux[j])
    for j in range(len(nsowave)):
        if nsowave[j] > 5169 and nsowave[j] < 5176:
            nsox.append(nsowave[j])
            nsoy.append(nsoflux[j])
    nsox = nsox[::-1]
    nsoy = nsoy[::-1]
    #This part (below) could easily be replaced with scipy.optimize.minimize
    for j in np.arange(0.999,1.0005,0.000001):
        xtemp.append(j)
        ytemp.append(fun(j))
    adjust.append(xtemp[np.asarray(ytemp).argmin()])
print(adjust)
plt.clf()
#plt.plot(expresx,expresy,label = 'expres RV Corrected')
#plt.plot(nsox, nsoy, label = 'NSO')
#plt.legend()
plt.plot(xtemp,ytemp)

print(len(nsox))
print(len(expresx))