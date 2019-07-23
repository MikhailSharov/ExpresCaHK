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
expresfiles = ['101501_190503.1087.fits','103095_190503.1088.fits','141004_180602.1050.fits','86728_190210.1119.fits','95128_190210.1126.fits','109358_190210.1158.fits','144579_190427.1085.fits','146233_190427.1090.fits','157214_180625.1150.fits','165341_190506.1112.fits']

nsofile = 'NSOspec.csv'                                         #NSO SOLAR .CSV

#RV Correction------------------------------------------------------------------------------------------------------------
"""
The following section will calculate the radial veclocity shift of the object compared to the NSO solar spectrum
"""
nsoflux = np.asarray((pd.read_csv(nsofile)['Flux']))
nsowave = np.asarray((pd.read_csv(nsofile)['Wavelength']))
#DEFINING A REQUIRED FUNCTION
def near_index(arr, num):
    arr = np.asarray(arr)
    index = (np.abs(arr - num)).argmin()
    return index

#MAIN
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
    
#S_HK Calculation---------------------------------------------------------------------------------------------------------
"""
The following section  will calculate the integrated values of the line cores and continuum
"""
avg1 = []
avg2 = []
srpex = []
for i,j in zip(expresfiles, adjust):
    numsrpex = []
    waveimage_file = fits.open(i)
    
#H Line Core Calculation
    xl = []                                         #array reset
    yl = []
    tempx = []
    tempy = []

    image_data = waveimage_file[1].data.copy()    #order import

    for i in range(0, 6500):
        yl.append(image_data['spectrum'][6][i]/image_data['continuum'][6][i]) #general region import
        xl.append(image_data['bary_wavelength'][6][i])

    for l in range(len(xl)):
            xl[l] = xl[l]/j
    for i in range(len(xl)):
        if xl[i] < 3970.69 and xl[i] > 3968.51:   #specific region definition
            if xl[i] < 3969.6:
                if yl[i] > (1.28440)*xl[i] - (5097.16880):
                    yl[i] = (1.28440)*xl[i] - (5097.16880)
                else:
                    yl[i] = yl[i]
            if xl[i] >= 3969.6:
                if yl[i] > (-1.28440)*xl[i] + (5099.96880):
                    yl[i] = (-1.28440)*xl[i] + (5099.96880)
                else:
                    yl[i] = yl[i] 
                    
            tempy.append(yl[i])
            tempx.append(xl[i])

    numsrpex.append(np.trapz(tempy, tempx))     #trapezoidal integration

#K Line Core Calculation
    xl = []                                         #array reset
    yl = []
    tempx = []
    tempy = []

    for i in range(1000, 6500):
        yl.append(image_data['spectrum'][4][i]/image_data['continuum'][4][i]) #general region import
        xl.append(image_data['bary_wavelength'][4][i])

    for l in range(len(xl)):
            xl[l] = xl[l]/j
    for i in range(len(xl)):
        if xl[i] < 3935.79 and xl[i] > 3933.61:
            if xl[i] < 3934.7:
                if yl[i] > (1.28440)*xl[i] - (5052.34311):
                    yl[i] = (1.28440)*xl[i] - (5052.34311)
                else:
                    yl[i] = yl[i]
            if xl[i] >= 3934.7:
                if yl[i] > (-1.28440)*xl[i] + (5055.14311):
                    yl[i] = (-1.28440)*xl[i] + (5055.14311)
                else:
                    yl[i] = yl[i]
            tempy.append(yl[i])
            tempx.append(xl[i])

    numsrpex.append(np.trapz(tempy, tempx))     #trapezoidal integration

#R Continuum Calculation
    xl = []                                         #array reset
    yl = []
    tempx = []
    tempy = []

    for i in range(0, 6500):
        yl.append(image_data['spectrum'][7][i]/image_data['continuum'][7][i]) #general region import
        xl.append(image_data['bary_wavelength'][7][i])

    for l in range(len(xl)):
            xl[l] = xl[l]/j
    for i in range(len(xl)):
        if xl[i] < 4011 and xl[i] > 3991:   #specific region definition
            tempy.append(yl[i])
            tempx.append(xl[i])

    numsrpex.append(np.trapz(tempy, tempx))     #trapezoidal integration

#V Continuum Calculation
    xl = []                                         #array reset
    yl = []
    tempx = []
    tempy = []

    for i in range(0, 6500):
        yl.append(image_data['spectrum'][3][i]/image_data['continuum'][3][i]) #general region import
        xl.append(image_data['bary_wavelength'][3][i])

    for l in range(len(xl)):
            xl[l] = xl[l]/j
    for i in range(len(xl)):
        if xl[i] < 3911 and xl[i] > 3891:   #specific region definition
            tempy.append(yl[i])
            tempx.append(xl[i])

    numsrpex.append(np.trapz(tempy, tempx))     #trapezoidal integration

    avg1.append(numsrpex[0]/numsrpex[1])
    avg2.append(numsrpex[2]/numsrpex[3])
    srpex.append((numsrpex[0]+(numsrpex[1]*0.8021054298156616))/(numsrpex[2]+(numsrpex[3]*1.0457836985441216)))

#print(np.mean(avg1),np.mean(avg2))  #Uncomment this line to see an updated alpha and beta variables that would replace the values above  
    
    
#(OPTIONAL) (Only if NEW FIT required)----------------------------------------------------------------------------------------------------------------
"""
Should a new linear fit be required, for example if many new objects are added and the linear fit would improve from an update, the following block 
can be uncommented.
This will require knowledge of MWO S-values, which would be placed in the array ys, in the same order as expresfiles previously.
This path, when uncommented will print out the linear variables that are the best fit (least squares) of data
"""

"""
ys = np.array([])
x = np.asarray(srpex)
y = np.asarray(ys)

#GENERAL FITTING EQUATION       #uncomment whichever equation works best, deafult = linear
def func(x, A, c, d):
#    return A*np.exp(c*x) + d
#    return (A*((x-c)**2))+d
     return (A*x) + c

x_lin = np.linspace(0, xs.max(), 100)
p0 = [20, -0.120, 0.12]                                        # guessed params
w, _ = opt.curve_fit(func, xs, ys, p0=p0, bounds=([-np.inf,-np.inf,-np.inf], [np.inf,np.inf,0.1]), maxfev = 10000000)     
print("Estimated Parameters", w)  

#Model
y_model = func(x_lin, *w)

#PLOT
#Visualize data and fitted curves
plt.clf()
plt.plot(xs, ys, "ko", label="Data")
plt.plot(x_lin, y_model, "k--", label="Fit")
plt.title("Least squares regression")
plt.legend(loc="upper left")
"""

#Final Calibration--------------------------------------------------------------------------------------------------
"""
The following section calibrates the values for S_EXPRES and prints these values
"""

for i in range(len(srpex)):
    srpex[i] = (6.10546*srpex[i]) + 0.024619             #These values are either known or calculated in previous section

print(srpex)


