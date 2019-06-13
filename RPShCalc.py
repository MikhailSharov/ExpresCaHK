#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 10:05:05 2019

@author: msharov
"""
import numpy as np
from astropy import units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Box1DKernel
from scipy import *

nums = []
image_file = fits.open('114783_190503.1098.fits')

image_file.info()
image_data = image_file[1].data.copy()    #order import
image_file.close()
signoi = []

#H Core Line

xl = []                                         #array reset
yl = []
noise = []
signal = []
tempx = []
tempy = []

for i in range(0, 6500):
    yl.append(image_data['spectrum'][6][i]/image_data['continuum'][6][i]) #general region import
    xl.append(image_data['wavelength'][6][i])
    noise.append(image_data['uncertainty'][6][i])
    signal.append(image_data['spectrum'][6][i])
    

for i in range(len(xl)):
    if xl[i] < 3970.445 and xl[i] > 3969.355:   #specific region definition
        tempy.append(yl[i])
        tempx.append(xl[i])
        signoi.append(signal[i]/noise[i])

nums.append(np.trapz(tempy, tempx))     #trapezoidal integration
print(np.average(signoi))
#K Core Line

xl = []                                         #array reset
yl = []
tempx = []
tempy = []

for i in range(4000, 6500):
    yl.append(image_data['spectrum'][4][i]/image_data['continuum'][4][i]) #general region import
    xl.append(image_data['wavelength'][4][i])

for i in range(len(xl)):
    if xl[i] < 3935.545 and xl[i] > 3934.455:   #specific region definition
        tempy.append(yl[i])
        tempx.append(xl[i])

nums.append(np.trapz(tempy, tempx))     #trapezoidal integration

#R Cont Line

xl = []                                         #array reset
yl = []
tempx = []
tempy = []

for i in range(0, 6500):
    yl.append(image_data['spectrum'][7][i]/image_data['continuum'][7][i]) #general region import
    xl.append(image_data['wavelength'][7][i])

for i in range(len(xl)):
    if xl[i] < 4011.9 and xl[i] > 3991.9:   #specific region definition
        tempy.append(yl[i])
        tempx.append(xl[i])

nums.append(np.trapz(tempy, tempx))     #trapezoidal integration

#V Cont Line

xl = []                                         #array reset
yl = []
tempx = []
tempy = []

for i in range(0, 6500):
    yl.append(image_data['spectrum'][3][i]/image_data['continuum'][3][i]) #general region import
    xl.append(image_data['wavelength'][3][i])

for i in range(len(xl)):
    if xl[i] < 3911 and xl[i] > 3891:   #specific region definition
        tempy.append(yl[i])
        tempx.append(xl[i])

nums.append(np.trapz(tempy, tempx))     #trapezoidal integration


print(nums[0],nums[1],nums[2],nums[3])
print(nums[0]/nums[1],nums[2]/nums[3])
print((nums[0]+(nums[1]*0.937))/(nums[2]+(nums[3]*1.096)))
#plt.clf()
#plt.plot(xl, yl)
#plt.show()