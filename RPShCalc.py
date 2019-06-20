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
s = []
avg1 = []
avg2 = []
for j in ['101501_190503.1087.fits','103095_190503.1088.fits','95128_190210.1126.fits','144579_190427.1085.fits','146233_190427.1090.fits','145675_190515.1111.fits']:
    nums = []
    image_file = fits.open(j)
    
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
        xl.append(image_data['bary_wavelength'][6][i])
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
        xl.append(image_data['bary_wavelength'][4][i])
    
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
        xl.append(image_data['bary_wavelength'][7][i])
    
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
        xl.append(image_data['bary_wavelength'][3][i])
    
    for i in range(len(xl)):
        if xl[i] < 3911 and xl[i] > 3891:   #specific region definition
            tempy.append(yl[i])
            tempx.append(xl[i])
    
    nums.append(np.trapz(tempy, tempx))     #trapezoidal integration
    
    
    print(nums[0],nums[1],nums[2],nums[3])
    print(nums[0]/nums[1],nums[2]/nums[3])
    s.append((nums[0]+(nums[1]*0.937))/(nums[2]+(nums[3]*1.096)))
    avg1.append(nums[0]/nums[1])
    avg2.append(nums[2]/nums[3])
print(np.mean(np.asarray(avg1)))
print(np.mean(np.asarray(avg2)))
print(s)
#plt.clf()
#plt.plot(xl, yl)
#plt.show()