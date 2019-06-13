#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:33:34 2019
@author: msharov
"""
import numpy as np
from astropy import units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Box1DKernel
from scipy import *


trix = []
triy = []
nums = []
limage_file = fits.open('190406_190524.1131.spec.fits')
#H Core Line

limage_file.info()
limage_data = limage_file[0].data.copy()    #order import
limage_flux =  limage_data[0,6,:]
limage_wl =  limage_data[1,6,:]
limage_cont =  limage_data[2,6,:]

xl = []                                         #array reset
yl = []
tempx = []
tempy = []

for i in range(1000, 6000):
    yl.append(limage_flux[i]/(limage_cont[i]*1000)) #general region import
    xl.append(limage_wl[i])


for i in range(len(xl)):
    if xl[i] < 3969.545 and xl[i] > 3968.455:   #specific region definition
        tempy.append(yl[i])
        tempx.append(xl[i])

nums.append(np.trapz(tempy, tempx))     #trapezoidal integration

#K Core Line

limage_data = limage_file[0].data.copy()
limage_flux =  limage_data[0,4,:]
limage_wl =  limage_data[1,4,:]
limage_cont =  limage_data[2,4,:]

xl = []
yl = []
tempx = []
tempy = []

for i in range(1000, 6000):
    yl.append(limage_flux[i]/(limage_cont[i]*1000))
    xl.append(limage_wl[i])


for i in range(len(xl)):
    if xl[i] < 3935.545 and xl[i] > 3934.455:
        tempy.append(yl[i])
        tempx.append(xl[i])

nums.append(np.trapz(tempy, tempx))

#R continuum

limage_data = limage_file[0].data.copy()
limage_flux =  limage_data[0,7,:]
limage_wl =  limage_data[1,7,:]
limage_cont =  limage_data[2,7,:]

xl = []
yl = []
tempx= []
tempy =[]

for i in range(1000, 6000):
    yl.append(limage_flux[i]/(limage_cont[i]*1000))
    xl.append(limage_wl[i])


for i in range(len(xl)):
    if xl[i] < 4011 and xl[i] > 3991:
        tempy.append(yl[i])
        tempx.append(xl[i])

nums.append(np.trapz(tempy, tempx))

#V Continuum

limage_data = limage_file[0].data.copy()
limage_flux =  limage_data[0,3,:]
limage_wl =  limage_data[1,3,:]
limage_cont =  limage_data[2,3,:]

xl = []
yl = []
tempx = []
tempy = []

for i in range(1000, 6000):
    yl.append(limage_flux[i]/(limage_cont[i]*1000))
    xl.append(limage_wl[i])


for i in range(len(xl)):
    if xl[i] < 3911 and xl[i] > 3891:
        tempy.append(yl[i])
        tempx.append(xl[i])

nums.append(np.trapz(tempy, tempx))

limage_file.close()

print(nums[0],nums[1],nums[2],nums[3])
print((nums[0]+(nums[1]*1.117))/(nums[2]+(nums[3]*1.832)))