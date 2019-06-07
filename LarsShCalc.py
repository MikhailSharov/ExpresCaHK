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


xl = []
yl = []
trix = []
triy = []
temp = []
nums = []
#H Core Line
limage_file = fits.open('157214_180625.1150.spec.fits')
limage_file.info()
limage_data = limage_file[0].data.copy()
limage_flux =  limage_data[0,6,:]
limage_wl =  limage_data[1,6,:]
limage_cont =  limage_data[2,6,:]

xl = []
yl = []
temp = []

for i in range(1000, 6000):
    yl.append(limage_flux[i]/(limage_cont[i]*1000))
    xl.append(limage_wl[i])


for i in range(len(xl)):
    if xl[i] < 3969.545 and xl[i] > 3968.455:
        temp.append(yl[i])

nums.append(np.mean(temp))
#K Core Line
limage_data = limage_file[0].data.copy()
limage_flux =  limage_data[0,4,:]
limage_wl =  limage_data[1,4,:]
limage_cont =  limage_data[2,4,:]

xl = []
yl = []
temp = []

for i in range(1000, 6000):
    yl.append(limage_flux[i]/(limage_cont[i]*1000))
    xl.append(limage_wl[i])


for i in range(len(xl)):
    if xl[i] < 3935.545 and xl[i] > 3934.455:
        temp.append(yl[i])

nums.append(np.mean(temp))
# R continuum
limage_data = limage_file[0].data.copy()
limage_flux =  limage_data[0,7,:]
limage_wl =  limage_data[1,7,:]
limage_cont =  limage_data[2,7,:]

xl = []
yl = []
temp = []

for i in range(1000, 6000):
    yl.append(limage_flux[i]/(limage_cont[i]*1000))
    xl.append(limage_wl[i])


for i in range(len(xl)):
    if xl[i] < 4011 and xl[i] > 3991:
        temp.append(yl[i])

nums.append(np.mean(temp))

#V Continuum
limage_data = limage_file[0].data.copy()
limage_flux =  limage_data[0,3,:]
limage_wl =  limage_data[1,3,:]
limage_cont =  limage_data[2,3,:]

xl = []
yl = []
temp = []

for i in range(1000, 6000):
    yl.append(limage_flux[i]/(limage_cont[i]*1000))
    xl.append(limage_wl[i])


for i in range(len(xl)):
    if xl[i] < 3911 and xl[i] > 3891:
        temp.append(yl[i])

nums.append(np.mean(temp))

limage_file.close()
print(nums[0],nums[1],nums[2],nums[3])

print((nums[0]+nums[1])/(nums[2]+nums[3]))