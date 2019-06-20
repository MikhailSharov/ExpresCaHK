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

files = ['101501_190518.1122.spec.fits','217014_180602.1082.spec.fits','115617_190515.1077.spec.fits','72905_190317.1090.spec.fits','126053_190316.1092.spec.fits','103095_190503.1088.spec.fits','89744_190317.1073.spec.fits','120136_190317.1097.spec.fits','141004_180602.1050.spec.fits','152391_180602.1056.spec.fits','86728_190210.1119.spec.fits','95128_190210.1126.spec.fits','109358_190210.1158.spec.fits','95735_190419.1119.spec.fits','144579_190427.1085.spec.fits','146233_190427.1090.spec.fits','114783_190503.1098.spec.fits','157214_180625.1150.spec.fits','149661_190505.1113.spec.fits','165341_190506.1112.spec.fits','145675_190515.1111.spec.fits','185144_190524.1126.spec.fits','186408_190524.1127.spec.fits','190406_190524.1131.spec.fits']
trix = []
triy = []
avg1 = []
avg2 = []
s= []
for  i in files:
    nums = []
    limage_file = fits.open(i)
    noise = []
    signal = []
    signoi = []
    #H Core Line
    
    limage_file.info()
    limage_data = limage_file[0].data.copy()    #order import
    limage_flux =  limage_data[0,6,:]
    limage_wl =  limage_data[1,6,:]
    limage_cont =  limage_data[2,6,:]
    limage_unct = limage_data[3,6,:]
    
    xl = []                                         #array reset
    yl = []
    tempx = []
    tempy = []
    
    for i in range(1000, 6000):
        yl.append(limage_flux[i]/(limage_cont[i]*1000)) #general region import
        xl.append(limage_wl[i])
    #    noise.append(limage_unct[i])
    #    signal.append(limage_flux[i])
    
    
    for i in range(len(xl)):
        if xl[i] < 3969.545 and xl[i] > 3968.455:   #specific region definition
            tempy.append(yl[i])
            tempx.append(xl[i])
    #        signoi.append(signal[i]/noise[i])
    
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
    
    print((nums[0]),(nums[1]),(nums[2]),(nums[3]))
    avg1.append(nums[0]/nums[1])
    avg2.append(nums[2]/nums[3])
    s.append(((nums[0]+(nums[1]*1.296))/(nums[2]+(nums[3]*2.557))))
print(np.mean(np.asarray(avg1)), np.mean(np.asarray(avg2)))
print(s)

#ys = np.array([0.311,0.149,0.162,0.367,0.165,0.188,0.137,0.133,0.155,0.393,0.156,0.172,0.177,0.440,0.163,0.174,0.238,0.156,0.293,0.343,0.149,0.194,0.148,0.165])
#xs = np.array([0.0103, 0.0228, 0.0379, 0.0398, 0.0395, 0.0103, 0.0235, 0.0288, 0.0122, 0.0330, 0.0093, 0.0070, 0.0079, 0.0307, 0.0106, 0.0085, 0.0374, 0.0088, 0.0194, 0.0162, 0.0136, 0.0236, 0.0061, 0.0074])