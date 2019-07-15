#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 16:07:07 2019

@author: msharov
"""

import numpy as np
from astropy import units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Box1DKernel
plt.rcParams.update({'font.size': 12})
"""
the following block extracts and plots the spectrum that was extracted using Ryan Petersburg's pipeline
"""
temp1 = []

o = 4
b = 1
m = 1000
n = 5100
plt.clf()
SNR = []
for i,j in zip(['101501_190503.1087.fits'],[0.999983]):
    xrp=[]
    yrp=[]
    image_file = fits.open(i)
    image_file.info()
    image_flux= image_file[0].data.copy()
    image_data =image_file[1].data.copy()
    image_wl = image_data['bary_wavelength']
    image_un = image_data['uncertainty']
    image_si = image_data['spectrum']
    signal = []
    uncertainty = []
    
    for k in range(len(image_wl)):
        image_wl[k] = image_wl[k]/j
        
    image_file.close()
    #print(type(image_wl))
    print(image_wl.shape)
    #image_header = image_file[1].header
    #print(image_header)
    
    for k in range(m,n):
        yrp.append((image_flux[o][k]))
    temp1.append(yrp)    
        
    #uncomment follwing two lines in order to smooth function
    box_kernal = Box1DKernel(2)
    yrp = convolve(yrp, box_kernal)
        
    for k in range(m,n):
        xrp.append(image_wl[o][k])
        signal.append(image_si[o][k])
        uncertainty.append(image_un[o][k])
#    for k in range(len(xrp)):
#        if xrp[k] > 3994.89 and xrp[k] < 3994.90:
#            SNR.append([signal[k]/uncertainty[k], i])
    plt.plot(xrp,yrp, c = 'k', label = 'RP')
    #plt.plot(image_data['wavelength'][6],image_data['spectrum'][6]/image_data['continuum'][6])

#print((np.correlate(temp1[1], temp1[0]))/len(temp1[1]))

"""
the following block extracts and plots the spectrum that was extracted using Lars' repack pipeline
"""

"""
for i,j,l in zip(['141004_180602.1050.spec.fits'],['141004_180602.1050.fits'],[1.0000631]):
    xl = []
    yl = []
    trix = []
    triy = []
    temp = []
#    #H Core Line
    limage_file = fits.open(i)
    image_file = fits.open(j)
    limage_file.info()
    image_data =image_file[1].data.copy()
    limage_data = limage_file[0].data.copy()
    limage_flux =  limage_data[0,o,:]
    limage_wl = image_data['bary_wavelength'][o][1:]
    for k in range(len(limage_wl)):
        limage_wl[k] = limage_wl[k]/l
    limage_cont =  limage_data[2,o,:]
    
    xl = []
    yl = []
    temp = []
    
    for i in range(m, n):
        yl.append((limage_flux[i]/(limage_cont[i]*1000))/b)
        xl.append(limage_wl[i])
    
    
    for i in range(len(xl)):
        if xl[i] < 3969.545 and xl[i] > 3968.455:
            temp.append(yl[i])
    
    print(np.mean(temp))
    #K Core Line
#    limage_data = limage_file[0].data.copy()
#    limage_flux =  limage_data[0,4,:]
#    limage_wl =  limage_data[1,4,:]
#    limage_cont =  limage_data[2,4,:]
#    
#    xl = []
#    yl = []
#    temp = []
#    
#    for i in range(1000, 6000):
#        yl.append(limage_flux[i]/(limage_cont[i]*1000))
#        xl.append(limage_wl[i])
#    
#    
#    for i in range(len(xl)):
#        if xl[i] < 3935.545 and xl[i] > 3934.455:
#            temp.append(yl[i])
#    
#    print(np.mean(temp))
#    # R continuum
#    limage_data = limage_file[0].data.copy()
#    limage_flux =  limage_data[0,7,:]
#    limage_wl =  limage_data[1,7,:]
#    limage_cont =  limage_data[2,7,:]
#    
#    xl = []
#    yl = []
#    temp = []
#    
#    for i in range(1000, 6000):
#        yl.append(limage_flux[i]/(limage_cont[i]*1000))
#        xl.append(limage_wl[i])
#    
#    
#    for i in range(len(xl)):
#        if xl[i] < 4011 and xl[i] > 3991:
#            temp.append(yl[i])
#    
#    print(np.mean(temp))
#    
#    #V Continuum
#    limage_data = limage_file[0].data.copy()
#    limage_flux =  limage_data[0,3,:]
#    limage_wl =  limage_data[1,3,:]
#    limage_cont =  limage_data[2,3,:]
#    
#    xl = []
#    yl = []
#    temp = []
#    
#    for i in range(1000, 6000):
#        yl.append(limage_flux[i]/(limage_cont[i]*1000))
#        xl.append(limage_wl[i])
#    
#    
#    for i in range(len(xl)):
#        if xl[i] < 3911 and xl[i] > 3891:
#            temp.append(yl[i])
#    
#    print(np.mean(temp))
#    
#    
#    limage_file.close()
#    #
    box_kernal = Box1DKernel(3)
    yl = convolve(yl, box_kernal)
"""
#    #
trix = [3933.61, 3969.9, 3970.69]
triy = [0, 1.4, 0]
#    #
#plt.plot(xl, yl, label = 'Lars')
#plt.legend()
plt.plot(trix,triy, c='red')