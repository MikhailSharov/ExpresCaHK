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

"""
the following block extracts and plots the spectrum that was extracted using Ryan Petersburg's pipeline
"""


plt.clf()
for i,j in zip(['101501_190503.1087.fits'],[1.000263, 0.999945, 1.0000631, 1.000467, 1.000357]):
    xrp=[]
    yrp=[]
    image_file = fits.open(i)
    image_file.info()
    image_flux= image_file[0].data.copy()
    image_data =image_file[1].data.copy()
    image_wl = image_data['bary_wavelength']
    for k in range(len(image_wl)):
        image_wl[k] = image_wl[k]/j
        
    image_file.close()
    #print(type(image_wl))
    print(image_wl.shape)
    #image_header = image_file[1].header
    #print(image_header)
    
    for i in range(1000,5000):
        yrp.append(image_flux[3][i])
        
    #uncomment follwing two lines in order to smooth function
    box_kernal = Box1DKernel(3)
    yrp = convolve(yrp, box_kernal)
        
    for i in range(1000,5000):
        xrp.append(image_wl[3][i])
        
    plt.plot(xrp,yrp)
    #plt.plot(image_data['wavelength'][6],image_data['spectrum'][6]/image_data['continuum'][6])


"""
the following block extracts and plots the spectrum that was extracted using Lars' repack pipeline
"""
"""
plt.clf()
for i,j in zip(['141004_180602.1050.spec.fits','86728_190210.1119.spec.fits','103095_190503.1088.spec.fits'], ['141004_180602.1050.fits','86728_190210.1119.fits', '101501_190503.1087.fits']):
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
    limage_flux =  limage_data[0,6,:]
    limage_wl = image_data['bary_wavelength'][6][1:]
    limage_cont =  limage_data[2,6,:]
    
    xl = []
    yl = []
    temp = []
    
    for i in range(1000, 5000):
        yl.append(limage_flux[i]/(limage_cont[i]*1000))
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
#trix = [3933.055,3933.6,3934.145]
#triy = [0,1.2,0]
#    #
#plt.plot(xl, yl)

#plt.plot(trix,triy, c='k')
