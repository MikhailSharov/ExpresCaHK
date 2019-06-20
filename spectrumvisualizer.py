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
for i in ['141004_180602.1051.fits']:
    xrp=[]
    yrp=[]
    image_file = fits.open(i)
    image_file.info()
    image_flux= image_file[0].data.copy()
    image_data =image_file[1].data.copy()
    image_wl = image_data['bary_wavelength']
    image_file.close()
    #print(type(image_wl))
    print(image_wl.shape)
    #image_header = image_file[1].header
    #print(image_header)
    
    for i in range(1000,5000):
        yrp.append(image_flux[6][i])
        
    #uncomment follwing two lines in order to smooth function
    box_kernal = Box1DKernel(6)
    yrp = convolve(yrp, box_kernal)
        
    for i in range(1000,5000):
        xrp.append(image_wl[6][i])
        
    plt.plot(xrp,yrp)
    #plt.plot(image_data['wavelength'][6],image_data['spectrum'][6]/image_data['continuum'][6])


"""
the following block extracts and plots the spectrum that was extracted using Lars' repack pipeline
"""
"""
plt.clf()
for i in ['186408_190524.1127.spec.fits','89744_190317.1073.spec.fits','101501_190518.1122.spec.fits']:
    xl = []
    yl = []
    trix = []
    triy = []
    temp = []
#    #H Core Line
    limage_file = fits.open(i)
    limage_file.info()
    limage_data = limage_file[0].data.copy()
    limage_flux =  limage_data[0,6,:]
    limage_wl =  limage_data[1,6,:]
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
    box_kernal = Box1DKernel(5)
    yl = convolve(yl, box_kernal)
#    #
#    #trix = [3991,3991,4011,4011]
#    #triy = [0,1.2,1.2,0]
#    #
    plt.plot(xl, yl)
#    #plt.plot(trix,triy, c='r')
"""