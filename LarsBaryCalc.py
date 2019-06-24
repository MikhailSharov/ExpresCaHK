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

larsfiles = ['101501_190503.1087.spec.fits','103095_190503.1088.spec.fits','141004_180602.1050.spec.fits','86728_190210.1119.spec.fits','95128_190210.1126.spec.fits','109358_190210.1158.spec.fits','144579_190427.1085.spec.fits','146233_190427.1090.spec.fits','157214_180625.1150.spec.fits','165341_190506.1112.spec.fits','145675_190515.1111.spec.fits']
rpexfiles = ['101501_190503.1087.fits','103095_190503.1088.fits','141004_180602.1050.fits','86728_190210.1119.fits','95128_190210.1126.fits','109358_190210.1158.fits','144579_190427.1085.fits','146233_190427.1090.fits','157214_180625.1150.fits','165341_190506.1112.fits','145675_190515.1111.fits']
adjust = [1.000263,0.999945, 1.0000631, 1.000467, 1.000357, 1.000321, 1.000124, 1.000357, 1.0000353, 1.000280, 1.000245]

avg1 = []
avg2 = []
avg3 = []
avg4 = []
slars = []
srpex = []

#calculating the s values from Lars Repack and RP Extraction
for  i,j,k in  zip(larsfiles, rpexfiles, adjust):
    numslars = []
    numsrpex = []
    specimage_file = fits.open(i)
    waveimage_file = fits.open(j)
    noise = []
    signal = []
    signoi = []
    #H Core Line-----------------------------------------------------------------

    #Lars Repack H Core
    specimage_data = specimage_file[0].data.copy()    #order import
    waveimage_data = waveimage_file[1].data.copy()
    image_flux =  specimage_data[0,6,:]
    image_wl =  waveimage_data['bary_wavelength'][6][1:]
    for l in range(len(image_wl)):
        image_wl[l] = image_wl[l]/k
        
    image_cont =  specimage_data[2,6,:]

    xl = []                                         #array reset
    yl = []
    tempx = []
    tempy = []

    for i in range(1000, 6000):
        yl.append(image_flux[i]/(image_cont[i]*1000)) #general region import
        xl.append(image_wl[i])


    for i in range(len(xl)):
        if xl[i] < 3969.045 and xl[i] > 3967.955:   #specific region definition
            tempy.append(yl[i])
            tempx.append(xl[i])

    numslars.append(np.trapz(tempy, tempx))     #trapezoidal integration
    #RP Extraction H Core
    xl = []                                         #array reset
    yl = []
    tempx = []
    tempy = []

    image_data = waveimage_file[1].data.copy()    #order import

    for i in range(0, 6500):
        yl.append(image_data['spectrum'][6][i]/image_data['continuum'][6][i]) #general region import
        xl.append(image_data['bary_wavelength'][6][i])
        noise.append(image_data['uncertainty'][6][i])
        signal.append(image_data['spectrum'][6][i])

    for l in range(len(xl)):
            xl[l] = xl[l]/k
    for i in range(len(xl)):
        if xl[i] < 3969.045 and xl[i] > 3967.955:   #specific region definition
            tempy.append(yl[i])
            tempx.append(xl[i])

    numsrpex.append(np.trapz(tempy, tempx))     #trapezoidal integration
   
    #K Core Line---------------------------------------------------------------

    specimage_data = specimage_file[0].data.copy()    #order import
    waveimage_data = waveimage_file[1].data.copy()
    image_flux =  specimage_data[0,4,:]
    image_wl =  waveimage_data['bary_wavelength'][4][1:]
    for l in range(len(image_wl)):
        image_wl[l] = image_wl[l]/k
    image_cont =  specimage_data[2,4,:]

    xl = []
    yl = []
    tempx = []
    tempy = []

    for i in range(1000, 6000):
        yl.append(image_flux[i]/(image_cont[i]*1000))
        xl.append(image_wl[i])


    for i in range(len(xl)):
        if xl[i] < 3934.145 and xl[i] > 3933.055:
            tempy.append(yl[i])
            tempx.append(xl[i])

    numslars.append(np.trapz(tempy, tempx))

    xl = []                                         #array reset
    yl = []
    tempx = []
    tempy = []

    for i in range(1000, 6500):
        yl.append(image_data['spectrum'][4][i]/image_data['continuum'][4][i]) #general region import
        xl.append(image_data['bary_wavelength'][4][i])

    for l in range(len(xl)):
            xl[l] = xl[l]/k
    for i in range(len(xl)):
        if xl[i] < 3934.145 and xl[i] > 3933.055:   #specific region definition
            tempy.append(yl[i])
            tempx.append(xl[i])

    numsrpex.append(np.trapz(tempy, tempx))     #trapezoidal integration

    #R continuum------------------------------------------------------------------

    specimage_data = specimage_file[0].data.copy()    #order import
    waveimage_data = waveimage_file[1].data.copy()
    image_flux =  specimage_data[0,7,:]
    image_wl =  waveimage_data['bary_wavelength'][7][1:]
    for l in range(len(image_wl)):
        image_wl[l] = image_wl[l]/k
    image_cont =  specimage_data[2,7,:]

    xl = []
    yl = []
    tempx= []
    tempy =[]

    for i in range(1000, 6000):
        yl.append(image_flux[i]/(image_cont[i]*1000))
        xl.append(image_wl[i])


    for i in range(len(xl)):
        if xl[i] < 4011 and xl[i] > 3991:
            tempy.append(yl[i])
            tempx.append(xl[i])

    numslars.append(np.trapz(tempy, tempx))

    xl = []                                         #array reset
    yl = []
    tempx = []
    tempy = []

    for i in range(0, 6500):
        yl.append(image_data['spectrum'][7][i]/image_data['continuum'][7][i]) #general region import
        xl.append(image_data['bary_wavelength'][7][i])

    for l in range(len(xl)):
            xl[l] = xl[l]/k
    for i in range(len(xl)):
        if xl[i] < 4011 and xl[i] > 3991:   #specific region definition
            tempy.append(yl[i])
            tempx.append(xl[i])

    numsrpex.append(np.trapz(tempy, tempx))     #trapezoidal integration

    #V Continuum

    specimage_data = specimage_file[0].data.copy()    #order import
    waveimage_data = waveimage_file[1].data.copy()
    image_flux =  specimage_data[0,3,:]
    image_wl =  waveimage_data['bary_wavelength'][3][1:]
    for l in range(len(image_wl)):
        image_wl[l] = image_wl[l]/k
    image_cont =  specimage_data[2,3,:]

    xl = []
    yl = []
    tempx = []
    tempy = []

    for i in range(1000, 6000):
        yl.append(image_flux[i]/(image_cont[i]*1000))
        xl.append(image_wl[i])


    for i in range(len(xl)):
        if xl[i] < 3911 and xl[i] > 3891:
            tempy.append(yl[i])
            tempx.append(xl[i])

    numslars.append(np.trapz(tempy, tempx))

    xl = []                                         #array reset
    yl = []
    tempx = []
    tempy = []

    for i in range(0, 6500):
        yl.append(image_data['spectrum'][3][i]/image_data['continuum'][3][i]) #general region import
        xl.append(image_data['bary_wavelength'][3][i])

    for l in range(len(xl)):
            xl[l] = xl[l]/k
    for i in range(len(xl)):
        if xl[i] < 3911 and xl[i] > 3891:   #specific region definition
            tempy.append(yl[i])
            tempx.append(xl[i])

    numsrpex.append(np.trapz(tempy, tempx))     #trapezoidal integration

    specimage_file.close()

    print(numslars[0],numslars[1],numslars[2],numslars[3])
    avg1.append(numslars[0]/numslars[1])
    avg2.append(numslars[2]/numslars[3])
    print(numsrpex[0],numsrpex[1],numsrpex[2],numsrpex[3])
    avg3.append(numsrpex[0]/numsrpex[1])
    avg4.append(numsrpex[2]/numsrpex[3])
    slars.append(((numslars[0]+(numslars[1]*1.374))/(numslars[2]+(numslars[3]*3.198))))
    srpex.append((numsrpex[0]+(numsrpex[1]*0.725))/(numsrpex[2]+(numsrpex[3]*1.030)))
print(np.mean(avg1))
print(np.mean(avg2))
print(np.mean(avg3))
print(np.mean(avg4))
print(slars)
print(srpex)
