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

larsfiles = ['101501_190503.1087.spec.fits','103095_190503.1088.spec.fits','141004_180602.1050.spec.fits','86728_190210.1119.spec.fits','95128_190210.1126.spec.fits','109358_190210.1158.spec.fits','144579_190427.1085.spec.fits','146233_190427.1090.spec.fits','157214_180625.1150.spec.fits','165341_190506.1112.spec.fits']
rpexfiles = ['101501_190503.1087.fits','103095_190503.1088.fits','141004_180602.1050.fits','86728_190210.1119.fits','95128_190210.1126.fits','109358_190210.1158.fits','144579_190427.1085.fits','146233_190427.1090.fits','157214_180625.1150.fits','165341_190506.1112.fits']
adjust = [0.999983, 0.999671, 0.999784, 1.000019, 1.000038, 1.000321, 1.000027, 0.999801, 1.000042, 0.999986]
H = [0.77066,0.40073,0.32554,0.71436,1.06734,0.61444,0.36296,0.56210,0.65065,0.26976]
K = [0.36544,0.20643,0.16917,0.36433,0.58846,0.34679,0.18579,0.26511,0.36165,0.12475]
R = [1.20782,0.64610,0.55024,1.09828,1.54777,0.91930,0.59523,0.90916,0.99681,0.51338]
V = [0.38729,0.22063,0.17349,0.37276,0.57043,0.32449,0.21318,0.25695,0.35578,0.13843]

avg1 = []
avg2 = []
avg3 = []
avg4 = []
slars = []
srpex = []

#calculating the s values from Lars Repack and RP Extraction
for  i,j,k,m,n,o,p in  zip(larsfiles, rpexfiles, adjust, H, K, R, V):
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
        yl.append((image_flux[i]/(image_cont[i]*1000))/m) #general region import
        xl.append(image_wl[i])

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
        yl.append((image_flux[i]/(image_cont[i]*1000))/n)
        xl.append(image_wl[i])


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
        yl.append((image_flux[i]/(image_cont[i]*1000))/o)
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

    #V Continuum-------------------------------------------------------------

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
        yl.append((image_flux[i]/(image_cont[i]*1000))/p)
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
    slars.append(((numslars[0]+(numslars[1]*0.7663168891068424))/(numslars[2]+(numslars[3]*1.052049758691577))))
    srpex.append((numsrpex[0]+(numsrpex[1]*0.795690677293081))/(numsrpex[2]+(numsrpex[3]*1.0499143494669443)))

print(np.mean(avg1))
print(np.mean(avg2))
print(np.mean(avg3))
print(np.mean(avg4))
print(slars)
print(srpex)
