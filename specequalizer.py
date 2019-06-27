#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 13:08:18 2019

@author: mikhailsharov
"""
import numpy as np
from astropy import units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Box1DKernel
from scipy import *


larsfiles = ['101501_190503.1087.spec.fits','103095_190503.1088.spec.fits','141004_180602.1050.spec.fits','86728_190210.1119.spec.fits','95128_190210.1126.spec.fits','109358_190210.1158.spec.fits','144579_190427.1085.spec.fits','146233_190427.1090.spec.fits','157214_180625.1150.spec.fits','165341_190506.1112.spec.fits']
rpexfiles = ['101501_190503.1087.fits','103095_190503.1088.fits','141004_180602.1050.fits','86728_190210.1119.fits','95128_190210.1126.fits','109358_190210.1158.fits','144579_190427.1085.fits','146233_190427.1090.fits','157214_180625.1150.fits','165341_190506.1112.fits']
adjust = [1.000263, 0.999945, 1.0000631, 1.000467, 1.000357, 1.000321, 1.000124, 1.000357, 1.0000353, 1.000280]




b = 0.26976

for  i,j,o in  zip(larsfiles, rpexfiles, adjust):
    for l,m,n in zip([6],[1600],[4000]):
        xrp=[]
        yrp=[]
        image_file = fits.open(j)
        image_flux= image_file[0].data.copy()
        image_data =image_file[1].data.copy()
        image_wl = image_data['bary_wavelength']
        image_un = image_data['uncertainty']
        image_si = image_data['spectrum']
        
        for k in range(len(image_wl)):
            image_wl[k] = image_wl[k]/o
            
        #print(type(image_wl))
        #image_header = image_file[1].header
        #print(image_header)
        
        for k in range(m,n):
            yrp.append((image_flux[l][k])) 
            
        #uncomment follwing two lines in order to smooth function
        box_kernal = Box1DKernel(3)
            
        for k in range(m,n):
            xrp.append(image_wl[l][k])
        
        
        limage_file = fits.open(i)
        image_file = fits.open(j)
        image_data =image_file[1].data.copy()
        limage_data = limage_file[0].data.copy()
        limage_flux =  limage_data[0,l,:]
        limage_wl = image_data['bary_wavelength'][l][1:]
        for k in range(len(limage_wl)):
            limage_wl[k] = limage_wl[k]/o
        limage_cont =  limage_data[2,l,:]
        image_file.close()
        limage_file.close()
        
        
        xl = []
        yl = []
        
        for i in range(m, n):
            yl.append((limage_flux[i]/(limage_cont[i]*1000))/b)
            xl.append(limage_wl[i])
            
        print(j)
        print(sum(yrp) - sum(yl))
            
        
