#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:00:16 2019

@author: msharov
"""
import numpy as np
from astropy import units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Box1DKernel
from scipy import *
from scipy.optimize import leastsq
from math import *
import scipy.optimize as opt


x = np.array([0.311,0.149,0.162,0.367,0.165,0.188,0.137,0.133,0.155,0.393,0.156,0.172,0.177,0.440,0.163,0.174,0.238,0.156,0.293,0.343,0.149,0.194,0.148,0.165])
y = np.array([0.0109,0.0246,0.0432,0.0461,0.0455,0.0109,0.0256,0.0336,0.0130,0.037,0.0101,0.0076,0.0086,0.034,0.0113,0.0091,0.0428,0.0101,0.0206,0.0173,0.0144,0.0256,0.0065,0.008])
chis = []
#GENERAL EQUATION ------------------------------------------------------------
def func(x, A, c, d):
    return (A*((x-c)**2))+d
#
# SURVEY ----------------------------------------------------------------------
# Plotting Sampling Data
#plt.clf()
#plt.plot(x, y, "ko", label="Data")
###
x_lin = np.linspace(0, x.max(), 50)                   # 50 evenly spaced digits between 0 and max
##
## Trials
#A, c, d = 2.3, 0.25, 0.001
#y_trial1 = func(x_lin,  A,     c, d)
#y_trial2 = func(x_lin, -1, -1e-3, 1)
#y_trial3 = func(x_lin, -1, -3e-3, 1)
#
#plt.plot(x_lin, y_trial1, "--", label="Trial 1")
#plt.plot(x_lin, y_trial2, "--", label="Trial 2")
#plt.plot(x_lin, y_trial3, "--", label="Trial 3")
#plt.legend()
#
## REGRESSION ------------------------------------------------------------------
#p0 = [2.3, 0.25, 0.001]                                        # guessed params
#w, _ = opt.curve_fit(func, x, y, p0=p0, bounds=(-np.inf, [7.,np.inf,np.inf]), maxfev = 10000000)     
#print("Estimated Parameters", w)  
#
## Model
#y_model = func(x_lin, *w)
#
## PLOT ------------------------------------------------------------------------
## Visualize data and fitted curves
#plt.clf()
#plt.plot(x, y, "ko", label="Data")
#plt.plot(x_lin, y_model, "k--", label="Fit")
#plt.title("Least squares regression")
#plt.legend(loc="upper left")

for i in range(len(y)):
    if x[i] > 0.2404:
        y[i] = 1*np.sqrt((y[i]-0.0002)/(2.767))+0.2404
    if x[i] < 0.2404:
        y[i] = -1*np.sqrt((y[i]-0.0002)/(2.767))+0.2404

plt.clf()
plt.scatter(x,y,c='k',label = 'Fitted Data')
plt.plot([0,0.4],[0,0.4],'k--',label='x=y')
plt.title('$S_{Expres}$ Fitted With Least Squares')
plt.xlabel('$S_{MW}$')
plt.ylabel('$S_{Expres}$')
plt.legend(loc='upper left')
#plt.gca().set_yscale('log')
#plt.gca().set_xscale('log')
plt.plot
for i in range(len(x)):
    chis.append(((y[i]-x[i])**2)/x[i])
print(sum(chis))

#plt.clf()
#plt.scatter(x,y)