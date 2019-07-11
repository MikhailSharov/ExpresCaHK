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

plt.rcParams.update({'font.size': 12})

y= []
x= []
SNRH= []
ys = np.array([0.287,0.188,0.155,0.156,0.172,0.177,0.163,0.174,0.156,0.343])
xs = np.array([0.031753941997664556, 0.02962062896181317, 0.025835582783012376, 0.02056515438317584, 0.019247908546459698, 0.02147274290723572, 0.027782226034097503, 0.02120407722984877, 0.020994957078146703, 0.05400447903383693])
xlars = np.array([0.030307364140276804, 0.03182197152813854, 0.03685367716229007, 0.02250276829710482, 0.01950520989659399, 0.02333987182394748, 0.03361299259751715, 0.0264265051297938, 0.023800457727525298, 0.05340623284517232])
SNRHs = np.array([0.0170, 0.0157, 0.0114, 0.0088, 0.0076, 0.0087, 0.0123, 0.0080, 0.0089, 0.0343])
for i in range(len(ys)):
    y.append(ys[i])
    x.append(xs[i])
    SNRH.append(SNRHs[i])
x = np.asarray(x)
y = np.asarray(y)
SNRH = np.asarray(SNRH)
SNRHs = np.array([0.388,0.280,0.229,0.231,0.233,0.321,0.239,0.238,0.309,0.275,0.416,0.445,0.358,0.251,0.298,0.342,0.237,0.358,0.271,0.348,0.274,0.254,0.418,0.668])
SNRK = np.array([0.355,0.264,0.234,0.222,0.225,0.291,0.230,0.223,0.294,0.258,0.288,0.305,0.276,0.240,0.289,0.269,0.226,0.340,0.246,0.281,0.254,0.237,0.368,0.523])
SNRR = np.array([4.022,1.910,1.018,0.924,0.954,3.393,1.582,1.312,3.042,1.370,4.095,4.958,3.928,1.342,3.171,3.805,1.064,4.117,2.102,2.810,2.569,1.679,5.591,7.547])
SNRV = np.array([2.214,1.211,0.930,0.900,0.912,1.888,1.053,0.984,1.741,1.087,2.274,2.865,2.268,1.023,1.826,1.998,0.915,2.383,1.154,1.512,1.410,1.083,3.226,4.457])
#xp = np.array([0.311,0.149,0.162,0.165,0.188,0.137,0.133,0.155,0.393,0.156,0.172,0.177,0.440,0.163,0.174,0.156,0.149])
#yp = np.array([0.0188,0.0216,0.0226,0.0469,0.0165,0.0121,0.0157,0.0183,0.0423,0.0106,0.0083,0.0095,0.0227,0.0199,0.0087,0.0156,0.0146])
#rpsigp = np.array([3.6,0.575,0.013,0.08,2.47,0.27,0.36,1.402,0.302,2.62,3.46,2.48,0.584,2.268,1.915,3.12,0.9])
chis = []
rms = []
##GENERAL EQUATION ------------------------------------------------------------
def func(x, A, c, d):
#    return A*np.exp(c*x) + d
#    return (A*((x-c)**2))+d
    return (A*x) + c
#
# SURVEY ----------------------------------------------------------------------
# Plotting Sampling Data
plt.clf()
plt.scatter(x, y, c = 'k', label="RP Extraction S-values")
plt.scatter(xlars, ys, c = 'orange', label = 'Lars Repack S-values')
####
x_lin = np.linspace(0, xs.max(), 50)                   # 50 evenly spaced digits between 0 and max
###
## Trials
#A, c, d = 20,-0.12,0.12
#y_trial1 = func(x_lin,  A,     c, d)
#y_trial2 = func(x_lin, -1, -1e-3, 1)
#y_trial3 = func(x_lin, -1, -3e-3, 1)
#
#plt.plot(x_lin, y_trial1, "--", label="Trial 1")
####plt.plot(x_lin, y_trial2, "--", label="Trial 2")
####plt.plot(x_lin, y_trial3, "--", label="Trial 3")
#plt.legend()
##
# REGRESSION ------------------------------------------------------------------
#p0 = [20, -0.120, 0.12]                                        # guessed params
#w, _ = opt.curve_fit(func, xs, ys, p0=p0, bounds=([-np.inf,-np.inf,-np.inf], [np.inf,np.inf,0.12]), maxfev = 1000000000)     
#print("Estimated Parameters", w)  
##
## Model
#y_model = func(x_lin, *w)
###
## PLOT ------------------------------------------------------------------------
## Visualize data and fitted curves
#plt.clf()
#plt.plot(xs, ys, "ko", label="Data")
#plt.plot(x_lin, y_model, "k--", label="Fit")
#plt.title("Least squares regression")
plt.legend(loc="upper left")

#for i in range(len(y)):
#    xs[i] = (5.528*xs[i]) + 0.0465
##
#plt.clf()
#plt.scatter(xs,ys,c = 'k',label = 'Calibrated $S_{EXPRES}$ Data')
#plt.plot([0,0.5],[0,0.5],'k--',label='Ideal $S_{MW} = S_{EXPRES}$ Fitting')
##plt.plot([0,0.5],[0.041,0.541],'b--',label='rms')
##plt.plot([0,0.5],[-0.041,0.459],'b--')
#plt.title('$S_{EXPRES}$ Calibrated to $S_{MW}$')
#plt.xlabel('$S_{EXPRES}$')
#plt.ylabel('$S_{MW}$')
#plt.legend(loc='upper left')
##plt.gca().set_yscale('log')
##plt.gca().set_xscale('log')
#plt.plot
#for i in range(len(xlars)):
#    chis.append(((ys[i]-xs[i])**2)/xs[i])
#    rms.append((ys[i]-xs[i])**2)
#print('chi squared = ', sum(chis))
#print('rms = ', str(np.sqrt((sum(rms))/(len(rms)))))

#plt.clf()
#plt.scatter(x,y)