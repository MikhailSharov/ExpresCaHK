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
x = [0.311,0.149,0.188,0.367,0.137,0.133,0.163,0.155]
y = [0.23116,0.473,0.2304,0.1144,0.5658,0.7630,1.122,0.2743]

ax = plt.gca()

ax.scatter(x,y)
ax.set_yscale('log')
ax.set_xscale('log')
plt.plot