# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 14:00:35 2018

@author: ahmad
"""

import numpy as np
import matplotlib.pyplot as plt

from basic_units import radians

ran = 6284
phi = -np.pi
r  = 0.9989
lamd = 0.000001550
Q = 6000000
b = 0.000025
l = 2*np.pi*b
n = 3.45
g = 1.0001

aa = (2*np.pi*n)/(Q*lamd)
a = np.exp((-1j*aa*l)/2)

Ii = np.ndarray(ran, float)
I1i = np.ndarray(ran, float)
Io = 2
phit = np.ndarray(ran, float)

for i in range(0,ran):
    
    I  = Io*np.exp(g*phi)
    I1  = Io*np.exp(-aa*phi)
    
    phi = phi + 0.001
    
    Ii[i] = I
    I1i[i] = I1
    
    phit[i] = phi


fig, axs = plt.subplots(2)

axs[0].plot(phit, Ii, xunits=radians)
axs[1].plot(phit, I1i, xunits=radians)

fig.tight_layout()
plt.show()