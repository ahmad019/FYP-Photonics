# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 14:18:16 2018

@author: ahmad
"""

import numpy as np
import matplotlib.pyplot as plt

ran = 6284
phi = -np.pi
r  = 0.9989
lamd = 0.000001550
Q = 6000000
b = 0.000025
l = 2*np.pi*b
n = 3.45

aa = (2*np.pi*n)/(Q*lamd)
a = np.exp((-1j*aa*l)/2)

PHIi = np.ndarray(ran, float)
phit = np.ndarray(ran, float)

for i in range(0,ran):
    
    PHI  = np.pi + phi + np.arctan((r*np.sin(phi))/(a-r*np.cos(phi))) + np.arctan((r*a*np.sin(phi))/(1-r*a*np.cos(phi)))
    
    phi = phi + 0.001
    
    PHIi[i] = abs(PHI)**2
    
    phit[i] = phi

#plt.xlim([-0.01,0.01])
#plt.ylim([0.6,1.05])

plt.title('All-pass ring resonator')
plt.xlabel('Single-pass phase')
plt.ylabel('Effective phase shift')
plt.scatter(phit,PHIi, color='m')

plt.show()