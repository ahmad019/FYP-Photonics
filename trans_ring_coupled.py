# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 15:08:10 2018

@author: ahmad
"""

import numpy as np
import matplotlib.pyplot as plt

ran = 6284
phi = -np.pi
r  = 0.9989
r1 = 0.7989
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
    
    #PHI = ((1-r**2)*(1-r1**2)*a)/(1-(2*r*r1*a*np.cos(phi))+(r*r1*a)**2)
    PHI = (r1**2*a**2)-(2*r*r1*a*np.cos(phi)+r**2)/(1-(2*r*r1*a*np.cos(phi))+(r*r1*a)**2)
    
    phi = phi + 0.001
    
    PHIi[i] = PHI
    
    phit[i] = phi

#plt.xlim([-0.01,0.01])
#plt.ylim([0.6,1.05])

plt.title('All-pass ring resonator')
plt.xlabel('Single-pass phase')
plt.ylabel('Effective phase shift')
plt.plot(phit,PHIi, 'm')

plt.show()