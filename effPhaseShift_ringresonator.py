# -*- coding: utf-8 -*-
"""
Created on Sun Nov 18 21:50:24 2018

@author: ahmad
"""


import numpy as np
import matplotlib.pyplot as plt

ran = 6284
phi = -np.pi
r = 0.9989
r2 = 0.7541
r3 = 0.5012
r4 = 0.2592
r5 = 0.0


PHIi = np.ndarray(ran, float)
PHI2i = np.ndarray(ran, float)
PHI3i = np.ndarray(ran, float)
PHI4i = np.ndarray(ran, float)
PHI5i = np.ndarray(ran, float)


phit = np.ndarray(ran, float)

for i in range(0,ran):
    
    PHI = np.pi + phi + 2*np.arctan((r*np.sin(phi))/(1-r*np.cos(phi)))
    PHI2 = np.pi + phi + 2*np.arctan((r2*np.sin(phi))/(1-r2*np.cos(phi)))
    PHI3 = np.pi + phi + 2*np.arctan((r3*np.sin(phi))/(1-r3*np.cos(phi)))
    PHI4 = np.pi + phi + 2*np.arctan((r4*np.sin(phi))/(1-r4*np.cos(phi)))
    PHI5 = np.pi + phi + 2*np.arctan((r5*np.sin(phi))/(1-r5*np.cos(phi)))
    
    phi = phi + 0.001
    
    PHIi[i] = PHI
    PHI2i[i] = PHI2 
    PHI3i[i] = PHI3 
    PHI4i[i] = PHI4 
    PHI5i[i] = PHI5 
    
    phit[i] = phi

#plt.xlim([-1,1])
#plt.ylim([0.6,1.05])


plt.title('all pass')
plt.xlabel('single-pass phase')
plt.ylabel('effective phase shift')
plt.plot(phit,PHIi, 'm')
plt.plot(phit,PHI2i, 'orange')
plt.plot(phit,PHI3i, 'purple')
plt.plot(phit,PHI4i, 'r')
plt.plot(phit,PHI5i, 'b')

plt.show()