# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 21:20:12 2018

@author: ahmad
"""
import numpy as np
import matplotlib.pyplot as plt

t = 21
t2 = 12
ra = 41
phi = -0.1
r = 11
r2 = 15
Ei = np.ndarray(ra, float)
E = np.ndarray(ra, float)
phit = np.ndarray(ra, float)

for i in range(0,ra):
    
    Et = ((1-t*t2)*np.exp(1j*phi/2))/(1-(r*r2)*np.exp(1j*phi))
    Et1 = Et.conjugate()
    
    phi = phi + 0.005
    Ei[i] = abs(Et*Et1)
    #E[i] = Et.real
    phit[i] = phi

plt.xlim([-2,2])
plt.plot(phit,Ei, 'm')
#plt.plot(phit,Ei, 'm')

plt.show()