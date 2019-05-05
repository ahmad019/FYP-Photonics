# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 14:00:35 2018

@author: ahmad
"""

import numpy as np
import matplotlib.pyplot as plt


ran = 1000
phi = 0
r  = 0.9989
lamd = 0.000001550
Q = 6000000
b = 0.000025
l = 2*np.pi*b
n = 3.45
g = 1.0001

a1, a2, a3 = 0.01, 0.1, 1

#aa = (2*np.pi*n)/(Q*lamd)
#aa = 1.0001
#a = np.exp((-1j*aa*l)/2)

Ii = np.ndarray(ran, float)
I1i = np.ndarray(ran, float)
I2i = np.ndarray(ran, float)
Io = 2
phit = np.ndarray(ran, float)

for i in range(0,ran):
    
    I  = Io*np.exp(-a1*phi)
    I1  = Io*np.exp(-a2*phi)
    I2  = Io*np.exp(-a3*phi)
    
    phi = phi + 0.01
    
    Ii[i] = I
    I1i[i] = I1
    I2i[i] = I2
    
    phit[i] = phi


fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(5,8))

axs[0].set_title("Beer's law plot with attenuation = 0.01 cm-1")
axs[0].plot(phit, Ii, 'orange')
axs[1].set_title("Beer's law plot with attenuation = 0.1 cm-1")
axs[1].plot(phit, I1i, 'y')
axs[2].set_title("Beer's law plot with attenuation = 1 cm-1")
axs[2].plot(phit, I2i, 'y')
fig.tight_layout()
plt.show()

fig.savefig('foo.png', dpi=300)#, transparent=True)

