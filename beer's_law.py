# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 14:00:35 2018

@author: ahmad
"""

import numpy as np
import matplotlib.pyplot as plt


ran = 200
phi = 0
r  = 0.9989
lamd = 0.000001550
Q = 6000000
b = 0.000025
l = 2*np.pi*b
n = 3.45
g = 1.0001

#aa = (2*np.pi*n)/(Q*lamd)
aa = 1.0001
a = np.exp((-1j*aa*l)/2)

Ii = np.ndarray(ran, float)
I1i = np.ndarray(ran, float)
Io = 2
phit = np.ndarray(ran, float)

for i in range(0,ran):
    
    I  = Io*np.exp(g*phi)
    I1  = Io*np.exp(-aa*phi)
    
    phi = phi + 0.01
    
    Ii[i] = I
    I1i[i] = I1
    
    phit[i] = phi


fig, axs = plt.subplots(2)

axs[0].set_title("Beer's law plot without gain")
axs[0].plot(phit, I1i, 'orange')
axs[1].set_title("Beer's law plot with gain")
axs[1].plot(phit, Ii, 'y')


fig.tight_layout()
plt.show()
#fig.savefig('foo.png', dpi=300, transparent=True)