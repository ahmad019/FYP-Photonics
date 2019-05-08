# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 21:08:47 2018

@author: ahmad
"""

import numpy as np
import matplotlib.pyplot as plt

ran = 2000
phi = -1
r = 0.98797
r2 = 0.97889
lamd = 0.000001550
Q = 6000000
b = 0.000025
l = 2*np.pi*b
n = 3.45

aa = (2*np.pi*n)/(Q*lamd)
aag = 500
a = np.exp((-aa*l)/2)
ag = np.exp(-(aa-aag)*l/2)

t = np.sqrt(1-(r**2))
t2 = np.sqrt(1-(r2**2))


Etagi = np.ndarray(ran, float)
Eragi = np.ndarray(ran, float)
Erai = np.ndarray(ran, float)
Etai = np.ndarray(ran, float)
phit = np.ndarray(ran, float)

for i in range(0,ran):
    
    Etag = (r-r2*ag*np.exp(1j*phi))/(1-r*r2*ag*np.exp(1j*phi)) #asymmeteric transmission
    
    Erag = -(t*t2*ag*np.exp(1j*phi/2))/(1-r*r2*ag*np.exp(1j*phi)) #asymmeteric reflection
    
    Eta = (r-r2*a*np.exp(1j*phi))/(1-r*r2*a*np.exp(1j*phi))
    
    Era = -(t*t2*a*np.exp(1j*phi/2))/(1-r*r2*a*np.exp(1j*phi))
    
    phi = phi + 0.001    
    
    Etagi[i] = abs(Etag)**2 
    Eragi[i] = abs(Erag)**2
    Etai[i] = abs(Eta)**2
    Erai[i] = abs(Era)**2
    phit[i] = phi



fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10,8))

axs[0, 0].set_title("Transmission without gain")
axs[0, 0].plot(phit, Etai, 'orange')
axs[0, 1].set_title("Reflection without gain")
axs[0, 1].plot(phit, Erai, 'y')
axs[1, 0].set_title("Transmission with gain")
axs[1, 0].plot(phit, Etagi, 'red')
axs[1, 1].set_title("Reflection with gain")
axs[1, 1].plot(phit, Etagi, 'purple')

fig.tight_layout()
plt.show()

fig.savefig('all-pall_gain.png', dpi=400, transparent=True)
