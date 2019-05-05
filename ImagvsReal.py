# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 14:43:18 2018

@author: ahmad
"""


import numpy as np
import matplotlib.pyplot as plt

ran = 6284 #for 2pi limit set to 6284
phi1 = 0
phi2 = 0
r1 = 0.98
r2 = 0.99

t = np.sqrt(1-(r1**2))
t2 = np.sqrt(1-(r2**2))

a1 = 0.88
a2 = 0.9999


Etii = np.ndarray(ran, float)
Etri = np.ndarray(ran, float)
Erri = np.ndarray(ran, float)
Erii = np.ndarray(ran, float)
phit = np.ndarray(ran, float)

for i in range(0,ran):
    
    Er12 = (r2-a2*np.exp(1j*phi2))/(1-r2*a2*np.exp(1j*phi2))    
    Eta = (r1-Er12*a1*np.exp(1j*phi1))/(1-r1*Er12*a1*np.exp(1j*phi1)) #asymmeteric trans
    
    #Era = (r1-r2*a1*np.exp(1j*phi1))/(1-r1*r2*a1*np.exp(1j*phi1)) #asymmeteric reflec
    
    
    
    phi1 = phi1 + 0.001
    phi2 = phi2 + 0.001    
    
    Etri[i] = Eta.real
    Etii[i] = Eta.imag
    Erri[i] = Er12.real
    Erii[i] = Er12.imag
    phit[i] = phi1
    
    
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8,4))

axs[0].set_xlim([0,1.2])
axs[0].set_ylim([-0.5,0.5])
axs[0].set_title('Transmision Imaginary vs Real')
axs[0].set_xlabel('Real axis')
axs[0].set_ylabel('Imaginary axis')
axs[0].grid()
axs[0].plot(Etri,Etii, 'c')


#axs[1].set_xlim([-0.5,1.2])
#axs[1].set_ylim([-1.0,1.0])
axs[1].set_title('Coupling r12 Imaginary vs Real')
axs[1].set_xlabel('Real axis')
axs[1].set_ylabel('Imaginary axis')
axs[1].grid()
axs[1].plot(Erri,Erii, 'b')

fig.tight_layout()
plt.show()

fig.savefig('coupled_ring_ImagvsReal.png', dpi=400)