# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 11:45:46 2018

@author: ahmad
"""

import numpy as np
import matplotlib.pyplot as plt
from basic_units import radians, degrees, cos

ran = 10000
phi = -0.5
r  = 0.9999
r2 = 0.9999

t = np.sqrt(1-(r**2))
t2 = np.sqrt(1-(r2**2))


Etai = np.ndarray(ran, float)
Erai = np.ndarray(ran, float)
Ersi = np.ndarray(ran, float)
phit = np.ndarray(ran, float)

for i in range(0,ran):
    
    Era = (r-r2*np.exp(1j*phi))/(1-r*r2*np.exp(1j*phi)) #asymmeteric reflection
    
    Eta = -(t*t2*np.exp(1j*phi/2))/(1-r*r2*np.exp(1j*phi)) #asymmeteric transmission
    
    PHI = np.angle(Era)
    
    phi = phi + 0.0001
    
    Etai[i] = abs(Eta)**2 #abs(Eta*EtaC)
    Erai[i] = abs(Era)**2 #abs(Era*EraC)
    Ersi[i] = PHI
    phit[i] = phi

#plt.xlim([-1,1])
#plt.ylim([-0.2,1.2])
#plt.legend(Eta,Ets, borderpad=2)
#plt.legend([Etai, Etsi], ["line 2", "line 1"])

#phit = [phit[i]*radians for i in range(0,ran)]

plt.title('Transmission')
plt.xlim([-0.1,0.1])
plt.xlabel('Real axis')
plt.ylabel('Imaginary axis')

plt.plot(phit,Erai, 'b') #transmittance
#plt.plot(phit,Erai, 'r') #reflectance

plt.show()

fig, axs = plt.subplots(1)

axs.set_title('Reflection Phase')
axs.set_xlim([-0.05,0.05])
axs.set_xlabel('Round Trip Phase')
axs.set_ylabel('Effective Phase')
axs.plot(phit,Ersi, 'b')

fig.tight_layout()
plt.show()
