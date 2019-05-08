# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 11:45:46 2018

@author: ahmad
"""

import numpy as np
import matplotlib.pyplot as plt
from basic_units import radians, degrees, cos

ran = 6284
phi = -np.pi
r = 0.98797
r2 = 0.97889

t = np.sqrt(1-(r**2))
t2 = np.sqrt(1-(r2**2))


Etai = np.ndarray(ran, float)
Erai = np.ndarray(ran, float)
Ersi = np.ndarray(ran, float)
phit = np.ndarray(ran, float)

for i in range(0,ran):
    
    Era = (r-r2*np.exp(1j*phi))/(1-r*r2*np.exp(1j*phi)) #asymmeteric reflection
    
    Eta = -(t*t2*np.exp(1j*phi/2))/(1-r*r2*np.exp(1j*phi)) #asymmeteric transmission

    
    phi = phi + 0.001    
    
    Etai[i] = abs(Eta)**2 #abs(Eta*EtaC)
    Erai[i] = abs(Era)**2 #abs(Era*EraC)
    Ersi[i] = np.angle(Era)
    phit[i] = phi

#plt.xlim([-1,1])
#plt.ylim([-0.2,1.2])
#plt.legend(Eta,Ets, borderpad=2)
#plt.legend([Etai, Etsi], ["line 2", "line 1"])

phit = [phit[i]*radians for i in range(0,ran)]

plt.title('Reflection + Transmission Fabry-Perot Resonator')
plt.xlabel('Round Trip Phase "Φ"')
plt.ylabel('Intensity')

plt.plot(phit,Etai, 'b') #transmittance
plt.plot(phit,Erai, 'r') #reflectance

plt.show()

fig, axs = plt.subplots(1)

axs.set_title('Reflection Phase')
axs.set_xlabel('Round Trip Phase "Φ"')
axs.set_ylabel('Effective Phase')
axs.plot(phit,Ersi, 'r')

fig.tight_layout()
plt.show()
