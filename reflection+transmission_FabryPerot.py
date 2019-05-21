# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 11:45:46 2018

@author: ahmad
"""

import numpy as np
import matplotlib.pyplot as plt
from basic_units import radians, degrees, cos

ran = 6284
phi = 0
r = 0.98797
r2 = 0.7889

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
    
    Etai[i] = Eta.real #abs(Eta*EtaC)
    Erai[i] = Eta.imag #abs(Era*EraC)
    Ersi[i] = Era.real
    phit[i] = Era.imag

#plt.xlim([-1,1])
#plt.ylim([-0.2,1.2])
#plt.legend(Eta,Ets, borderpad=2)
#plt.legend([Etai, Etsi], ["line 2", "line 1"])

#phit = [phit[i]*radians for i in range(0,ran)]

plt.title('Transmission')
plt.xlabel('Real axis')
plt.ylabel('Imaginary axis')

plt.plot(Etai,Erai, 'b') #transmittance
#plt.plot(phit,Erai, 'r') #reflectance

plt.show()

fig, axs = plt.subplots(1)

axs.set_title('Reflection')
axs.set_xlabel('Real axis')
axs.set_ylabel('Imaginary axis')
axs.plot(Ersi,phit, 'r')

fig.tight_layout()
plt.show()
