# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 14:43:18 2018

@author: ahmad
"""


import numpy as np
import matplotlib.pyplot as plt

ran = 200
phi = -1
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

    
    phi = phi + 0.01    
    
    Etai[i] = Eta.real #abs(Eta*EtaC)
    Erai[i] = Eta.imag #abs(Era*EraC)
    Ersi[i] = Etai[i] + Erai[i]
    phit[i] = phi

#plt.xlim([-1,1])
#plt.ylim([-2,2])


plt.title('Imaginary vs Real')
plt.xlabel('Real axis')
plt.ylabel('Imaginary axis')

plt.plot(Etai,Erai, 'b')
plt.show()
