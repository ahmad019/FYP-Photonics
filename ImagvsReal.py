# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 14:43:18 2018

@author: ahmad
"""


import numpy as np
import matplotlib.pyplot as plt

ran = 62840 #for 2pi limit set to 6284
phi = 0
r = 0.98797
r2 = 0.97889

t = np.sqrt(1-(r**2))
t2 = np.sqrt(1-(r2**2))


Etii = np.ndarray(ran, float)
Etri = np.ndarray(ran, float)
Erri = np.ndarray(ran, float)
Erii = np.ndarray(ran, float)
phit = np.ndarray(ran, float)

for i in range(0,ran):
    
    Era = (r-r2*np.exp(1j*phi))/(1-r*r2*np.exp(1j*phi)) #asymmeteric reflection
    
    Eta = -(t*t2*np.exp(1j*phi/2))/(1-r*r2*np.exp(1j*phi)) #asymmeteric transmission

    
    phi = phi + 0.001
    
    Etri[i] = Eta.real
    Etii[i] = Eta.imag
    Erri[i] = Era.real
    Erii[i] = Era.imag
    phit[i] = phi
    
    

plt.xlim([-1,1])
plt.ylim([-0.7,0.7])


plt.title('Transmision Imaginary vs Real (for limit 0 to 2pi)')
plt.xlabel('Real axis')
plt.ylabel('Imaginary axis')
plt.grid()

plt.plot(Etri,Etii, 'c')
plt.show()

plt.xlim([-0.5,1.2])
plt.ylim([-0.5,0.5])

plt.title('Reflection Imaginary vs Real (for limit 0 to 2pi)')
plt.xlabel('Real axis')
plt.ylabel('Imaginary axis')
plt.grid()

plt.plot(Erri,Erii, 'b')
plt.show()