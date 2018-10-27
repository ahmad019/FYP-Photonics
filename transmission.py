# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 21:20:12 2018

@author: ahmad
"""
import numpy as np
import matplotlib.pyplot as plt

ran = 200
phi = -0.8
r = 0.97
r2 = 0.989

t = np.sqrt(1-(r**2))
t2 = np.sqrt(1-(r2**2))


Etai = np.ndarray(ran, float)
Etsi = np.ndarray(ran, float)
phit = np.ndarray(ran, float)

for i in range(0,ran):
    
    Eta = -(t*t2*np.exp(1j*phi/2))/(1-r*r2*np.exp(1j*phi)) #asymmeteric transmission
    #Ets = (t*t)*np.exp(1j*phi)/(1-r*r*np.exp(1j*(phi))) #symmeteric
    
    EtaC = Eta.conjugate()
    #EtsC = Ets.conjugate()
    
    phi = phi + 0.01
    
    Etai[i] = abs(Eta)**2
    #Etsi[i] = abs(Ets)**2
    phit[i] = phi

#plt.xlim([-2,2])
#plt.legend(Eta,Ets, borderpad=2)
#plt.legend([Etai, Etsi], ["line 2", "line 1"])

plt.title('Transmission Curve of Fabry-Perot Resonator')
plt.xlabel('Round Trip Phase')
plt.ylabel('Transmittance')
plt.plot(phit,Etai, 'm') #asymmeteric
#plt.plot(phit,Etsi, 'r') #symmeteric

plt.show()