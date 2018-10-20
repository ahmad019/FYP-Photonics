# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 21:20:12 2018

@author: ahmad
"""
import numpy as np
import matplotlib.pyplot as plt

t = 21
t2 = 12
ran = 410
phi = -1
r = 11
r2 = 15
Etai = np.ndarray(ran, float)
Etsi = np.ndarray(ran, float)
phit = np.ndarray(ran, float)

for i in range(0,ran):
    
    Eta = ((1-t*t2)*np.exp(1j*phi/2))/(1-(r*r2)*np.exp(1j*phi)) #asymmeteric
    Ets = ((1-t*t)*np.exp(1j*phi/2))/(1-(r*r)*np.exp(1j*phi)) #symmeteric
    EtaC = Eta.conjugate()
    EtsC = Ets.conjugate()
    
    phi = phi + 0.05
    Etai[i] = abs(Eta*EtaC)
    Etsi[i] = abs(Ets*EtsC)
    phit[i] = phi

#plt.xlim([-2,2])
#plt.legend('r=r')
plt.title('Transmission Curve of Fabry-Perot Resonator')
plt.xlabel('Round trip/Resonant Phase')
plt.ylabel('Transmitted Energy (mod squared)')
plt.plot(phit,5.70*Etai, 'm') #asymmeteric
plt.plot(phit,Etsi, 'r') #symmeteric

plt.show()