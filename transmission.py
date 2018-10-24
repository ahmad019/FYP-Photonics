# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 21:20:12 2018

@author: ahmad
"""
import numpy as np
import matplotlib.pyplot as plt

t = 0.1
#t2 = 0.9
ran = 410
phi = -1
r = 0.95
#r2 = 0.1
Etai = np.ndarray(ran, float)
Etsi = np.ndarray(ran, float)
phit = np.ndarray(ran, float)

for i in range(0,ran):
    
   # Eta = ((t*t2)*np.exp(1j*phi/2))/(1-(r*r2)*np.exp(1j*phi)) #asymmeteric
    Ets = ((t*t)*np.exp(1j*phi/2))/(1-(r*r)*np.exp(1j*phi)) #symmeteric
    #*Ets = ((-t*t)*np.exp(1j*phi/2))/(1-(r*r)) #actualequation    
    
    #EtaC = Eta.conjugate()
    EtsC = Ets.conjugate()
    
    phi = phi + 0.05
    #Etai[i] = abs(Eta*EtaC)
    Etsi[i] = abs(Ets*EtsC)
    phit[i] = phi

#plt.xlim([-2,2])
#plt.legend(Eta,Ets, borderpad=2)
#plt.legend([Etai, Etsi], ["line 2", "line 1"])

plt.title('Transmission Curve of Fabry-Perot Resonator')
plt.xlabel('Round Trip Phase')
plt.ylabel('Transmittance')
#plt.plot(phit,Etai, 'm') #asymmeteric
plt.plot(phit,Etsi, 'r') #symmeteric

plt.show()