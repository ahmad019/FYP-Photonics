# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 11:45:46 2018

@author: ahmad
"""

import numpy as np
import matplotlib.pyplot as plt

ran = 200
phi = -0.8
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
    
    
    #EtaC = Eta.conjugate()
    #EraC = Era.conjugate()
    #ErsC = Ers.conjugate()
    
    
    phi = phi + 0.01
    Etai[i] = Eta #abs(Eta*EtaC)
    Erai[i] = Era #abs(Era*EraC)
    #Ersi[i] = abs(Ers*ErsC)
    phit[i] = phi

#plt.xlim([-2,2])
#plt.ylim([0,2])
#plt.legend(Eta,Ets, borderpad=2)
#plt.legend([Etai, Etsi], ["line 2", "line 1"])

plt.title('Reflection Curve of Fabry-Perot Resonator')
plt.xlabel('Round Trip Phase')
plt.ylabel('Reflectance')
#plt.plot(phit,Erai, 'm') #asymmeteric
plt.plot(phit,Etai, 'b') #asymmeteric
plt.plot(phit,Erai, 'r') #symmeteric

plt.show()