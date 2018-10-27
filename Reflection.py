# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 12:42:27 2018

@author: ahmad
"""

import numpy as np
import matplotlib.pyplot as plt

ran = 200
phi = -1
r = 0.9989
r2 = 0.9889

t = np.sqrt(1-(r**2))
t2 = np.sqrt(1-(r2**2))


Erai = np.ndarray(ran, float)
Ersi = np.ndarray(ran, float)
phit = np.ndarray(ran, float)

for i in range(0,ran):
    
    Era = (r-r2*np.exp(1j*phi))/(1-r*r2*np.exp(1j*phi)) #asymmeteric
    #Ers = (r*r*np.exp(1j*phi))/(1-r*r*np.exp(1j*phi)) #symmeteric
    
    
    EraC = Era.conjugate()
    #ErsC = Ers.conjugate()
    
    phi = phi + 0.01
    
    Erai[i] = abs(Era)**2 #np.power(Era.real,2) #
    #Ersi[i] = abs(Ers)**2 #*ErsC)
    phit[i] = phi

plt.xlim([-1,1])
plt.ylim([0.6,1.05])
#plt.legend(Eta,Ets, borderpad=2)
#plt.legend([Etai, Etsi], ["line 2", "line 1"])

plt.title('Reflection Curve of Fabry-Perot Resonator')
plt.xlabel('Round Trip Phase')
plt.ylabel('Reflectance')
plt.plot(phit,Erai, 'm') #asymmeteric
#plt.plot(phit,Ersi, 'r') #symmeteric

plt.show()