# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 12:42:27 2018

@author: ahmad
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 21:20:12 2018

@author: ahmad
"""
import numpy as np
import matplotlib.pyplot as plt

t = 0.1
#t2 = 12
ran = 410
phi = -1
r = 0.9
#r2 = 0.9
Erai = np.ndarray(ran, float)
Ersi = np.ndarray(ran, float)
phit = np.ndarray(ran, float)

for i in range(0,ran):
    
   # Era = ((1-t*t2)*np.exp(1j*phi/2))/(1-(r*r2)*np.exp(1j*phi)) #asymmeteric
    Ers = ((t*t)*np.exp(1j*phi)*r)/(1-(np.power(r,4)*np.exp(1j*(2*phi)))) #symmeteric
    #EraC = Era.conjugate()
    ErsC = Ers.conjugate()
    
    phi = phi + 0.05
   # Erai[i] = abs(Era*EraC)
    Ersi[i] = abs(Ers*ErsC)
    phit[i] = phi

#plt.xlim([-2,2])
#plt.ylim([0,2])
#plt.legend(Eta,Ets, borderpad=2)
#plt.legend([Etai, Etsi], ["line 2", "line 1"])

plt.title('Reflection Curve of Fabry-Perot Resonator')
plt.xlabel('Round Trip Phase')
plt.ylabel('Reflectance')
#plt.plot(phit,5.70*Erai, 'm') #asymmeteric
plt.plot(phit,Ersi, 'r') #symmeteric

plt.show()