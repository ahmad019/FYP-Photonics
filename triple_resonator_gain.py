# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 15:00:19 2019

@author: ahmad
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 13:39:29 2018

@author: ahmad
"""

import cmath
import numpy as np
import matplotlib.pyplot as plt

ran = 2000
phi1 = -1
phi2 = -1
phi3 = -1
r1 =  0.9465476797
r2 = 0.99189867767
r3 = 0.9798389
lamd = 0.000001550
Q = 6000000
b = 0.000025
l = 2*np.pi*b
n = 3.45

#aa = (2*np.pi*n)/(Q*lamd)
aa = -1.0001
#aa2 = (2*np.pi*n)/(Q*lamd)
aa2 = -1.1001
aa3 = -1.00001
a1 = np.exp((-aa*l)/2)
a2 = np.exp((-aa2*l)/2)
a3 = np.exp((-aa3*l)/2)

#r1 = r2*a1 #critical

Etai = np.ndarray(ran, float)
Erai = np.ndarray(ran, float)
PHI = np.ndarray(ran, float)

phi1t = np.ndarray(ran, float)
phi2t = np.ndarray(ran, float)

for i in range(0,ran):
    
    
    Er23 = (r3-np.exp(1j*phi3)*np.exp(-aa3*l/2))/(1-r3*np.exp(1j*phi3)*np.exp(-aa3*l/2)) #coupling r23
    
    Er12 = (r2-Er23*np.exp(1j*phi2)*np.exp(-aa2*l/2))/(1-r2*Er23*np.exp(1j*phi2)*np.exp(-aa3*l/2)) #coupling r12
    
    Eta = (r1-Er12*a1*np.exp(1j*phi1)*np.exp(-aa*l/2))/(1-r1*Er12*np.exp(1j*phi1)*np.exp(-aa*l/2)) #transmission
    
    #Era = (r1-r2*a1*np.exp(1j*phi1))/(1-r1*r2*a1*np.exp(1j*phi1))
    
    
    phi1 = phi1 + 0.001    
    phi2 = phi2 + 0.001    
    phi3 = phi3 + 0.001
    
    Etai[i] = abs(Eta)**2 
    #Erai[i] = abs(Era)**2
    PHI[i] = np.angle(Eta)

    phi1t[i] = phi1
    phi2t[i] = phi2
    
    
fig = plt.figure()
fig2 = plt.figure()

plt.title("Transmitted field // EIA")
plt.plot(phi1t,Etai, 'b')
plt.show()
fig.savefig('qstemp.png', dpi=200)
plt.title("Phase")
plt.plot(phi2t,PHI, 'r')
plt.show()

#fig2.savefig('qtemp.png', dpi=200)