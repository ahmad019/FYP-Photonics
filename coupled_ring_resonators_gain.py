# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 13:39:29 2018

@author: ahmad
"""


import numpy as np
import matplotlib.pyplot as plt

ran = 2000
phi1 = -1
phi2 = -1
r1 =  0.94565476797
r2 = 0.999189867767
lamd = 0.000001550
Q = 6000000
b = 0.000025
l = 2*np.pi*b
n = 3.45

#aa = (2*np.pi*n)/(Q*lamd)
aa = -1.001
#aa2 = (2*np.pi*n)/(Q*lamd)
aa2 = -1.0001
a1 = np.exp((-aa*l)/2)
a2 = np.exp((-aa2*l)/2)

#r1 = r2*a1 #critical

Etai = np.ndarray(ran, float)
Erai = np.ndarray(ran, float)
PHI = np.ndarray(ran, float)

phi1t = np.ndarray(ran, float)
phi2t = np.ndarray(ran, float)

for i in range(0,ran):
    
    
    Er12 = (r2-a2*np.exp(1j*phi2))/(1-r2*a2*np.exp(1j*phi2)) #coupling r12
    
    Eta = (r1-Er12*a1*np.exp(1j*phi1))/(1-r1*Er12*a1*np.exp(1j*phi1)) #transmission
    
    #Era = (r1-r2*a1*np.exp(1j*phi1))/(1-r1*r2*a1*np.exp(1j*phi1))
    
    
    phi1 = phi1 + 0.001    
    phi2 = phi2 + 0.001    
    
    Etai[i] = abs(Eta)**2 
    #Erai[i] = abs(Era)**2
    #PHI[i] = np.angle(Eta)
    PHI[i] = np.arctan(r1*a1*abs(Er12)*np.sin(phi2+phi1))/(1-r1*a1*abs(Er12)*np.cos(phi2+phi1)) - np.arctan((a1*abs(Er12)*np.sin(phi2+phi1))/(r1-a1*abs(Er12)*np.cos(phi2+phi1)))
    
    phi1t[i] = phi1
    phi2t[i] = phi2
    
    
fig = plt.figure()

#plt.ylim(top=1.02,bottom=0.9)

plt.title("Transmitted field // EIA")
plt.plot(phi1t,Etai, 'b')
plt.show()

#plt.xlim([-0.2,0.2])
plt.title("Phase")
plt.plot(phi2t,PHI, 'r')
plt.show()

#fig.savefig('tempds.png', dpi=200)