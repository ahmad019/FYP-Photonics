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
r1 =  0.98797
r2 = 0.97889
lamd = 0.000001550
Q = 6000000
b = 0.000025
l = 2*np.pi*b
n = 3.45

aa = (2*np.pi*n)/(Q*lamd)
aa2 = (2*np.pi*n)/(Q*lamd)
a1 = np.exp((-aa*l)/2)
a2 = np.exp((-aa2*l)/2)

#r1 = r2*a1 #critical

Etai = np.ndarray(ran, float)
Erai = np.ndarray(ran, float)


phi1t = np.ndarray(ran, float)
phi2t = np.ndarray(ran, float)

for i in range(0,ran):
    
    
    Er12 = (r2-a2*np.exp(1j*phi2))/(1-r2*a2*np.exp(1j*phi2)) #coupling r12
    
    Eta = (r1-Er12*a1*np.exp(1j*phi1))/(1-r1*Er12*a1*np.exp(1j*phi1)) #transmission
    
    Eta = (r1-r2*a1*np.exp(1j*phi1))/(1-r1*r2*a1*np.exp(1j*phi1))
    
    
    phi1 = phi1 + 0.001    
    phi2 = phi2 + 0.001    
    
    Etai[i] = abs(Eta)**2 
    Erai[i] = abs(Er12)**2

    phi1t[i] = phi1
    phi2t[i] = phi2

plt.plot(phi1t,Etai, 'b')
plt.show()
plt.plot(phi2t,Erai, 'r')
plt.show()
