# -*- coding: utf-8 -*-
"""
Created on Sat Dec 15 19:44:25 2018

@author: ahmad
"""


import numpy as np
import matplotlib.pyplot as plt

ran = 6284
phi1 = -np.pi
phi2 = -np.pi
r1 =  0.998897
r2 =  0.979889
lamd = 0.000001550
Q1 = 5000000
Q2 = 6000000
b = 0.000025
l = 2*np.pi*b
n = 3.45
t = np.sqrt(1-r2**2)

aa = (2*np.pi*n)/(Q1*lamd)
aa2 = (2*np.pi*n)/(Q2*lamd)
a1 = np.exp((-aa*l)/2)
a2 = np.exp((-aa2*l)/2)
'''
r2 = a2
r1 = r2*a2 #critical
t = np.sqrt(1-r2**2)
'''
Etai = np.ndarray(ran, float)
Erai = np.ndarray(ran, float)
PHI = np.ndarray(ran, float)
PHIt = np.ndarray(ran, float)


phi1t = np.ndarray(ran, float)
phi2t = np.ndarray(ran, float)

for i in range(0,ran):
    
    
    Er12 = (r2-r1*a2*np.exp(1j*phi2))/(1-r2*a2*np.exp(1j*phi2)) #coupling r12
    
    Eta = (1j*t*Er12*a1*np.exp(1j*phi1))/(1-Er12*a1*np.exp(1j*phi1)) #transmission
    
    
    phi1 = phi1 + 0.001    
    phi2 = phi2 + 0.001    
    
    Etai[i] = abs(Eta)**2
    Erai[i] = abs(Er12)**2
    PHI[i] = np.angle(Eta)
    PHIt[i] = np.angle(Er12)

    phi1t[i] = phi1
    phi2t[i] = phi2


fig = plt.figure()

#plt.xlim([-0.5,0.5])
#plt.title("Intracavity phase")
#plt.xlabel('phase')
#plt.ylabel('Intensity')
plt.plot(phi1t,PHI, 'b')
plt.show()

plt.plot(phi1t,Etai, 'b')
plt.show()

plt.plot(phi2t,Erai, 'r')
plt.show()

plt.title("Intracavity phase")
plt.plot(phi2t,PHIt, 'r')
plt.show()

fig.savefig('temp.png', dpi=200)