# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 13:39:29 2018

@author: ahmad
"""


import numpy as np
import matplotlib.pyplot as plt

ran = 6284
phi1 = -np.pi
phi2 = -np.pi
#r1 =  0.90
r2 = 0.999
lamd = 0.000001550
Q1 = 5000000
Q2 = 7000000
c = 299792458
b = 0.000025
l = 2*np.pi*b
n = 3.45

g1 = 0*0.091
g2 = 0.120

aa = (2*np.pi*n)/(Q1*lamd)
aa2 = (2*np.pi*n)/(Q2*lamd)
a1 = np.exp((g1-aa*l)/2)
a2 = np.exp((g2-aa2*l)/2)

#a1 = 0.88
#a2 = 1

#r2 = a2
r1 = r2*a1 #critical

#a1 = a1.round(decimals=2)
#a2 = a2.round(decimals=2)

Etai = np.ndarray(ran, float)
dPHIt = np.ndarray(ran, float)
PHIt = np.ndarray(ran, float)
PHI12 = np.ndarray(ran, float)

phi1t = np.ndarray(ran, float)
phi2t = np.ndarray(ran, float)

with open('Coupled_resonator_phi.phase.trans.csv', 'w') as fp:

    for i in range(0,ran):
        
        
        Er12 = (r2-a2*np.exp(1j*phi2))/(1-r2*a2*np.exp(1j*phi2)) #coupling r12
        
        Eta = (r1-Er12*a1*np.exp(1j*phi1))/(1-r1*Er12*a1*np.exp(1j*phi1)) #transmission
        
        #Era = (r1-r2*a1*np.exp(1j*phi1))/(1-r1*r2*a1*np.exp(1j*phi1))
        
        phi12 =  np.arctan((r2*a2*np.sin(phi2))/(1-r2*a2*np.cos(phi2))) - np.arctan((a2*np.sin(phi2))/(r2-a2*np.cos(phi2)))
        
        #PHIt[i] = np.arctan((a1*abs(Er12)*np.sin(phi12+phi1))/(r1-a1*abs(Er12)*np.cos(phi12+phi1))) - np.arctan((r1*a1*abs(Er12)*np.sin(phi12+phi1))/(1-r1*a1*abs(Er12)*np.cos(phi12+phi1)))
        #PHIt[i] = phi1 - (np.arctan((a2*np.sin(phi2))/(r2-a2*np.cos(phi2))) - np.arctan((a2*r2*np.sin(phi2))/(1-a2*r2*np.cos(phi2))))
        PHIt[i] = np.angle(Eta)
        
        PHI12[i] = np.angle(Er12)
        
        phi1 = phi1 + 0.001    
        phi2 = phi2 + 0.001    
        
        Etai[i] = abs(Eta)**2 
        
        dPHIt[i] = ((1-r1*abs(Er12)*a1*np.cos(phi12+phi1))**2/((1-r1*abs(Er12)*a1*np.cos(phi12+phi1))**2 + (r1*abs(Er12)*a1*np.sin(phi1+phi12))**2))  - ((r1-a1*abs(Er12)*np.cos(phi12+phi1))**2 / ((r1-a1*abs(Er12)*np.cos(phi12+phi1))**2  + (abs(Er12)*a1*np.sin(phi12+phi1))**2))
        
        fp.write("{},{},{}\n".format(phi1,PHIt[i], Etai[i]))
    
        phi1t[i] = phi1
        phi2t[i] = phi2
    

#dPHIr = ((1-r1*np.sin(phi))**2/((1-r1*np.sin(phi))**2 + (np.sin(phi))**2))  + ((r1-np.cos(phi))**2 / ((r1-np.cos(phi))**2  + (np.sin(phi))**2)) #singapore paper phase derivt

phi1t = phi1t.round(decimals=4)

vgt = (1/dPHIt)*c/n
ngt = dPHIt * n

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,5))
    
axs[0,0].set_title("Transmitted field")
#axs[0,0].set_xlim([-1,1])
axs[0,0].set_xlabel("detuning")
axs[0,0].set_ylabel("Transmittance")
axs[1,0].text(0.001,0.001,"Gain1 =%.2f" %g1 + "\nGain2 =%.2f" %g2 + "\nCritical coupled\nr1=%.2f" %r1 + "\nr2=%.2f"%r2,fontsize=8)
axs[0,0].plot(phi1t,Etai, 'b')

#axs[0,1].set_xlim([-0.3,0.3])
axs[0,1].set_title("Effective Phase")
axs[0,1].set_xlabel("detuning")
axs[0,1].set_ylabel("Effective phase")
axs[0,1].plot(phi1t,PHIt, 'r')

#axs[1,0].set_xlim([-0.01,0.01])
axs[1,0].set_title("Phi12 Phase")
axs[1,0].set_xlabel("detuning")
axs[1,0].set_ylabel("Effective phase")
axs[1,0].plot(phi1t,PHI12, 'y')

#axs[1,1].set_xlim([-0.25,0.25])
axs[1,1].set_title("Group Index")
axs[1,1].set_xlabel("detuning")
axs[1,1].set_ylabel("Effective phase")
axs[1,1].plot(phi1t,c/vgt, 'g')

fig.tight_layout()
plt.show()

fig.savefig("coupled_ring_EIA.png",dpi=400)
