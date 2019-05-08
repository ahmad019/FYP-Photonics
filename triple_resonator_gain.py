# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 15:00:19 2019

@author: ahmad
"""


import numpy as np
import matplotlib.pyplot as plt
from basic_units import radians, degrees, cos


ran = 6284
phi1 = -np.pi
phi2 = -np.pi
phi3 = -np.pi
r1 = 0.8999
r2 = 0.9889
r3 = 0.9989
lamd = 1550e-9
Q1 = 1e5
Q2 = 1e6
Q3 = 1e7
b = 10e-6
l = 2*np.pi*b
n = 3.45
c = 299792458

tau = (n*2*np.pi*b)/c

aa = (2*np.pi*n)/(Q1*lamd)
g1 = 0
aa2 = (2*np.pi*n)/(Q2*lamd)
g2 = 0
aa3 = (2*np.pi*n)/(Q3*lamd)
g3 = 0
a1 = np.exp(-(aa-g1)*l/2)
a2 = np.exp(-(aa2-g2)*l/2)
a3 = np.exp(-(aa3-g3)*l/2)

g=g1+g2+g3
#r2 = r3*a2
#r1 = r2*a1 #critical


Etai = np.ndarray(ran, float)
Erai = np.ndarray(ran, float)
PHIt = np.ndarray(ran, float)
EtaImag = np.ndarray(ran, float)
EtaReal = np.ndarray(ran, float)
vgt = np.ndarray(ran, float)

phi1t = np.ndarray(ran, float)
phi2t = np.ndarray(ran, float)

with open('Triple_resonator_phi.phase.trans.csv', 'w') as fp:

    for i in range(0,ran):
    
    
        Er23 = (r3-a3*np.exp(1j*phi3))/(1-r3*a3*np.exp(1j*phi3)) #coupling r23
    
        Er12 = (r2-a2*Er23*np.exp(1j*phi2))/(1-r2*Er23*a2*np.exp(1j*phi2)) #coupling r12
        
        Eta = (r1-Er12*a1*np.exp(1j*phi1))/(1-r1*Er12*a1*np.exp(1j*phi1)) #transmission
        
        #phi23 = np.pi + phi3 + 2*np.arctan((r3*np.sin(phi3))/(1-r3*np.cos(phi3))) #phase of r23
        
        phi23 = np.arctan(a3*np.sin(phi3)/(r3-a3*np.cos(phi3))) + np.arctan(r3*a3*np.sin(phi3)/(1-r3*a3*np.cos(phi3)))
        
        #phi23 = np.angle(Er23)
        
        phi12 = np.arctan((a2*abs(Er23)*np.sin(phi2+phi23))/(r2-a2*abs(Er23)*np.cos(phi2+phi23))) + np.arctan((r2*abs(Er23)*a2*np.sin(phi2+phi23))/(1-r2*abs(Er23)*a2*np.cos(phi2+phi23))) #phase of r12
        
        #phi12 = np.angle(Er12)
        
        phi1 = phi1 + 0.001    
        phi2 = phi2 + 0.001    
        phi3 = phi3 + 0.001
        
        Etai[i] = abs(Eta)**2 
        
        EtaImag[i] = Eta.imag
        EtaReal[i] = Eta.real
    
        PHIt[i] = np.angle(Eta)
        #PHIt[i] = np.arctan((a1*abs(Er12)*np.sin(phi1+phi12))/(r1-a1*abs(Er12)*np.cos(phi1+phi12))) + np.arctan((r1*a1*abs(Er12)*np.sin(phi1+phi12))/(1-r1*a1*abs(Er12)*np.cos(phi1+phi12))) #total eff phase
        
        rho = phi1 + phi12
        
        vg = 1/(n/c)*(((r1-a1*abs(Er12)*np.cos(rho))**2/((r1-a1*abs(Er12)*np.cos(rho))**2 + (a1*abs(Er12)*np.sin(rho))**2)) + ((1-r1*a1*abs(Er12)*np.cos(rho))**2 / ((1-r1*a1*abs(Er12)*np.cos(rho))**2 + (r1*a1*abs(Er12)*np.sin(rho))**2) ))
        
        #vgt[i] = vg
        
        
        #fp.write("{},{},{}\n".format(phi1,PHI[i], Etai[i]))        
        
        phi1t[i] = phi1
        phi2t[i] = phi12
        
        
dydx = np.diff(PHIt)/np.diff(phi1t)
vgt = (1/dydx)*c/n
ngt = dydx * n

phi1t = [phi1t[i]*radians for i in range(0,ran)]
phi2t = [phi2t[i]*radians for i in range(0,ran)]

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,5))

axs[0,0].set_title("Transmitted field")
axs[0,0].set_xlabel("Round Trip phase")
axs[0,0].set_ylabel("Transmittance")
axs[0,0].plot(phi1t,Etai, 'b', xunits=radians)
#plt.show()
#fig.savefig('Triple_resonator_trans.png', dpi=400)

#axs[0,1].set_xlim([-0.02,0.02])
axs[0,1].set_title("Phase")
axs[0,1].set_xlabel("Round Trip phase")
axs[0,1].set_ylabel("Total phase")
axs[0,1].plot(phi1t,PHIt, 'r', xunits=radians)
#plt.show()

#axs[1,0].set_xlim([-0.2,0.2])
axs[1,0].set_title("Transmission Imag vs Real")
axs[1,0].set_xlabel("Real axis")
axs[1,0].set_ylabel("Imaginart axis")
axs[1,0].plot(EtaReal,EtaImag, 'y')

#axs[1,1].set_xlim([-0.5,0.5])
axs[1,1].set_title("Group Velocity")
axs[1,1].set_xlabel("detuning")
axs[1,1].set_ylabel("Vg")
axs[1,1].plot(phi1t,np.append([1],vgt), 'g')
#axs[1,1].plot(phi1t,vgt, 'g')

fig.tight_layout()
plt.show()

#fig.savefig('Triple_resonator_phase.trans%d.png'%g, dpi=400)
