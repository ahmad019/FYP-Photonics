# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 13:39:29 2018

@author: ahmad
"""


import numpy as np
import matplotlib.pyplot as plt
from basic_units import radians, degrees, cos

ran = 10000
phi1 = -0.5
phi2 = -0.5
r1 = 0.8998
r2 = 0.9998
lamd = 0.000001550
Q1 = 1e5
Q2 = 1e6
c = 299792458
b = 0.000025
l = 2*np.pi*b
n = 3.45

g1 = 0
g2 = 60
tau = (n*2*np.pi*b)/c

aa = (2*np.pi*n)/(Q1*lamd)
aa2 = (2*np.pi*n)/(Q2*lamd)
a1 = np.exp(-(aa-g1)*l/2)
a2 = np.exp(-(aa2-g2)*l/2)

#a1 = 0.88
#a2 = 0.9999 boyd

#r2 = a2
#r1 = r2*a1 #critical

#a1 = a1.round(decimals=2)
#a2 = a2.round(decimals=2)

Etai = np.ndarray(ran, float)
EtaImag = np.ndarray(ran, float)
EtaReal = np.ndarray(ran, float)
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
        
        #phi12 =  np.arctan((r2*a2*np.sin(phi2))/(1-r2*a2*np.cos(phi2))) - np.arctan((a2*np.sin(phi2))/(r2-a2*np.cos(phi2)))
        phi12 = np.arctan((a2*((r2**2) - 1)*np.sin(phi2))/(((a2**2) + 1)*r2 - a2*(r2**2 + 1)*np.cos(phi2)))
        
        #PHIt[i] = -np.arctan((a1*abs(Er12)*np.sin(phi12+phi1))/(r1-a1*abs(Er12)*np.cos(phi12+phi1))) + np.arctan((r1*a1*abs(Er12)*np.sin(phi12+phi1))/(1-r1*a1*abs(Er12)*np.cos(phi12+phi1)))
        #PHIt[i] = phi1 - (np.arctan((a2*np.sin(phi2))/(r2-a2*np.cos(phi2))) - np.arctan((a2*r2*np.sin(phi2))/(1-a2*r2*np.cos(phi2))))
        PHIt[i] = np.arctan((a1*abs(Er12)*((r1**2) - 1)*np.sin(phi1+phi12))/((a1**2)*abs(Er12)**2*r1-a1*abs(Er12)*((r1**2)+1)*np.cos(phi1+phi12)+r1))
        
        PHI12[i] = phi12
        
        phi1 = phi1 + 0.0001
        phi2 = phi2 + 0.0001
        
        Etai[i] = abs(Eta)**2 
        
        EtaImag[i] = Eta.imag
        EtaReal[i] = Eta.real
        
        dPHIt[i] = ((1-r1*abs(Er12)*a1*np.cos(phi12+phi1))**2/((1-r1*abs(Er12)*a1*np.cos(phi12+phi1))**2 + (r1*abs(Er12)*a1*np.sin(phi1+phi12))**2))  - ((r1-a1*abs(Er12)*np.cos(phi12+phi1))**2 / ((r1-a1*abs(Er12)*np.cos(phi12+phi1))**2  + (abs(Er12)*a1*np.sin(phi12+phi1))**2))
        
        fp.write("{},{},{}\n".format(phi1,PHIt[i], Etai[i]))
    
        phi1t[i] = phi1  #*2*np.pi/tau #(frequency detuning)
        phi2t[i] = phi2
      

over = r1 - phi2t[499]*a1
ove = r1 - r2*a1
ovr = r2 - a2

#dPHIr = ((1-r1*np.sin(phi))**2/((1-r1*np.sin(phi))**2 + (np.sin(phi))**2))  + ((r1-np.cos(phi))**2 / ((r1-np.cos(phi))**2  + (np.sin(phi))**2)) #singapore paper phase derivt

phi1t = phi1t.round(decimals=4)

dydx = np.diff(PHIt)/np.diff(phi1t)
#dydx = np.append(dydx, dydx[99])

vgt = (1/dydx)*c/n
ngt = dydx * n


#phi1t = [phi1t[i]*radians for i in range(0,ran)]
#phi2t = [phi2t[i]*radians for i in range(0,ran)]

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,5))
    
axs[0,0].set_title("Transmitted field")
#axs[0,0].set_xlim([-0.25,0.25])
#axs[0,0].set_ylim([0,1.2])
axs[0,0].set_xlabel("Frequency detuning (THz)")
axs[0,0].set_ylabel("Transmittance")
#axs[1,0].text(0.001,0.0000,"Gain1 =%f" %g1 + "\nGain2 =%f" %g2 + "\nUnder coupled\nr1=%.8f" %r1 + "\nr2=%.8f"%r2,fontsize=8)
axs[0,0].plot(phi1t,Etai, 'b')

#axs[0,1].set_xlim([-np.pi/8,np.pi/8])
axs[0,1].set_title("Phasor plot")
axs[0,1].set_xlabel("Frequency detuning (THz)")
axs[0,1].set_ylabel("Effective phase")
#axs[0,1].plot(Etai,np.append([1],ngt), 'm')
axs[0,1].plot(EtaReal,EtaImag, 'm')

#axs[1,0].set_xlim([-0.1,0.1])
#axs[1,0].set_ylim([-0.1,0.1])
axs[1,0].set_title("Effective Phase")
axs[1,0].set_xlabel("Frequency detuning (THz)")
axs[1,0].set_ylabel("Effective phase")
axs[1,0].plot(phi1t,PHIt, 'r')

#axs[1,1].set_xlim([-0.05,0.05])
#axs[1,1].set_ylim([-2000,2500])
#axs[1,1].set_title(" (n * dΦ/dφ)")
axs[1,1].set_xlabel("Frequency detuning (THz)")
axs[1,1].set_ylabel("Group Index")
#axs[1,1].plot(phi1t, n*(np.append([1],dydx)) , 'g')
axs[1,1].plot(phi1t, np.append([1],ngt) , color='g')

fig.tight_layout()
plt.show()

fig.savefig("coupled_ring_EIA.png",dpi=400)


print(over,ove,ovr)