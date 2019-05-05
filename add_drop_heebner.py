# -*- coding: utf-8 -*-
"""
Created on Sat Feb 23 18:08:21 2019

@author: ahmad

Add drop micro ring resonator equations as from heebner book
"""


import numpy as np
import matplotlib.pyplot as plt

ran = 200
phi = -1
#r1 = 0.9797
r2 = 0.899
c = 299792458
lamd = 0.000001550
Q = 6000000
b = 0.000025
l = 2*np.pi*b
n = 3.45

g = 0 #zero for passive
aa = (2*np.pi*n)/(Q*lamd)
a = np.exp(((g-aa)*l)/2)

#r2 = a
r1 = r2*a #critically coupled

t = np.sqrt(1-(r1**2))
t2 = np.sqrt(1-(r2**2))


Etai = np.ndarray(ran, float)
Erai = np.ndarray(ran, float)
PHIt = np.ndarray(ran, float)
PHIr = np.ndarray(ran, float)
phit = np.ndarray(ran, float)
dPHIr = np.ndarray(ran, float)

for i in range(0,ran):
    
    Era = ((r2*a)**2 - 2*r1*r2*a*np.cos(phi) + r1**2) / (1 - 2*r1*r2*a*np.cos(phi) + (r1*r2*a)**2)
    
    Eta = ((1-r1**2)*(1-r2**2)*a)/(1 - 2*r1*r2*a*np.cos(phi) + (r1*r2*a)**2)

    PHIt[i] = phi/2    
    PHIr[i] = np.arctan((r2*np.sin(phi))/(r1-r2*np.cos(phi))) + np.arctan((r1*r2*np.sin(phi))/(1-r1*r2*np.cos(phi)))
    
    phi = phi + 0.01    
    
    Etai[i] = abs(Eta)**2
    Erai[i] = abs(Era)**2
    
    dPHIr[i] = ((r1-r2*np.cos(phi))**2/((r1-r2*np.cos(phi))**2 + (r2*np.sin(phi))**2))  + ((1-r1*r2*np.cos(phi))**2 / ((1-r1*r2*np.cos(phi))**2  + (r1-r2*np.sin(phi))**2))

    
    phit[i] = phi

dPHIt = 1/2
#dPHIr = ((r1-r2*np.cos(phi))**2/((r1-r2*np.cos(phi))**2 + (r2*np.sin(phi))**2))  + ((1-r1*r2*np.cos(phi))**2 / ((1-r1*r2*np.cos(phi))**2  + (r1-r2*np.sin(phi))**2))

vgt = (2*c)/n
ngt = n/2

vgr = (1/dPHIr)*c/n
ngr = dPHIr * n

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8,4))

axs[0].set_title("Reflection + Transmission")
axs[0].set_xlabel("Round Trip phase")
axs[0].set_ylabel("Intensity")
axs[0].plot(phit,Etai, 'b')
axs[0].plot(phit,Erai, 'r')

#axs[1].set_xlim([-0.2,0.2])
axs[1].set_title("Transmission Phase")
axs[1].set_xlabel("detuning")
axs[1].set_ylabel("Effective phase")
axs[1].plot(phit,PHIt, 'r')

fig2, axs = plt.subplots(nrows=1, ncols=2, figsize=(8,4))

axs[0].set_title("Reflection Phase")
axs[0].set_xlabel("detuning")
axs[0].set_ylabel("Effective phase")
#axs[0].set_xlim([-0.01,0.01])
#axs.set_ylim([-1,1])
#axs.legend("gain = 1.001")#,loc="upper right")
axs[0].plot(phit,PHIr, 'r')

axs[1].set_title("Group Velocity")
axs[1].set_xlabel("detuning")
axs[1].set_ylabel("Velocity")
axs[1].plot(phit,vgr, 'r')

#plt.text(0.5, 2,"Gain =%.2f" %g + "\nCritical coupled\nr1=%.2f" %r1 + "\nr2=%.2f"%r2,fontsize=12,transform=axs[1].transAxes, withdash=True)
plt.grid()
fig.tight_layout()
fig2.tight_layout()
plt.show()

#fig2.savefig('add_drop_phase(Heebner)_over.png', dpi=400)