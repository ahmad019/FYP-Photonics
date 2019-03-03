# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 15:00:19 2019

@author: ahmad
"""


import numpy as np
import matplotlib.pyplot as plt

ran = 20000
phi1 = -1
phi2 = -1
phi3 = -1
#r1 =  0.9695465476797
#r2 = 0.898795189867767
r3 = 0.8968896789798389
lamd = 0.000001550
Q = 6000000
b = 0.000025
l = 2*np.pi*b
n = 3.45
c = 299792458

tau = (n*2*np.pi*b)/c

aa = (2*np.pi*n)/(Q*lamd)
g1 = 1.90001
aa2 = (2*np.pi*n)/(Q*lamd)
g2 = 2.981001
aa3 = (2*np.pi*n)/(Q*lamd)
g3 = 4.01
a1 = np.exp(((g1-aa)*l)/2)
a2 = np.exp(((g2-aa2)*l)/2)
a3 = np.exp(((g3-aa3)*l)/2)

r2 = r3*a2
r1 = r2*a1 #critical


Etai = np.ndarray(ran, float)
Erai = np.ndarray(ran, float)
PHI = np.ndarray(ran, float)

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
        
        phi1 = phi1 + 0.0001    
        phi2 = phi2 + 0.0001    
        phi3 = phi3 + 0.0001
        
        Etai[i] = abs(Eta)**2 
        
    
        #PHI[i] = np.angle(Eta)
        
        PHI[i] = np.arctan((a1*abs(Er12)*np.sin(phi1+phi12))/(r1-a1*abs(Er12)*np.cos(phi1+phi12))) + np.arctan((r1*a1*abs(Er12)*np.sin(phi1+phi12))/(1-r1*a1*abs(Er12)*np.cos(phi1+phi12))) #total eff phase
        
        #fp.write("{},{},{}\n".format(phi1,PHI[i], Etai[i]))        
        
        phi1t[i] = phi1
        phi2t[i] = phi12
        
        rho = phi1 + phi12
        
    
vg = 1/(n/c)*(((r1-a1*abs(Er12)*np.cos(rho))**2/((r1-a1*abs(Er12)*np.cos(rho))**2 + (a1*abs(Er12)*np.sin(rho))**2)) + ((1-r1*a1*abs(Er12)*np.cos(rho))**2 / ((1-r1*a1*abs(Er12)*np.cos(rho))**2 + (r1*a1*abs(Er12)*np.sin(rho))**2) ))


fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8,4))

axs[0].set_title("Transmitted field // EIA")
axs[0].set_xlabel("Round Trip phase")
axs[0].set_ylabel("Transmittance")
axs[0].plot(phi1t,Etai, 'b')
#plt.show()
#fig.savefig('Triple_resonator_trans.png', dpi=400)

#axs[1].set_xlim([-0.2,0.2])
axs[1].set_title("Phase")
axs[1].set_xlabel("Round Trip phase")
axs[1].set_ylabel("Total phase")
axs[1].plot(phi1t,PHI, 'r')
#plt.show()

fig.tight_layout()
plt.show()

#fig.savefig('Triple_resonator_phase.trans.png', dpi=400)


#np.savetxt("fsoo.csv", PHI, delimiter=",")
#np.savetxt("fsoo.csv", Etai, delimiter=",")


#with open('fsoo.csv', 'w') as fp:
    #np.savetxt(PHI, Etai, '%f', ',')
 #   fp.write("{},{},{}\n".format(phi1t,PHI, Etai))
    