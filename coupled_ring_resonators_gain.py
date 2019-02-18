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
r1 =  0.943565476797
r2 = 0.998189867767
lamd = 0.000001550
Q = 6000000
b = 0.000025
l = 2*np.pi*b
n = 3.45

aa = (2*np.pi*n)/(Q*lamd)
g1 = 0.9001
aa2 = (2*np.pi*n)/(Q*lamd)
g2 = 0.91001
aa3 = (2*np.pi*n)/(Q*lamd)
g3 = 0.9101
a1 = np.exp(((g1-aa)*l)/2)
a2 = np.exp(((g2-aa2)*l)/2)
a3 = np.exp(((g3-aa3)*l)/2)

#r1 = r2*a1 #critical

Etai = np.ndarray(ran, float)
Erai = np.ndarray(ran, float)
PHI = np.ndarray(ran, float)

phi1t = np.ndarray(ran, float)
phi2t = np.ndarray(ran, float)

with open('Coupled_resonator_phi.phase.trans.csv', 'w') as fp:

    for i in range(0,ran):
        
        
        Er12 = (r2-a2*np.exp(1j*phi2))/(1-r2*a2*np.exp(1j*phi2)) #coupling r12
        
        Eta = (r1-Er12*a1*np.exp(1j*phi1))/(1-r1*Er12*a1*np.exp(1j*phi1)) #transmission
        
        #Era = (r1-r2*a1*np.exp(1j*phi1))/(1-r1*r2*a1*np.exp(1j*phi1))
        
        
        phi1 = phi1 + 0.001    
        phi2 = phi2 + 0.001    
        
        Etai[i] = abs(Eta)**2 
        #Erai[i] = abs(Era)**2
        #PHI[i] = np.arctan(Eta.imag/Eta.real)
        
        PHI[i] = np.arctan(r1*a1*abs(Er12)*np.sin(phi2+phi1))/(1-r1*a1*abs(Er12)*np.cos(phi2+phi1)) - np.arctan((a1*abs(Er12)*np.sin(phi2+phi1))/(r1-a1*abs(Er12)*np.cos(phi2+phi1)))
        
        fp.write("{},{},{}\n".format(phi1,PHI[i], Etai[i]))
        
        phi1t[i] = phi1
        phi2t[i] = phi2
        
    


fig, axs = plt.subplots(2)

axs[0].set_title("Transmitted field // EIA")
axs[0].plot(phi1t,Etai, 'b')
#plt.show()
#fig.savefig('Triple_resonator_trans.png', dpi=400)

axs[1].set_xlim([-0.15,0.15])
axs[1].set_title("Phase")
axs[1].plot(phi2t,PHI, 'r')
#plt.show()

fig.tight_layout()
plt.show()

fig.savefig('Coupled_resonator_phase.trans.png', dpi=400)

#plt.ylim(top=1.02,bottom=0.9)

#plt.title("Transmitted field // EIA")
#plt.plot(phi1t,Etai, 'b')
#plt.show()

#plt.xlim([-0.12,0.12])
#plt.title("Phase")
#plt.plot(phi2t,PHI, 'r')
#plt.show()

#fig.savefig('tempds.png', dpi=200)