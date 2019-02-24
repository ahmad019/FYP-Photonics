# -*- coding: utf-8 -*-
"""
Created on Sun Feb 24 12:59:36 2019

@author: ahmad

Code for solving equations of different systems using sympy symbolic expressions
"""

import sympy as s
import numpy as np
import matplotlib.pyplot as plt


ran = 200
PHI = np.ndarray(ran, float)
phi12t = np.ndarray(ran, float)
Etai = np.ndarray(ran, float)
phit = np.ndarray(ran, float)

phi = -1
lamd = 0.000001550
Q = 6000000
b = 0.000025
l = 2*np.pi*b
n = 3.45
c = 299792458
r = 0.989
rr = 0.978
aat = (2*np.pi*n)/(Q*lamd)
g1 = 1.9001
aa2t = (2*np.pi*n)/(Q*lamd)
g2 = 2.91001
a1t = np.exp(((g1-aat)*l)/2)
a2t = np.exp(((g2-aa2t)*l)/2)

phi1 = s.symbols(r'\phi1')
phi2 = s.symbols(r'\phi2')
phi12 = s.symbols(r'\phi12')
r1 =  s.symbols('r1')
r2 = s.symbols('r2')
r12 =  s.symbols('r12')
lamd = s.symbols(r'\lambda')
Q = s.symbols('Q')
b = s.symbols('b')
l = s.symbols('l')
n = s.symbols('n')
aa = s.symbols(r'\alpha1')
aa2 = s.symbols(r'\alpha2')
a1 = s.symbols('a1')
a2 = s.symbols('a2')
s.init_printing(use_unicode=True)

Er12 = (r2-a2*s.exp(s.I*phi2))/(1-r2*a2*s.exp(s.I*phi2))

Eta = (r1-r12*a1*s.exp(s.I*phi1))/(1-r1*r12*a1*s.exp(s.I*phi1))

phi12 = s.arg(Er12)

phie = s.arg(Eta)

for i in range(0,ran):
    
    r12t = Er12.subs({a2:a2t,phi2:phi,r2:rr})
    
    r12tc = s.N(r12t)
    
    comp = s.simplify(Eta.subs({a1:a1t,a2:a2t,phi1:phi,phi2:phi,r2:rr,r1:r,r12:r12tc}))
    
    comp2 = s.N(comp)
    
    Etai[i] = abs(comp2)**2
    
    phi12t[i] = phi12.subs({a2:1,phi2:phi,r2:rr})
    
    phiec = s.simplify(phie)
    
    PHI[i] = s.N(phiec.subs({a1:a1t,a2:a2t,phi2:phi,a2:a2t,phi1:phi,r1:r,r2:rr,r12:r12tc,phi12:phi12t[i]}))
    
    phit[i] = phi
    
    phi = phi + 0.01
    
    
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8,4))

axs[0].set_title("Transmission")
axs[0].set_xlabel("Round Trip phase")
axs[0].set_ylabel("Intensity")
axs[0].plot(phit,Etai, 'b')


axs[1].set_title("Coupling r12")
axs[1].set_xlabel("detuning")
axs[1].set_ylabel("Effective phase")
axs[1].plot(phit,phi12t, 'r')

fig2, axs = plt.subplots(1)

axs.set_title("Total Phase")
axs.set_xlabel("detuning")
axs.set_ylabel("Effective phase")
axs.plot(phit,PHI, 'r')

plt.text(-0.75,1.0,"Gain = 0\nUnder coupled\nr1=%f" %r + "\nr2=%f"%rr,fontsize=12, withdash=True)
#plt.grid()
fig.tight_layout()
fig2.tight_layout()
plt.show()

#fig.savefig('Coupled_sympy(trans).png', dpi=400)
#fig2.savefig('Coupled_sympy(phase).png', dpi=400)

#phie = s.atan((a1*abs(r12)*s.sin(phi1+phi12))/(r1-a1*abs(r12)*s.cos(phi1+phi12))) + s.atan((r1*a1*abs(r12)*s.sin(phi1+phi12))/(1-r1*a1*abs(r12)*s.cos(phi1+phi12)))


#Er12_E = Er12.expand(complex=True)

#c = s.atan(s.im(Et)/s.re(Et))
    

#print(c)
