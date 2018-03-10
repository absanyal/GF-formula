# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 15:46:02 2017

@author: AB Sanyal
"""

import numpy as np
import matplotlib.pyplot as plt

N = 100 #No of sites to be evaluated

eta = 0.05
epsilon = np.zeros(N)
t = 1

#for i in range(N):
#    if ((i % 50) == 0 and i>0):
#        epsilon[i] = 10

epsilon[5] = 0.5

def z(omega):
    return complex(omega, eta)

def Gcpp(p, omega):
    if (p == N-1):
        return 1 / ( z(omega) - epsilon[N-1] )
    else:
        return 1 / ( z(omega) - epsilon[p - 1] - \
        t * Gcpp(p + 1, omega) * t )

omega_list = np.linspace(-4, 4, num = 1000)

A_list = []

for omega in omega_list:
    A_list.append( (-1/np.pi) * np.imag( Gcpp(0, omega) ) )




plt.plot(omega_list, A_list)
plt.title(str(N) + ' sites with formula')
#plt.savefig(str(N)+'sites_formula_full.pdf')
plt.show()
