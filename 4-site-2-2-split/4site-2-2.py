# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 19:22:20 2018

@author: amit
"""

import numpy as np
import matplotlib.pyplot as plt

I = complex(0, 1)
eta = 0.1
t = 1

e = [0, 0, 0, 0]

tau01 = np.zeros((2,2))
tau01[1][0] = t

tau10 = np.zeros((2,2))
tau10[0][1] = t

def z(omega):
    return omega + I * eta

def G1(omega):
    h = np.zeros((2,2), dtype = np.complex)
    h[0][0] = z(omega) - e[0]
    h[1][1] = z(omega) - e[1]
    h[1][0] = t
    h[0][1] = t
    return h

def G2(omega):
    h = np.zeros((2,2), dtype = np.complex)
    h[0][0] = z(omega) - e[2]
    h[1][1] = z(omega) - e[3]
    h[1][0] = t
    h[0][1] = t
    return h

omega_list = np.linspace(-6, 6, 5000)
A_list = []

for omega in omega_list:
    A = np.linalg.inv( G1(omega) - \
    np.dot( tau01, np.dot( np.linalg.inv(G2(omega)), tau10 ) ) )
    A_list.append( (-1 / np.pi) * np.imag( A[0][0] ) )

plt.plot(omega_list, A_list)
plt.title('4 sites with 2+2 partition')
#plt.savefig('4-sites_2p2.pdf')
plt.show()
