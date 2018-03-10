# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 15:46:02 2017

@author: AB Sanyal
"""

import numpy as np
import matplotlib.pyplot as plt

n = 100 #number of sites
eta = 0.05
epsilon = np.zeros(n)
epsilon[50] = 10

def z(omega):
    return complex(omega, eta)

t = 1.0

omega_list = np.linspace(-4, 4, num = 1000)

A_mat = []
for omega in omega_list:
    H = np.zeros((n, n), dtype = np.complex)
    for i in range(n):
        H[i][i] = z(omega) - epsilon[i]
    for i in range(0, n - 1):
        H[i][i + 1] = t
    for i in range(1, n):
        H[i][i-1] = t
#    H[0][n - 1] = 0
#    H[n - 1][0] = 0
    Hi = np.linalg.inv(H)
    A_mat.append( -1/(np.pi) * np.imag(Hi[0][0]) )
    #A_mat.append( -1/(np.pi) * np.imag(np.trace(Hi)) / n )

plt.plot(omega_list, A_mat)
plt.title(str(n) + ' site matrix inverted directly with no PBC')
#plt.savefig(str(n)+'sitematrixinvert-no-pbc.pdf')
#plt.savefig(str(n)+'sitematrixinvert-with-pbc.pdf')
plt.show()
