# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 15:46:02 2017

@author: AB Sanyal
"""

import numpy as np
import matplotlib.pyplot as plt

N = 3 #No of sites to be evaluated
epsilon00 = 0

n = N - 1 #size of the submatrix
eta = 0.1
epsilon = np.zeros(n)

def z(omega):
    return complex(omega, eta)

t = 1.0

omega_list = np.linspace(-8, 8, num = 1000)

A_list = []

#Populate the submatrix
for omega in omega_list:
    G_nn = np.zeros((n, n), dtype = np.complex)
    for i in range(n):
        G_nn[i][i] = z(omega) - epsilon[i]
    for i in range(0, n - 1):
        G_nn[i][i + 1] = t
    for i in range(1, n):
        G_nn[i][i-1] = t

    tau = np.zeros(n)
    tau[0] = t
    tau[n-1] = t

    G_nn = np.array(G_nn)
    
    G_nn_i = np.linalg.inv(G_nn)
    
    G_nn_i_tau = np.dot(G_nn_i, tau)
    tGt = np.dot(tau, G_nn_i_tau)
    G00 = pow( (z(omega) - epsilon00) - tGt, -1 )
    
    A_list.append( - (1/np.pi) * np.imag( G00 ) )



plt.plot(omega_list, A_list)
plt.title(str(n+1) + ' sites_mixed')
plt.savefig(str(n+1)+'sites_mixed.pdf')
plt.show()