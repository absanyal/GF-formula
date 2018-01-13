# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 08:41:51 2018

@author: amit
"""

import numpy as np
import matplotlib.pyplot as plt

I = complex(0, 1)
eta = 0.1
t = 1

eL = 0
eR = 0

U = 3

e = [2 * eL, 2 * eR, eL + eR, eL + eR]

tau01 = t * np.ones((2, 2))

tau10 = t * np.ones((2,2))

def z(omega):
    return omega + I * eta

#HB is the entire block Hamiltonian

def G00(omega):
    h = np.zeros((2,2), dtype = np.complex)
    h[0][0] = z(omega) - e[0] + U
    h[1][1] = z(omega) - e[1] + U
    return h

def G11(omega):
    h = np.zeros((2,2), dtype = np.complex)
    h[0][0] = z(omega) - e[2]
    h[1][1] = z(omega) - e[3]
    return h

omega_list = np.linspace(-6, 6, 5000)
A_list = []

for omega in omega_list:
    Gc = np.linalg.inv( G00(omega) - \
    np.dot( tau01, np.dot( np.linalg.inv( G11(omega) ), tau10 ) ) )
    A_list.append( (-1/np.pi) * np.imag( Gc[0][0] ) )

plt.plot(omega_list, A_list)
plt.title('2 sites with 2 particles')
plt.savefig('2-sites_2-particles.pdf')
plt.show()


#A_list = []
#for omega in omega_list:
#    Gc = np.linalg.inv( \
#    np.block([[G00(omega), tau01],
#               [tau10, G11(omega)]]))
#    A_list.append( (-1/np.pi) * np.imag( Gc[0][0] ) )
#
#plt.plot(omega_list, A_list)
#plt.title('2 sites with 2 particles by inversion')
##plt.savefig('2-sites_2-particles.pdf')
#plt.show()

#x = np.diag([e[0], e[1]])
#y = np.diag([e[2], e[3]])
#
#Hamil = np.block([[x, tau01],[tau10, y]])
#
#evals = np.sort( np.linalg.eigvals( Hamil ) )
#
#for ev in evals:
#    print( round( ev, 3 ), end = '\t' )
#
#print('')
