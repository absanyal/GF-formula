# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 15:46:02 2017

@author: AB Sanyal
"""

import numpy as np
import matplotlib.pyplot as plt

eta = 0.1
epsilon = [0, 0, 0, 0]
2
def z(omega):
    return complex(omega, eta)

t = 1.0

t01 = t
t02 = t
t10 = t
t20 = t

def Gc11(omega):
    return pow( ( z(omega) - epsilon[1] ) - t * pow( z(omega) - epsilon[2], -1 ) * t, -1)

def Gc12(omega):
    return -Gc11(omega) * t * pow( z(omega) - epsilon[2], -1 )

def Gc21(omega):
    return pow( z(omega) - epsilon[2], -1 ) * t * Gc11(omega)

def Gc22(omega):
    return pow( z(omega) - epsilon[2], -1 ) - Gc21(omega) * t * pow( z(omega) - epsilon[2], -1 )

def tGct(omega):
    return t01 * Gc11(omega) * t10 + t01 * Gc12(omega) * t20 + t02 * Gc21(omega) * t10 + t02 * Gc22(omega) * t20

def Gcc00(omega):
    return pow( (z(omega) - epsilon[0]) - tGct(omega), -1 )

#A_list = []
omega_list = np.linspace(-8, 8, num = 1000)
#
#for omega in omega_list:
#    A_list.append( -1/(np.pi) * np.imag(Gcc00(omega)) )
#
#plt.plot(omega_list, A_list, 'g')

A_mat = []
for omega in omega_list:
    h = [ [ z(omega) - epsilon[0], t, 0, t ],
          [ t, z(omega) - epsilon[1], t, 0 ],
          [ 0, t, z(omega) - epsilon[2], t ],
          [ t, 0, t, z(omega) - epsilon[3] ]
        ]
    H = np.array(h)
    Hi = np.linalg.inv(H)
    A_mat.append( -1/(np.pi) * np.imag(Hi[0][0]) )

plt.plot(omega_list, A_mat)
plt.title('4 site matrix inverted directly')
plt.savefig('4sitematrixdirect.pdf')
plt.show()