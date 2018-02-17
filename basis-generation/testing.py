# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 16:47:09 2018

@author: amit
"""

import basisgeneration as bg
import progbar as pb
import numpy as np
import matplotlib.pyplot as plt

N = 4
n = 4

p = 1

U = 8

spin = (0.5 * n) % 1
#spin = 1

t = 1

startpoint = -20
stoppoint = 20

print("\033c", end = '')

print("Calculating for", N, "sites with", n, "fermions and Sz =", spin)

print("Generating the basis...this may take a moment...")

basis = bg.createbasis(N, n, spin )

I = complex(0, 1)
eta = 0.1

print("The basis has", len(basis), "states.")

#i = 0
#for s in basis:
#    print(i, s.getstate())
#    i += 1
#"""
H = []

print("Generating the Hamiltonian...")

for s1 in basis:
    for s2 in basis:

        ta = 0
        for sigma in [-1, +1]:
            for i in range(0, N-1):
                s2a = bg.clonestate(s2)
                s2a.move(i, i+1, sigma)
                ta += bg.innerproduct(s1, s2a)

        tb = 0
        for sigma in [-1, +1]:
            for i in range(0, N-1):
                s2b = bg.clonestate(s2)
                s2b.move(i+1, i, sigma)
                ta += bg.innerproduct(s1, s2b)

        term = t * (ta + tb)

        H.append(term)

        i = len(H)
        pb.progressbar(i, 0, pow(len(basis),2))

H = np.array(H, dtype = np.complex)

H = np.reshape(H, ( len(basis), len(basis) ) )

for i in range(len(basis)):
    a = basis[i]
    particles = np.array(a.upconfig) + np.array(a.downconfig)
    for nump in particles:
        if (nump == 2):
            H[i][i] += U

#print(H)

def z(omega):
    return omega + I * eta

def G(omega):
    return np.linalg.inv(z(omega) \
    * np.eye(len(basis), dtype = np.complex) - H)


w_list = np.linspace(startpoint, stoppoint, 5000)

#A_list = []
#for w in w_list:
#    A_list.append( -(1/np.pi) * np.imag( G(w)[p][p] ) )
#
#plt.plot(w_list, A_list)
#plt.show()

total_A = np.zeros(len(w_list))

print("\nThe eigenvalues of the Hamiltonian are:")

ev = np.linalg.eigvalsh(H)

for e in ev:
    print( round(e, 2), sep = '\t', end = ' ' )

print("\nGenerating density of states...")

i = 0
for slabel in range(len(basis)):
    Ap_list = []
    for w in w_list:
        Ap_list.append(\
        (1/len(basis)) * (-1/np.pi) * np.imag( G(w)[slabel][slabel] ) )
    Ap_list = np.array(Ap_list)
    total_A += Ap_list

    pb.progressbar(i, 0, len(basis) - 1)
    i += 1

plt.plot(w_list, total_A)
plt.title(
        "DOS for " + str(N) + " sites with " + str(n) \
            + " fermions, Sz = " + str(spin) + ", U = " + str(U)
        )

plt.savefig(str(N) + "_" + str(n) + "_" + str(int(spin*10)) + str(U) + ".pdf")
plt.show()

print('')