# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 12:25:35 2018

@author: amit
"""

import basisgeneration as bg
import progbar as pb
import numpy as np
import matplotlib.pyplot as plt
import time
import os
os.system('cls')
os.system('clear')

N = 4
n = 2

p = 10

U = 8

eta = 0.1

spin = (0.5 * n) % 1
#spin = 0.5

t = 0
tprime = 2

I = complex(0, 1)

#startpoint = 7.9
#stoppoint = 8.2

#print("\033c", end = '')

if (n == 1):
    fs = "fermion"
else:
    fs = "fermions"

print("Calculating for", N, "sites with", n, fs, "and Sz =", spin)

waitmessagelist = [
                "this may take a while...",
                "please wait while this finishes...",
                "this will be done in a moment...",
                "this may take some time...",
                "this will take some time..."
                ]

waitmsg = np.random.choice(waitmessagelist)

print("Generating the basis...", waitmsg, sep = '')

t_basis_start = time.perf_counter()

uobasis = bg.createbasis(N, n, spin)

basis = []

i = 0
for n_l in range(n+1)[::-1]:
    spins = 0.5 * np.array(list(range(-n_l, n_l + 1)))
    for Sz_l in spins:
        tempbasis = bg.createsubbasis(uobasis, n_l, Sz_l)
        basis += tempbasis

t_basis_stop = time.perf_counter()

if (len(basis) == 1):
    statesm = "state."
else:
    statesm = "states."

print("The basis has", len(basis), statesm)
print("Basis generated and arranged in", \
    round((t_basis_stop - t_basis_start), 5), 's.' )

#i = 0
#for s in basis:
#    print(i, s.getstate())
#    i += 1

def mel(state1, state2):

    #calculate the hopping to right
    term = 0
    for sigma in [-1, 1]:
        for i in range(N):
            for j in range(N):
                if (i != j and abs(i-j) == 1):
                    s2 = bg.clonestate(state2)
                    s2.move(i, j, sigma)
                    termtemp = bg.innerproduct(state1, s2)
                    #print(termtemp)
                    term += termtemp
#                    print((i, j, sigma), state1.getstate(), state2.getstate(),
#                          s2.getstate(), termtemp, sep = '\t')
#        print("*" * 80)

    return term

H = np.zeros( (len(basis), len(basis)), dtype = np.float )

a = basis[0]
b = basis[1]

mel(a, b)

for bi in range(len(basis)):
    for bj in range(len(basis)):

        state1 = basis[bi]
        state2 = basis[bj]

        if (bi <= bj):
            if (state1.getleftnum() == state2.getleftnum()):
                H[bi][bj] = t * mel(state1, state2)
                H[bj][bi] = H[bi][bj]
            if (state1.getleftnum() != state2.getleftnum()):
                H[bi][bj] = tprime * mel(state1, state2)
                H[bj][bi] = H[bi][bj]

        if (bi == bj):
            a = basis[bi]
            particles = np.array(a.upconfig) + np.array(a.downconfig)
            for nump in particles:
                if (nump == 2):
                    H[bi][bj] += U

        Hprog = len(basis) * (bi + 1) + bj + 1
        pb.progressbar(Hprog, 0, len(basis) * len(basis))

print('\nThe Hamiltonian matrix is:')
for i in range(len(H)):
    for j in range(len(H)):
        if (H[i, j] >= 0):
            print(' ', H[i, j], end = ' ', sep = '')
        else:
            print(H[i, j], end = ' ', sep = '')
    print('')