# Wed Sep 26 12:24:15 IST 2018

import os
import time

import matplotlib.pyplot as plt
import numpy as np

import basisgeneration as bg
import progbar as pb

from math import factorial

os.system('clear')

N = 4
n = 2
Sz = 1

U = 0
t = 1
tprime = t

b = bg.createbasis(N, n, Sz)

i = 0
for state in b:
    print(i, state.getstate())
    i += 1


def mel(state1, state2):

    # calculate the hopping to right
    term = 0
    for sigma in [-1, 1]:
        for i in range(N):
            for j in range(N):
                if (i != j and abs(i - j) == 1):
                    s2 = bg.clonestate(state2)
                    s2.move(i, j, sigma)
                    termtemp = bg.innerproduct(state1, s2)
                    # print(termtemp)
                    term += termtemp
    return term


def overlap(basis1, basis2):
    H = np.zeros((len(basis1), len(basis2)))

    for bi in range(len(basis1)):
        for bj in range(len(basis2)):

            state1 = basis1[bi]
            state2 = basis2[bj]

            if (state1.getleftnum() == state2.getleftnum()):
                H[bi][bj] = t * mel(state1, state2)
                # H[bj][bi] = H[bi][bj]
            if (state1.getleftnum() != state2.getleftnum()):
                H[bi][bj] = tprime * mel(state1, state2)
                # H[bj][bi] = H[bi][bj]

            if (basis1 == basis2 and bi == bj):
                a = basis1[bi]
                particles = np.array(a.upconfig) + np.array(a.downconfig)
                for nump in particles:
                    if (nump == 2):
                        H[bi][bj] += U

    H = np.array(H)
    return H


def hamilprint(H):
    # print('\nThe Hamiltonian matrix is:')
    for i in range(len(H)):
        for j in range(len(H)):
            if (H[i, j] >= 0):
                print(' ', H[i, j], end='\t', sep='')
            else:
                print(H[i, j], end='\t', sep='')
        print('')


H = overlap(b, b)
hamilprint(H)
