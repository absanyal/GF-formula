###### Sat Mar 31 03:35:13 IST 2018

import basisgeneration as bg
import progbar as pb
import numpy as np
import matplotlib.pyplot as plt
import time
import os
os.system('cls')
os.system('clear')

N = 4
n = 4
Sz = (0.5 * n) % 1

l_n1 = 2
l_n2 = 2

t = 1
tprime = 2
U = 4

basis1 = bg.createlfsbasis(N, n, Sz, l_n1)
basis2 = bg.createlfsbasis(N, n, Sz, l_n2)


def mel(state1, state2):

    #calculate the hopping to right
    term = 0
    for sigma in [-1, 1]:
        for i in range(N):
            for j in range(N):
                if (i != j and abs(i - j) == 1):
                    s2 = bg.clonestate(state2)
                    s2.move(i, j, sigma)
                    termtemp = bg.innerproduct(state1, s2)
                    #print(termtemp)
                    term += termtemp

    return term


H = np.zeros((len(basis1), len(basis2)), dtype=np.float)

H_counter = 0
for bi in range(len(basis1)):
    for bj in range(len(basis2)):

        state1 = basis1[bi]
        state2 = basis2[bj]

        if (state1.getleftnum() == state2.getleftnum()):
            H[bi][bj] = t * mel(state1, state2)
            #H[bj][bi] = H[bi][bj]
        if (state1.getleftnum() != state2.getleftnum()):
            H[bi][bj] = tprime * mel(state1, state2)
            #H[bj][bi] = H[bi][bj]

        if ( l_n1 == l_n2 and bi == bj):
            a = basis1[bi]
            particles = np.array(a.upconfig) + np.array(a.downconfig)
            for nump in particles:
                if (nump == 2):
                    H[bi][bj] += U

        H_counter += 1
        pb.progressbar(H_counter, 0, len(basis1) * len(basis2))

print('\nThe Hamiltonian matrix is:')
for i in range(len(basis1)):
    for j in range(len(basis2)):
        if (H[i, j] >= 0):
            print(' ', H[i, j], end=' ', sep='')
        else:
            print(H[i, j], end=' ', sep='')
    print('')