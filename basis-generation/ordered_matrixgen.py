# Fri Mar 30 23:52:24 IST 2018

import itertools
import os
import time

import matplotlib.pyplot as plt
import numpy as np

from numpy import dot, around
from numpy.linalg import inv

import basisgeneration as bg
import progbar as pb

os.system('cls')
os.system('clear')

N = 4
n = 4
Sz = (0.5 * n) % 1

U = 0
tprime = 1
t = 1

w = 0
eta = 0.1

l_n = 2
l_Sz = (0.5 * l_n) % 1

b = bg.createlbsbasis(N, n, Sz, l_n, l_Sz)

# i = 0
# for s in b:

#     print(i, s.getstate())

#     i += 1

# print("*" * 20)


def gethalfsector1(state):
    p = list(state.getstate())
    p = p[:6]
    p = ''.join(p)
    # p.strip(' ')
    return p


def gethalfsector2(state):
    p = list(state.getstate())
    p = p[6:]
    p = ''.join(p)
    # p.strip(' ')
    return p


leftsectorset = set([])

for state in b:
    p = gethalfsector1(state)
    leftsectorset.add(p)

leftsectorsetold = sorted(leftsectorset)

leftsectorset = [0, 0, 0, 0]

leftsectorset[0] = leftsectorsetold[3]
leftsectorset[1] = leftsectorsetold[0]
leftsectorset[2] = leftsectorsetold[1]
leftsectorset[3] = leftsectorsetold[2]

i = 0
for s in leftsectorset:
    print(i, s)
    i += 1

print("*" * 20)

ordb = []

for ls in leftsectorset:
    tempb = []
    for s in b:
        if (gethalfsector1(s) == ls):
            # print(s.getstate())
            tempb.append(s)
    # print('-' * 20)
    ordb += tempb

# print("*" * 20)
# i = 0
# for s in ordb:
#     print(i, s.getstate(), sep='\t')
#     i += 1
# print('*' * 20)


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
#                    print((i, j, sigma), state1.getstate(), state2.getstate(),
#                          s2.getstate(), termtemp, sep = '\t')
#        print("*" * 80)

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


H = overlap(ordb, ordb)


def hamilprint(H):
    # print('\nThe Hamiltonian matrix is:')
    for i in range(len(H)):
        for j in range(len(H)):
            if (H[i, j] >= 0):
                print(' ', H[i, j], end='\t', sep='')
            else:
                print(H[i, j], end='\t', sep='')
        print('')


np.savetxt('matrixfixed.txt', H, fmt='%4.1f')
# hamilprint(H)
# print('*'*20)


def z(w):
    return w + complex(0, 1) * eta


def G(H, w):
    return np.linalg.inv(z(w) * np.eye(len(H)) - H)


def tmm(a, b, c):
    return dot(a, dot(b, c))


invGF = G(H, w)
np.savetxt('fullyinverted.txt', invGF, fmt='%5.2f')

# hamilprint(around(invGF, 3))
# print('*'*20)


#############################################################################
H1 = H[0:8, 0:8]
H2 = H[8:16, 8:16]
tau12 = H[0:8, 8:16]
tau21 = H[8:16, 0:8]

g1 = G(H1, w)
g2 = G(H2, w)

lG11 = G(H1, w)
rG11 = inv(inv(g1) - tmm(tau12, g2, tau21))

lG22 = inv(inv(g2) - tmm(tau21, g1, tau12))
rG22 = G(H2, w)

fG11 = inv(inv(g1) - tmm(tau12, rG22, tau21))
fG22 = inv(inv(g2) - tmm(tau21, lG11, tau12))

fG12 = tmm(lG11, tau12, fG22)
fG21 = tmm(rG22, tau21, fG11)

fGF = np.block([[fG11, fG12], [fG21, fG22]])

np.savetxt('rgf.txt', fGF, fmt='%5.2f')

truth = around(invGF, 3) == around(fGF, 3)

# print(truth)

i = 0
for a in truth:
    for b in a:
        if (b == True):
            i += 1

print(i, 'matches out of', pow(len(invGF), 2))
