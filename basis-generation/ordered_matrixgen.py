# Fri Mar 30 23:52:24 IST 2018

import itertools
import os
import time

import matplotlib.pyplot as plt
import numpy as np

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

l_n = 2
l_Sz = (0.5 * l_n) % 1

b = bg.createlbsbasis(N, n, Sz, l_n, l_Sz)

i = 0
for s in b:

    print(i, s.getstate())

    i += 1

print("*" * 20)


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

leftsectorset = sorted(leftsectorset)

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
            print(s.getstate())
            tempb.append(s)
    print('-' * 20)
    ordb += tempb

print("*" * 20)
i = 0
for s in ordb:
    print(i, s.getstate(), sep='\t')
    i += 1
print('*' * 20)


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


H = np.zeros((len(ordb), len(ordb)), dtype=np.float)

for bi in range(len(ordb)):
    for bj in range(len(ordb)):

        state1 = ordb[bi]
        state2 = ordb[bj]

        if (bi <= bj):
            if (state1.getleftnum() == state2.getleftnum()):
                H[bi][bj] = t * mel(state1, state2)
                H[bj][bi] = H[bi][bj]
            if (state1.getleftnum() != state2.getleftnum()):
                H[bi][bj] = tprime * mel(state1, state2)
                H[bj][bi] = H[bi][bj]

        if (bi == bj):
            a = ordb[bi]
            particles = np.array(a.upconfig) + np.array(a.downconfig)
            for nump in particles:
                if (nump == 2):
                    H[bi][bj] += U

        Hprog = len(ordb) * (bi + 1) + bj + 1
        pb.progressbar(Hprog, 0, len(ordb) * len(ordb))

print('\nThe Hamiltonian matrix is:')
for i in range(len(H)):
    for j in range(len(H)):
        if (H[i, j] >= 0):
            print(' ', H[i, j], end=' ', sep='')
        else:
            print(H[i, j], end=' ', sep='')
    print('')
