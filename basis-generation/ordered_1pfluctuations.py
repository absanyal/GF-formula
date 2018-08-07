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

if not os.path.isdir('ordmatplots'):
    os.makedirs('ordmatplots')

N = 4
n = 4
Sz = (0.5 * n) % 1

U = 0
tprime = 1
t = 1

w = 0
eta = 0.1


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

def z(w):
    return w + complex(0, 1) * eta


def G(H, w):
    return np.linalg.inv(z(w) * np.eye(len(H)) - H)


def tmm(a, b, c):
    return dot(a, dot(b, c))

#############################################################################
#############################################################################

# Generate the 2-2 l_Sz = 0 basis

l_n = 2
l_Sz = 0

b_22_00 = bg.createlbsbasis(N, n, Sz, l_n, l_Sz)

leftsectorset = set([])

for state in b_22_00:
    p = gethalfsector1(state)
    leftsectorset.add(p)

leftsectorsetold = sorted(leftsectorset)

leftsectorset = leftsectorsetold

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

ordb_22_00 = []

for ls in leftsectorset:
    tempb = []
    for s in b_22_00:
        if (gethalfsector1(s) == ls):
            # print(s.getstate())
            tempb.append(s)
    # print('-' * 20)
    ordb_22_00 += tempb

H_22_00 = overlap(ordb_22_00, ordb_22_00)

##############################################################################

l_n = 2
l_Sz = 1

b_22_10 = bg.createlbsbasis(N, n, Sz, l_n, l_Sz)

leftsectorset = set([])

for state in b_22_10:
    p = gethalfsector1(state)
    leftsectorset.add(p)

leftsectorsetold = sorted(leftsectorset)

leftsectorset = leftsectorsetold

i = 0
for s in leftsectorset:
    print(i, s)
    i += 1

print("*" * 20)

ordb_22_10 = []

for ls in leftsectorset:
    tempb = []
    for s in b_22_10:
        if (gethalfsector1(s) == ls):
            # print(s.getstate())
            tempb.append(s)
    # print('-' * 20)
    ordb_22_10 += tempb

H_22_10 = overlap(ordb_22_10, ordb_22_10)

##############################################################################

l_n = 2
l_Sz = -1

b_22_01 = bg.createlbsbasis(N, n, Sz, l_n, l_Sz)

leftsectorset = set([])

for state in b_22_01:
    p = gethalfsector1(state)
    leftsectorset.add(p)

leftsectorsetold = sorted(leftsectorset)

leftsectorset = leftsectorsetold

i = 0
for s in leftsectorset:
    print(i, s)
    i += 1

print("*" * 20)

ordb_22_01 = []

for ls in leftsectorset:
    tempb = []
    for s in b_22_01:
        if (gethalfsector1(s) == ls):
            # print(s.getstate())
            tempb.append(s)
    # print('-' * 20)
    ordb_22_01 += tempb

H_22_01 = overlap(ordb_22_01, ordb_22_01)

##############################################################################

l_n = 3

b_31 = bg.createlfsbasis(N, n, Sz, l_n)

leftsectorset = set([])

for state in b_31:
    p = gethalfsector1(state)
    leftsectorset.add(p)

leftsectorsetold = sorted(leftsectorset)

leftsectorset = leftsectorsetold

i = 0
for s in leftsectorset:
    print(i, s)
    i += 1

print("*" * 20)

ordb_31 = []

for ls in leftsectorset:
    tempb = []
    for s in b_31:
        if (gethalfsector1(s) == ls):
            # print(s.getstate())
            tempb.append(s)
    # print('-' * 20)
    ordb_31 += tempb

H_31 = overlap(ordb_31, ordb_31)

##############################################################################

l_n = 1

b_13 = bg.createlfsbasis(N, n, Sz, l_n)

leftsectorset = set([])

for state in b_13:
    p = gethalfsector1(state)
    leftsectorset.add(p)

leftsectorsetold = sorted(leftsectorset)

leftsectorset = leftsectorsetold

i = 0
for s in leftsectorset:
    print(i, s)
    i += 1

print("*" * 20)

ordb_13 = []

for ls in leftsectorset:
    tempb = []
    for s in b_31:
        if (gethalfsector1(s) == ls):
            # print(s.getstate())
            tempb.append(s)
    # print('-' * 20)
    ordb_13 += tempb

H_13 = overlap(ordb_13, ordb_13)

##############################################################################

# d = len(ordb_22_p)

# H1 = H_22_p[0:d, 0:d]
# H2 = H_22_p[d:2 * d, d:2 * d]
# tau12 = H_22_p[0:d, d:2 * d]
# tau21 = H_22_p[d:2 * d, 0:d]

# w_min = min(eig) - 2
# w_max = max(eig) + 2

# w_list = np.linspace(w_min * t, w_max * t, 2000)
# A = []
# # for w in w_list:
# #     invGF = G(H, w)
# #     A.append((-1 / np.pi) * np.imag(np.trace(invGF)))
# #     pb.progressbar(w, w_list[0], w_list[-1])

# # plt.plot(w_list, A)
# # plt.show()


# print('Calculating DOS...')
# for w in w_list:

#     g1 = G(H1, w)
#     g2 = G(H2, w)

#     lG11 = G(H1, w)
#     rG11 = inv(inv(g1) - tmm(tau12, g2, tau21))

#     lG22 = inv(inv(g2) - tmm(tau21, g1, tau12))
#     rG22 = G(H2, w)

#     fG11 = inv(inv(g1) - tmm(tau12, rG22, tau21))
#     fG22 = inv(inv(g2) - tmm(tau21, lG11, tau12))

#     fG12 = tmm(lG11, tau12, fG22)
#     fG21 = tmm(rG22, tau21, fG11)

#     fGF = np.block([[fG11, fG12], [fG21, fG22]])

#     A.append((-1 / np.pi) * np.imag(np.trace(fGF)))
#     pb.progressbar(w, w_list[0], w_list[-1])

# plt.plot(w_list, A)

# plt.title('Formula inversion, U = ' + str(U) + ', l_n = ' + str(l_n))
# # plt.savefig('ordmatplots/formula_inversion_U4.pdf')
# plt.show()
