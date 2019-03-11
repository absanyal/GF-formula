import basisgeneration as bg
import progbar as pb
import numpy as np
import matplotlib.pyplot as plt
import time
import os
os.system('clear')

N = 8
n = int(N/2)
S_z = 0.5 * n

U = 0

Unn = 8

eta = 0.1

t = 1
tprime = t
I = complex(0, 1)


def getunn(s1):
    unn = 0
    c = s1.upconfig[:]
    for i in range(len(c)-1):
        if (c[i] == 1 and c[i+1] == 1):
            unn += 1
    return unn


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
                H[bi][bj] += getunn(state1)
    
    H = np.array(H)
    return H

def z(omega):
    return omega + I * eta


def G(omega, H):
    return np.linalg.inv(z(omega)
                         * np.eye(np.shape(H)[0], dtype=np.complex) - H)

# Generate required sectors
sectorlist = []

basis = bg.createbasis(N, n, S_z)

for l_n in range(0, n+1):
    l_Sz = 0.5 * l_n
    sector = bg.createsubbasis(basis, l_n, l_Sz)
    sectorlist.append(sector)

# for sector in sectorlist:
#     for s in sector:
#         print(s.getstate(), '\t', getunn(s))
#     print('*'*50)

b1 = sectorlist[1]

H_list = []

for sector in sectorlist:
    H = overlap(sector, sector)
    H_list.append(H)


