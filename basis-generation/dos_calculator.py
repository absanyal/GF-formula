# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 23:01:20 2018


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

N = 6
n = 3

p = 10

U = 0

Unn = 8

eta = 0.1

# spin = (0.5 * n) % 1
spin = 0.5 * n

t = 1
tprime = t
I = complex(0, 1)

#startpoint = 7.9
#stoppoint = 8.2

#print("\033c", end = '')


def getunn(s1):
    unn = 0
    c = s1.upconfig[:]
    for i in range(len(c)-1):
        if (c[i] == 1 and c[i+1] == 1):
            unn += 1
    return unn


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

print("Generating the basis...", waitmsg, sep='')

t_basis_start = time.perf_counter()

uobasis = bg.createbasis(N, n, spin)

basis = []

i = 0
for n_l in range(n + 1)[::-1]:
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
print("Basis generated and arranged in",
      round((t_basis_stop - t_basis_start), 5), 's.')

# Print the basis states set
i = 0
for s in basis:
    print(i, s.getstate(), getunn(s))
    i += 1

unnlist = []
for s in basis:
    unnlist.append(getunn(s))

print("Generating the Hamiltonian...")

# Alternative Hamiltonian construction without squaring

t_H_start = time.perf_counter()

basis1 = basis[:]
basis2 = basis[:]

H = np.zeros((len(basis1), len(basis2)))


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


for bi in range(len(basis1)):
    for bj in range(len(basis2)):

        state1 = basis1[bi]
        state2 = basis2[bj]

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

    Hprog = len(basis1) * (bi + 1) + bj + 1
    pb.progressbar(Hprog, 0, len(basis1) * len(basis2))

for i in range(len(H)):
    H[i, i] += Unn * unnlist[i]

t_H_stop = time.perf_counter()

print("\nHamiltonian matrix generated in",
      round(t_H_stop - t_H_start, 5), 's.')

print('The Hamiltonian matrix is:')
for i in range(len(H)):
    for j in range(len(H)):
        if (H[i, j] >= 0):
            print(' ', H[i, j], end=' ', sep='')
        else:
            print(H[i, j], end=' ', sep='')
    print('')


def z(omega):
    return omega + I * eta


def G(omega):
    return np.linalg.inv(z(omega)
                         * np.eye(len(basis), dtype=np.complex) - H)


#p = input("Enter the state number to calculate the spectral weight for: ")
p = int(p)

ev = np.linalg.eigvalsh(H)
startpoint = np.floor(min(ev)) - 2
stoppoint = np.ceil(max(ev)) + 2

print("The eigenvalues of the Hamiltonian are:")
for e in ev:
    print(round(e, 2), sep='\t', end=' ')

print('')

w_list = np.linspace(startpoint, stoppoint, 5000)

#t_lsw_start = time.perf_counter()
# print("Generating local spectral weight function for state:\n", \
#        basis[p].getstate() )
#
#A_list = []
#i = 0
# for w in w_list:
#    A_list.append( -(1/np.pi) * np.imag( G(w)[p][p] ) )
#    pb.progressbar(i, 0, len(w_list) - 1)
#    i += 1
#
#t_lsw_stop = time.perf_counter()
# print("\nLocal spectral weight calculated in", \
#    round(t_lsw_stop - t_lsw_start, 5), 's.')
#plt.xlim(startpoint, stoppoint)
#plt.plot(w_list, A_list)
# plt.title( "Local spectral weight function for the state " +\
#             str(p) )
# plt.show()

print("Generating density of states...")

t_DOS_start = time.perf_counter()
i = 0
Ap_list = []
for w in w_list:
    Ap_list.append((1 / len(basis)) * (-1 / np.pi) * np.imag(np.trace(G(w))))

    pb.progressbar(i, 0, len(w_list) - 1)
    i += 1

t_DOS_stop = time.perf_counter()

print("\nDensity of states generated in",
      round(t_DOS_stop - t_DOS_start, 5), 's.')

plt.xlim(startpoint, stoppoint)
plt.plot(w_list, Ap_list)
plt.ylabel('Density of states')
plt.xlabel('Frequency')

plt.title(
    "DOS for " + str(N) + " sites with " + str(n)
    + " " + fs + ", Sz = " + str(spin) + r", $U_{NN} = $" + str(Unn)
)

plt.savefig(
    str(N) + "_" + str(n) + "_" + str(int(spin * 10)) +
    "_" + str(U) + ".pdf")

plt.show()

f = open('dos_unn.dat', 'w')
for i in range(len(w_list)):
    f.write(str(w_list[i]) + '\t' + str(Ap_list[i]) + '\n')
f.close()
