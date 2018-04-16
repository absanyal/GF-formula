# Sun Apr 1 18:19:47 IST 2018

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

reqblocknums = [int(n / 2)]
reqblocknums = set(reqblocknums)
reqblocknums = list(reqblocknums)
reqblocknums.sort(reverse=True)

U = 8

eta = 0.1

spin = (0.5 * n) % 1
#spin = 0.5

t = -1
tprime = t
I = complex(0, 1)

#startpoint = 7.9
#stoppoint = 8.2

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
# i = 0
# for s in basis:
#     print(i, s.getstate())
#     i += 1

# Select the indices with required block numbers
print("Selecting indices of the matching block...")
req_indices = []
j = 0
for blocknum in reqblocknums:
    for i in range(len(basis)):
        if (basis[i].getleftnum() == blocknum):
            req_indices.append(i)
        j += 1
        pb.progressbar(j, 0, len(basis) - 1)

print('')


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
                    term += termtemp

    return term


print("Generating the Hamiltonian...")

# Full Hamiltonian construction

t_H_start = time.perf_counter()

basis1 = basis[:]
basis2 = basis[:]

H = np.zeros((len(basis1), len(basis2)))

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

t_H_stop = time.perf_counter()

print("\nHamiltonian matrix generated in",
      round(t_H_stop - t_H_start, 5), 's.')

# print('The Hamiltonian matrix is:')
# for i in range(len(H)):
#     for j in range(len(H)):
#         if (H[i, j] >= 0):
#             print(' ', H[i, j], end = ' ', sep = '')
#         else:
#             print(H[i, j], end = ' ', sep = '')
#     print('')


def z(omega):
    return omega + I * eta


def G(omega):
    return np.linalg.inv(z(omega)
                         * np.eye(len(basis), dtype=np.complex) - H)


ev = np.linalg.eigvalsh(H)
startpoint = np.floor(min(ev)) - 2
stoppoint = np.ceil(max(ev)) + 2

#print("The eigenvalues of the Hamiltonian are:")
# for e in ev:
#    print( round(e, 2), sep = '\t', end = ' ' )
#
# print('')

w_list = np.linspace(startpoint, stoppoint, 2000)

blocktext = ''
for a in reqblocknums:
    blocktext += str(a) + ' '

t_lsw_start = time.perf_counter()
print("Generating density of states for block",
      blocktext)

A_list = []
i = 0
for w in w_list:
    bdos = 0
    for snum in req_indices:
        bdos += -(1 / np.pi) * np.imag(G(w)[snum][snum] / len(req_indices))
    A_list.append(bdos)
    pb.progressbar(i, 0, len(w_list) - 1)
    i += 1

t_lsw_stop = time.perf_counter()

print("\nDensity of states calculated in",
      round(t_lsw_stop - t_lsw_start, 5), 's.')
plt.xlim(startpoint, stoppoint)
plt.plot(w_list, A_list)
plt.title("Density of states for the block " +
          blocktext)
plt.show()
print(startpoint, stoppoint)
