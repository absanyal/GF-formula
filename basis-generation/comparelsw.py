###### Tue Mar 20 11:44:42 IST 2018

import basisgeneration as bg
import progbar as pb
import numpy as np
import matplotlib.pyplot as plt
import time
import os
os.system('cls')
os.system('clear')

N = 6
n = 6

p = 272

U = 6

eta = 0.05

spin = (0.5 * n) % 1
#spin = 0.5

t = -1

I = complex(0, 1)

# startpoint = -6.0
# stoppoint = 10.0

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

# #Print the basis states set
# i = 0
# for s in basis:
#     print(i, s.getstate())
#     i += 1

#Hamiltonian construction with interface hopping

print("Generating the Hamiltonian with interface hopping...")

tprime = t
t_H_start = time.perf_counter()

basis1 = basis[:]
basis2 = basis[:]

H = np.zeros((len(basis1), len(basis2)))

def mel(state1, state2):

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

t_H_stop = time.perf_counter()

print("Hamiltonian matrix with interface hopping generated in",
      round(t_H_stop - t_H_start, 5), 's.')


def z(omega):
    return omega + I * eta


def G(omega):
    return np.linalg.inv(z(omega)
                         * np.eye(len(basis), dtype=np.complex) - H)

p = int(p)

ev = np.linalg.eigvalsh(H)
startpoint = np.floor(min(ev)) - 2
stoppoint = np.ceil(max(ev)) + 2

# print("The eigenvalues of the Hamiltonian are:")
# for e in ev:
#    print( round(e, 2), sep = '\t', end = ' ' )

# print('')

w_list = np.linspace(startpoint, stoppoint, 2000)

t_lsw_start = time.perf_counter()
print("Generating local spectral weight function for state:\n", \
       basis[p].getstate() )

A_list = []
i = 0
for w in w_list:
   A_list.append( -(1/np.pi) * np.imag( G(w)[p][p] ) )
   pb.progressbar(i, 0, len(w_list) - 1)
   i += 1

t_lsw_stop = time.perf_counter()
print("\nLocal spectral weight calculated in", \
   round(t_lsw_stop - t_lsw_start, 5), 's.')
plt.xlim(startpoint, stoppoint)
plt.plot(w_list, A_list, color = 'green', label = 'With interface hopping')
plt.title( "Local spectral weight function for the state " +\
            str(p) )
#plt.show()

#########################################################################

#Hamiltonian construction without interface hopping

print("Generating the Hamiltonian without interface hopping...")
tprime = 0

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

t_H_stop = time.perf_counter()

print("Hamiltonian matrix with interface hopping generated in",
      round(t_H_stop - t_H_start, 5), 's.')


p = int(p)

# ev = np.linalg.eigvalsh(H)
# startpoint = np.floor(min(ev)) - 2
# stoppoint = np.ceil(max(ev)) + 2

# print("The eigenvalues of the Hamiltonian are:")
# for e in ev:
#    print( round(e, 2), sep = '\t', end = ' ' )

# print('')

w_list = np.linspace(startpoint, stoppoint, 2000)

t_lsw_start = time.perf_counter()
print("Generating local spectral weight function for state:\n",
      basis[p].getstate())

A_list = []
i = 0
for w in w_list:
   A_list.append(-(1 / np.pi) * np.imag(G(w)[p][p]))
   pb.progressbar(i, 0, len(w_list) - 1)
   i += 1

t_lsw_stop = time.perf_counter()
print("\nLocal spectral weight calculated in",
      round(t_lsw_stop - t_lsw_start, 5), 's.')
plt.xlim(startpoint, stoppoint)
plt.plot(w_list, A_list, color='red', label = 'Without interface hopping')
plt.title("Local spectral weight function for the state " +
          str(p))
plt.legend(loc = 'best', framealpha = 0.5, fancybox = True)
plt.show()
