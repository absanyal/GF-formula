# Mon Apr 16 00:01:29 IST 2018

import os
import time

import matplotlib.pyplot as plt
import numpy as np

import basisgeneration as bg
import progbar as pb

os.system('cls')
os.system('clear')

N = 6
n = 6
Sz = (0.5 * n) % 1

startpoint = -4
stoppoint = 28


t = -1
tprime = t
U = 8

eta = 0.05

print(
    "Calculating for", N, "sites with", n, "particles with total spin", Sz,
    " for the left block", int(n / 2), "with one-particle fluctuations."
)


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


print("Generating basis for block", int(n / 2))
basisn = bg.createlfsbasis(N, n, Sz, n / 2)

p = int(n / 2 + 1)
print("Generating basis for block", p)
basisp = bg.createlfsbasis(N, n, Sz, p)

m = int(n / 2 - 1)
print("Generating basis for block", m)
basism = bg.createlfsbasis(N, n, Sz, m)

print("Generating Hamiltonian for block", n)
H_n_n = overlap(basisn, basisn)

print("Generating Hamiltonian for block", p)
H_p_p = overlap(basisp, basisp)

print("Generating Hamiltonian for block", m)
H_m_m = overlap(basism, basism)

print("Generating overlap between blocks", p, "and", n)
tau_p_n = overlap(basisp, basisn)
tau_n_p = np.transpose(tau_p_n)

print("Generating overlap between blocks", m, "and", n)
tau_m_n = overlap(basism, basisn)
tau_n_m = np.transpose(tau_m_n)

# c = 0
# for i in basisn:
#     print(c, i.getstate())
#     c += 1

# print(H_n_n)

print("Generating DOS for block", int(n / 2))
omega_list = np.linspace(startpoint, stoppoint, 2000)
# omega_list = [0]
wc = 0
A_list = []
for omega in omega_list:
    z = omega + complex(0, 1 * eta)
    Gi_n = z * np.eye(len(H_n_n)) - H_n_n
    Gi_p = z * np.eye(len(H_p_p)) - H_p_p
    Gi_m = z * np.eye(len(H_m_m)) - H_m_m

    G = np.linalg.inv(
        Gi_n
        - np.dot(tau_n_m, np.dot(np.linalg.inv(Gi_m), tau_m_n))
        - np.dot(tau_n_p, np.dot(np.linalg.inv(Gi_p), tau_p_n))
    )

    A = (-1 / np.pi) * np.imag(np.trace(G)) / len(G)
    A_list.append(A)

    wc += 1
    pb.progressbar(wc, 0, len(omega_list))

print('')

plt.plot(omega_list, A_list)
plt.show()
