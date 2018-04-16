# Mon Apr 16 00:01:29 IST 2018

import os
import time

import matplotlib.pyplot as plt
import numpy as np

import basisgeneration as bg
import progbar as pb

import multiprocessing

os.system('cls')
os.system('clear')

N = 6
total_n = 6
Sz = (0.5 * total_n) % 1

startpoint = -4
stoppoint = 28

n = int(total_n / 2)
p = n + 1
m = n - 1

t = -1
tprime = t
U = 8

eta = 0.05

print(
    "Calculating for", N, "sites with", total_n, "particles with total spin", Sz,
    " for the left block", int(total_n / 2), "with one-particle fluctuations."
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


print("Generating the necessary bases...")
t1 = time.perf_counter()
largs = [n, p, m]
args = [(N, total_n, Sz, i) for i in largs]
pool = multiprocessing.Pool()
listofbasis = pool.starmap(bg.createlfsbasis, args)

basisn = listofbasis[0]
basisp = listofbasis[1]
basism = listofbasis[2]

t2 = time.perf_counter()
print("Bases generated in", t2 - t1, "seconds.")

print("Generating the necessary matrices...")
t1 = time.perf_counter()
args = [
    (basisn, basisn), (basisp, basisp), (basism, basism),
    (basisp, basisn), (basism, basisn)
]

pool = multiprocessing.Pool()
listofoverlaps = pool.starmap(overlap, args)

H_n_n = listofoverlaps[0]
H_p_p = listofoverlaps[1]
H_m_m = listofoverlaps[2]

tau_p_n = listofoverlaps[3]
tau_n_p = np.transpose(tau_p_n)

tau_m_n = listofoverlaps[4]
tau_n_m = np.transpose(tau_m_n)

t2 = time.perf_counter()
print("Matrices generated in", t2 - t1, "seconds.")

# c = 0
# for i in basisn:
#     print(c, i.getstate())
#     c += 1

# print(H_n_n)

print("Generating DOS for block", n)
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
