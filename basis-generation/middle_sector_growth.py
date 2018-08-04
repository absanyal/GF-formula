# Mon Apr 9 07:39:18 IST 2018

import multiprocessing
import os
import time

import matplotlib.pyplot as plt

import basisgeneration as bg

os.system('cls')
os.system('clear')

N_list = list(range(2, 6))
len_list = []


def createlfsbasiswm(N, n, Sz, l_n):
    print("Calculating basis size", N)
    return bg.createlfsbasis(N, n, Sz, l_n)

args = [(n, n, (0.5 * n) % 1, int(n / 2)) for n in N_list]
# print(args)

pool = multiprocessing.Pool()
listofbasis = pool.starmap(createlfsbasiswm, args)


for b in listofbasis:
    len_list.append(len(b))

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_yscale('log')
plt.plot(N_list, len_list)
for i, j in zip(N_list, len_list):
    ax.annotate(str(j), xy=(i, j + 0.5))
plt.grid()
plt.show()
# plt.savefig('middle.pdf')
