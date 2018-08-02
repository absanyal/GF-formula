# Mon Apr 9 07:39:18 IST 2018

import multiprocessing
import os
import time

import matplotlib.pyplot as plt

import basisgeneration as bg

os.system('cls')
os.system('clear')

N_list = []
len_list = []

for N in range(2, 17):
    n = N
    Sz = (0.5 * n) % 1
    l_n = int(n / 2)
    l_Sz = (0.5 * l_n) % 1

    b = bg.createlbsbasis(N, n, Sz, l_n, l_Sz)

    print(N, len(b))

    N_list.append(N)
    len_list.append(len(b))

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_yscale('log')
plt.plot(N_list, len_list)
for i, j in zip(N_list, len_list):
    ax.annotate(str(j), xy=(i, j + 0.5))
plt.grid()
# plt.show()
plt.savefig('middle.pdf')
