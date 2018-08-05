# Mon Apr 9 07:39:18 IST 2018

import multiprocessing
import os
import time

import matplotlib.pyplot as plt

import basisgeneration as bg

os.system('cls')
os.system('clear')

N_list = list(range(2, 13))
len_list = []


def createlfsbasiswm(N, n, Sz, l_n):
    p = bg.createlfsbasis(N, n, Sz, l_n)
    print("Finished calculating basis size", N)
    return p


args = [(n, n, (0.5 * n) % 1, int(n / 2 - 1)) for n in N_list]
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
plt.title("n-1 n+1 splitting")
# plt.show()
plt.savefig('middle_nminus1.pdf')

f = open("middle_nminus1.dat", 'w')

for i in range(len(N_list)):
    f.write(str(N_list[i]) + '\t' + str(len_list[i]) + '\n')

f.close()
