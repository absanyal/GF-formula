# Sun Aug 5 21:02:49 IST 2018

import os
import numpy as np
from math import factorial
import matplotlib.pyplot as plt

os.system('clear')

N_max = 20


def C(N, ln, pu):
    pd = ln - pu
    qu = N / 2 - pu
    qd = N / 2 - ln + pu
    return pow(factorial(N / 2), 4) / (
        factorial(pu) * factorial(pd) * factorial(qu) * factorial(qd) *
        factorial(N / 2 - pu) * factorial(N / 2 - pd) *
        factorial(N / 2 - qu) * factorial(N / 2 - qd)
    )


sa = []
N_range = range(2, N_max + 1, 2)
for N in N_range:
    size = 0
    ln = N / 2
    for pu in range(0, int(ln + 1)):
        # print(N)
        size += C(N, ln, pu)
    sa.append(np.log(size))

plt.plot(N_range, sa)

sa = []
N_range = range(2, N_max + 1, 2)
for N in N_range:
    size = 0
    ln = N / 2 - 1
    for pu in range(0, int(ln + 1)):
        # print(N)
        size += C(N, ln, pu)
    sa.append(np.log(size))

plt.plot(N_range, sa)

plt.show()
