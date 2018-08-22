# Wed Aug 22 14:26:23 IST 2018

import os
import time

import matplotlib.pyplot as plt
import numpy as np

import basisgeneration as bg
import progbar as pb

from math import factorial

os.system('cls')
os.system('clear')

seed = [[0], [1]]

N = 10
n = 5

oldb = seed[:]
newb = []

# print(oldb)

expectedbasissize = factorial(N) / (factorial(N - n) * factorial(n))
print("Basis size should be", expectedbasissize)

finalbasis = []

counter = 0
while (len(finalbasis) <= expectedbasissize and len(oldb) > 0):
    newb = []
    for state in oldb:
        # print("Looking at", state)
        no_of_p_assigned = sum(state)
        no_of_sites_assigned = len(state)
        no_of_p_unassigned = n - sum(state)
        no_of_sites_unassigned = N - len(state)

        if (no_of_p_assigned == n and no_of_sites_assigned == N):
            finalbasis.append(state)
            # todelete = state[:]
            # oldb.remove(state)
            # print("Adding to basis", state)
        elif (no_of_sites_unassigned < no_of_p_unassigned):
            # todelete = state[:]
            # oldb.remove(state)
            # print("Deleting", state)
            continue
        elif (no_of_p_unassigned == 0 and no_of_sites_unassigned > 0):
            state += [0] * int(no_of_sites_unassigned)
            newb.append(state)
            # print("Padding with zeros", state)
        else:
            temps0 = state + [0]
            temps1 = state + [1]
            newb.append(temps0)
            newb.append(temps1)
            # print("Branching", state)
        # print("Temp basis is now at\n", oldb)
        # print('-' * 50)
        # time.sleep(1)
        counter += 1

    oldb = newb[:]
    # print("#" * 50)
    # time.sleep(1)
    # print(oldb)

# print("Final length:", len(finalbasis))

# if (len(finalbasis) == expectedbasissize):
#     print("SIZE MATCH")
# else:
#     print("SIZE MISMATCH")

print("Total number of steps required:", counter, "to generate",
      expectedbasissize, "elements, efficiency = ",
      round(expectedbasissize / counter, 3))

# print("Final basis is\n", finalbasis)
