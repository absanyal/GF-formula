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

dm = 0

seed = [[0], [1]]

N = 16
n = 8

oldb = seed[:]
newb = []

# print(oldb)

expectedbasissize = factorial(N) / (factorial(N - n) * factorial(n))
print("Calculating for", N, "sites with", n, "particles")
print("Sector size should be", expectedbasissize)

finalbasis = []

t1 = time.perf_counter()

counter = 0
while (len(finalbasis) <= expectedbasissize and len(oldb) > 0):
    newb = []
    for state in oldb:
        if (dm == 1):
            print('-' * 50)
            print("Looking at", state)
            time.sleep(1)
        counter += 1
        no_of_p_assigned = sum(state)
        no_of_sites_assigned = len(state)
        no_of_p_unassigned = n - sum(state)
        no_of_sites_unassigned = N - len(state)

        if (no_of_p_assigned == n and no_of_sites_assigned == N):
            finalbasis.append(state)
            # todelete = state[:]
            # oldb.remove(state)
            if (dm == 1):
                print("Adding to basis", state)
        elif (no_of_sites_unassigned < no_of_p_unassigned):
            # todelete = state[:]
            # oldb.remove(state)
            if (dm == 1):
                print("Skipping", state)
            continue
        elif (no_of_p_unassigned == no_of_sites_unassigned):
            state += [1] * int(no_of_sites_unassigned)
            # newb.append(state)
            finalbasis.append(state)
            if (dm == 1):
                print("Padding with ones and adding to final", state)
        elif (no_of_p_unassigned == 0 and no_of_sites_unassigned > 0):
            state += [0] * int(no_of_sites_unassigned)
            # newb.append(state)
            finalbasis.append(state)
            if (dm == 1):
                print("Padding with zeros and adding to final", state)
        else:
            temps0 = state + [0]
            temps1 = state + [1]
            newb.append(temps0)
            newb.append(temps1)
            if (dm == 1):
                print("Branching", state)
                print("Created", temps0, "and", temps1)

    oldb = newb[:]
    # print("#" * 50)
    # print("Basis is now at", oldb)
    # print("#" * 50)
    # if (oldb != []):
    #     p = len(oldb[0])
    # pb.progressbar(p, 0, N)

t2 = time.perf_counter()

print("Finished calculations in", round(t2 - t1, 3), "seconds.")

print("Final length:", len(finalbasis))

if (len(finalbasis) == expectedbasissize):
    print("SIZE MATCH")
else:
    print("SIZE MISMATCH")

print("Total number of steps required:", counter, "to generate",
      expectedbasissize, "elements, overcounting factor = ",
      round(counter / expectedbasissize, 3))

fullbasis = []

i = 0
for s1 in finalbasis:
    for s2 in finalbasis:
        s = bg.state(s1, s2)
        print(i, s.getstate())
        i += 1
