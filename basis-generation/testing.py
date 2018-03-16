# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 16:47:09 2018

@author: amit
"""

import basisgeneration as bg
import progbar as pb
import numpy as np
import matplotlib.pyplot as plt
import itertools

import time

N = 4
n = 3
Sz = (0.5 * n) % 1

nu = int(Sz + 0.5 * n)
nd = int(n - nu)

#print(nu, nd)

us = [1 for i in range(nu)] + [0 for i in range(N-nu)]
ds = [1 for i in range(nd)] + [0 for i in range(N-nd)]

us_list = set(list(itertools.permutations(us)))
ds_list = set(list(itertools.permutations(ds)))

uobasis = []

for us_i in us_list:
    for ds_i in ds_list:
        uobasis.append(bg.state(us_i, ds_i))

basis = []

i = 0
for n_l in range(n+1)[::-1]:
    j = 0
    spins = 0.5 * np.array(list(range(-n_l, n_l + 1)))
    for Sz_l in spins:
        tempbasis = bg.createsubbasis(uobasis, n_l, Sz_l)
        basis += tempbasis

        for state in tempbasis:
            print(i, j, (n_l, Sz_l), state.getstate(),\
            state.intequiv(), sep = '\t')
            i += 1
            j += 1
#        #print(n_l, '\t', Sz_l)
#
    print("*" * 80)