# -*- coding: utf-8 -*-
"""
Created on Sat Mar 10 22:49:48 2018

@author: amit
"""

import basisgeneration as bg
import progbar as pb
import numpy as np
import matplotlib.pyplot as plt

N = 4
n = 4
spin = (0.5 * n) % 1

uobasis = bg.createbasis(N, n, spin)

#i = 0
#for state in uobasis:
#    print(i, '\t', state.getstate())
#    i += 1
#
#print("*" * 80)

basis = []

i = 0
for n_l in range(n+1):
    spins = 0.5 * np.array(list(range(-n_l, n_l + 1)))
    for Sz_l in spins:
        tempbasis = bg.createsubbasis(uobasis, n_l, Sz_l)
        basis += tempbasis

        for state in tempbasis:
            print(i, '\t', state.getstate())
            i += 1
#        #print(n_l, '\t', Sz_l)
#
print("*" * 80)

i= 0
for state in basis:
    print(i, '\t', state.getstate())
    i += 1