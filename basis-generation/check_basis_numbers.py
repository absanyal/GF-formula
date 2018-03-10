# -*- coding: utf-8 -*-
"""
Created on Sat Mar 10 19:45:19 2018

@author: amit
"""

import basisgeneration as bg
import progbar as pb
import numpy as np
import matplotlib.pyplot as plt

N = 11
n = 2
#spin = (0.5 * n) % 1
spin = 1

us = [0 for i in range(N)]
ds = [0 for i in range(N)]
us[0] = 1
us[-1] = 0
ds[0] = 0
ds[-1] = 0
ws = bg.state(us, ds)

wsn = 0

basis = bg.createbasis(N, n, spin)

i = 0
for state in basis:
    print(i, state.getstate())
    if (state.getstate() == ws.getstate()):
        wsn = i
    i += 1

print("Required state", '\n', ws.getstate(), '\n',  "is at", wsn)

#print("*" * 80)

#sb = bg.createsubbasis(basis, 1, 0.5)
#
#i = 0
#for state in sb:
#    print(i, state.getstate())
#    i += 1