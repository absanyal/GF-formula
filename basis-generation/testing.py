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

N = 6
n = 4
Sz = (0.5 * n) % 1

l_n = 1
l_Sz = (0.5 * l_n) % 1

b = bg.createlbsbasis(N, n, Sz, l_n, l_Sz)

i = 0
for s in b:

    print(i, s.getstate(), s.intequiv(), (s.getleftnum(),  \
    s.getleftSz()), sep = '\t')

    i += 1
