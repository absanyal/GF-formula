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
n = 6
Sz = (0.5 * n) % 1

l_n = 3

b = bg.createlfsbasis(N, n, Sz, l_n)

i = 0
for s in b:

    print(i, s.getstate(), s.intequiv(), (s.getleftnum(),  \
    s.getleftSz()), sep = '\t')

    i += 1
