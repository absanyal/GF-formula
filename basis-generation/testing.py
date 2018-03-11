# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 16:47:09 2018

@author: amit
"""

import basisgeneration as bg
import progbar as pb
import numpy as np
import matplotlib.pyplot as plt

import time

N = 6
n = 6
spin = (0.5 * n) % 1

reptimes = 1

t1 = time.perf_counter()
for i in range(reptimes):
    bg.createbasis(N, n, spin)
t2 = time.perf_counter()

print((t2 - t1)/reptimes)