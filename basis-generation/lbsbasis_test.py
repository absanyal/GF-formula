###### Fri Mar 30 23:52:24 IST 2018

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

    print(i, s.getstate(), s.intequiv(), (s.getleftnum(),
                                          s.getleftSz()), sep='\t')

    i += 1