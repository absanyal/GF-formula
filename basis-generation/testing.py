# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 16:47:09 2018

@author: amit
"""

import basisgeneration as bg
import numpy as np

basis = bg.createbasis(4, 4, 0)

#i = 0
#
#for b in basis:
#    print(b.getstate())

#a = basis[0]
#b = basis[1]
#
#print(a.getstate())
#tempa = bg.clonestate(a)
#print(tempa.getstate())
##a.destroy(0, -1)
#print(a.getstate())
#print(tempa.getstate())
#a.phase *= -1
#
#print(bg.innerproduct(a, tempa))

#
#print(bg.innerproduct(a, a))

#a = bg.state([1, 0, 1, 0], [0, 1, 1, 0])
##print(a.getoccupation(1, 1))
#print(a.getstate())
#a.move(0, 0, 1)
#print(a.getstate())
###print(a.binequiv())
###print(a.intequiv())
##
#a.destroy(0, 1)
#print(a.getstate())
#print(a.phase)
#a.create(0, 1)
#print(a.getstate())
#print(a.phase)
#a.destroy(2, 1)
#print(a.getstate())
#print(a.phase)
#a.create(2, 1)
#print(a.getstate())
#print(a.phase)

N = 4
n = 8

basis = bg.createbasis(N, n, 0)

for s in basis:
    print(s.getstate())

#H = []
#
#U = 1
#t = 0
#
#for s1 in basis:
#    for s2 in basis:
#
#        ta = 0
#        for sigma in [-1, +1]:
#            for i in range(0, N-1):
#                s2a = bg.clonestate(s2)
#                s2a.move(i, i+1, sigma)
#                ta += bg.innerproduct(s1, s2a)
#
#        tb = 0
#        for sigma in [-1, +1]:
#            for i in range(0, N-1):
#                s2b = bg.clonestate(s2)
#                s2b.move(i+1, i, sigma)
#                ta += bg.innerproduct(s1, s2b)
#
##        tc = 0
##        for sigma in [-1, +1]:
##            for i in range(0, N):
##                s2c = bg.clonestate(s2)
##                tc += bg.innerproduct(s1, s2c)
#
#        term = t * (ta + tb)
#
#        H.append(term)

#H = np.array(H)
#
#H = np.reshape(H, ( len(basis), len(basis) ) )
#
#for i in range(len(basis)):
#    a = basis[i]
#    particles = np.array(a.upconfig) + np.array(a.downconfig)
#    for n in particles:
#        if (n == 2):
#            H[i][i] += U
#
#print(H)
#
#ev = np.linalg.eigvalsh(H)
#
#for e in ev:
#    print( round(e, 2), end = ' ' )