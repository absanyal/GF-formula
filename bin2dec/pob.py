import numpy as np
import math
from math import factorial
import matplotlib.pyplot as plt
#import scipy

def C(n,r):
    if (r>n or n<0 or r<0):
        return 0
    else:
        return int(factorial(n) / factorial(r) / factorial(n-r))

def binlist(num, sites):
    b = [int(x) for x in str(num)]
    while(len(b) < sites):
        b = [0] + b
    return b


def b(x, j):
    return x[len(x)-j-1]

def posval(x, j):
    pj = 0
    for i in range(j+1):
        pj += b(x, i)
    return int(b(x, j) * C(j, pj))

def V(x):
    v = 0
    for j in range(len(x)):
        v += posval(x, j)
    return int(v)

def dec(x):
    s = 0
    for j in range(len(x)):
        s += pow(2, j) * b(x, j)
    return s

def valtobin(n, r, val):
    binnum = [0 for i in range(n)]
    j = n
    temp = val
    for k in range(r, 0, -1):
        bj = 0
        while (bj != 1):
            j -= 1
            p = C(j, k)
            if (temp >= p):
                temp -= p
                bj = 1
        binnum[j] = bj
    binnum.reverse()
    return binnum



#x = [0, 0, 1, 1, 0, 1, 0, 1, 1, 0]
#n = len(x)
#r = sum(x)
#print(x)
#print(dec(x))
#
##for j in range(len(x)):
##    print(j, b(x,j), posval(x,j))
#
#v = V(x)
#
#print(v)
#
#print(valtobin(n, r, v))

n = 16
r = 8
combos = C(n, r)

for i in range(combos):
    v=valtobin(n, r, i)
    print(i, v, dec(v), sep='\t')
