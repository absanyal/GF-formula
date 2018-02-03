# -*- coding: utf-8 -*-
"""
Created on Sat Feb  3 17:10:25 2018

@author: AB Sanyal
"""

import numpy as np

def makeint(binnum):
    return int(binnum[::-1], 2)

def makebin(decnum):
    return (bin(decnum)[2:])[::-1]


class state:

    def __init__(self, upconfig, downconfig, phase = 0):
        self.upconfig = upconfig
        self.downconfig = downconfig
        self.N = len(upconfig)
        self.phase = phase

    def getstate(self):
        p = ''
        for i in range(self.N):
            if (self.upconfig[i] == 1 and self.downconfig[i] == 0):
                p += chr(8593)
            if (self.upconfig[i] == 0 and self.downconfig[i] == 1):
                p += chr(8595)
            if (self.upconfig[i] == 1 and self.downconfig[i] == 1):
                p += chr(8645)
            if (self.upconfig[i] == 0 and self.downconfig[i] == 0):
                p += '_'
        return p

    def binequiv(self):
        binnum = ''
        for bit in self.upconfig:
            binnum += str(bit)
        for bit in self.downconfig:
            binnum += str(bit)
        return binnum

    def intequiv(self):
        return makeint(self.binequiv())


def makestatefromint(N, intrep):
    binrep = makebin(intrep)
    while (len(binrep) < 2*N ):
        binrep += '0'

    upsector = [int(binrep[i]) for i in range(N)]
    downsector = [int(binrep[i]) for i in range(N,(2*N))]

    return state(upsector, downsector)

def makestatefrombin(N, binrep):
    binrep = str(binrep)
    upsector = [int(binrep[i]) for i in range(N)]
    downsector = [int(binrep[i]) for i in range(N,(2*N))]
    return state(upsector, downsector)


a = state([1, 1, 1, 1], [1, 1, 0, 0])
print(a.getstate())
print(a.binequiv())
print(a.intequiv())
b = makestatefromint(4, a.intequiv())
print(b.getstate())