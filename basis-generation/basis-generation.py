# -*- coding: utf-8 -*-
"""
Created on Sat Feb  3 17:10:25 2018

@author: AB Sanyal
"""

import matplotlib.pyplot as plt

#Reverse a binary and convert to decimal integer
def makeint(binnum):
    return int(binnum[::-1], 2)

#Convert to binary and reverse
def makebin(decnum):
    return (bin(decnum)[2:])[::-1]


class state:

    """
    A state is characterized by the configurations of the up and down
    sectors and the phase of the state. '1' = up spin, '-1' = down spin.
    """

    def __init__(self, upconfig, downconfig, phase = 1):
        self.upconfig = upconfig
        self.downconfig = downconfig
        self.N = len(upconfig)
        self.phase = phase

    #Generate the state in common electronic notationb
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
                p += ' _ '
        return p

    #Given a state, generate its unique binary representation
    def binequiv(self):
        binnum = ''
        for bit in self.upconfig:
            binnum += str(bit)
        for bit in self.downconfig:
            binnum += str(bit)
        return binnum

    #Given a state, generate its unique decimal integer representation
    def intequiv(self):
        return makeint(self.binequiv())

    #Returns the total number of particles in the state
    def getnumparticles(self):
        n = 0
        for i in self.upconfig:
            n += i
        for j in self.downconfig:
            n += j
        return n

    #Returns the total spin of the state
    def getSz(self):
        s = 0
        for i in self.upconfig:
            s += i * 1
        for j in self.downconfig:
            s += j * (-1)
        return s

    #Definition of the c, c_dagger and n operators

    #Creation operator
    def create(self, site, sigma):
        if (sigma == 1 and self.upconfig[site] == 0):

            n1 = 0
            for i in range(0, site):
                n1 += self.upconfig[i]

            n2 = 0
            for i in range(0, site):
                n2 += self.downconfig[i]

            self.phase *= pow(-1, n1 + n2)
            self.upconfig[site] = 1

        elif (sigma == -1 and self.downconfig[site] == 0):

            n1 = 0
            for i in range(0, site + 1):
                n1 += self.upconfig[i]

            n2 = 0
            for i in range(0, site + 1):
                n2 += self.downconfig[i]

            self.phase *= pow(-1, n1 + n2)
            self.downconfig[site] = 1

        elif (sigma == 1 and self.upconfig[site] == 1):
            self.upconfig = [0 for i in range(self.N)]
            self.downconfig = [0 for i in range(self.N)]
            self.phase = 0

        elif (sigma == -1 and self.downconfig[site] == 1):
            self.upconfig = [0 for i in range(self.N)]
            self.downconfig = [0 for i in range(self.N)]
            self.phase = 0

    #Destruction operator
    def destroy(self, site, sigma):
        if (sigma == 1 and self.upconfig[site] == 1):

            n1 = 0
            for i in range(0, site+1):
                n1 += self.upconfig[i]

            n2 = 0
            for i in range(0, site+1):
                n2 += self.downconfig[i]

            self.phase *= pow(-1, n1 + n2)
            self.upconfig[site] = 0

        elif (sigma == -1 and self.downconfig[site] == 1):

            n1 = 0
            for i in range(0, site):
                n1 += self.upconfig[i]

            n2 = 0
            for i in range(0, site + 1):
                n2 += self.downconfig[i]

            self.phase *= pow(-1, n1 + n2)
            self.downconfig[site] = 0

        elif (sigma == 1 and self.upconfig[site] == 0):
            self.upconfig = [0 for i in range(self.N)]
            self.downconfig = [0 for i in range(self.N)]
            self.phase = 0

        elif (sigma == -1 and self.downconfig[site] == 0):
            self.upconfig = [0 for i in range(self.N)]
            self.downconfig = [0 for i in range(self.N)]
            self.phase = 0

    #Occupation number operator for a given site and sigma
    def getoccupation(self, site, sigma):
        self.destroy(site, sigma)
        self.create(site, sigma)

###########################################################################

#Given an integer representation and the number of sites,
#returns a unique state.

def makestatefromint(N, intrep):
    binrep = makebin(intrep)
    while (len(binrep) < 2*N ):
        binrep += '0'

    upsector = [int(binrep[i]) for i in range(N)]
    downsector = [int(binrep[i]) for i in range(N,(2*N))]

    return state(upsector, downsector)

#Given a binary representation, returns a unique state.
def makestatefrombin(binrep):
    N = len(binrep)/2
    binrep = str(binrep)
    upsector = [int(binrep[i]) for i in range(N)]
    downsector = [int(binrep[i]) for i in range(N,(2*N))]
    return state(upsector, downsector)

#Given the number of sites, number of particles and total spin,
#returns the complete basis set.

def createbasis(N, n_particles, S_z = 0):
    basis = []

    #make max sector
    upsector = [0 for i in range(N)]
    downsector = [0 for i in range(N)]
    for i in range(1, n_particles+1):
        downsector[-i] = 1
    maxstate = state(upsector, downsector)
    #print(maxstate.getstate())
    maxint = maxstate.intequiv()
    #print(maxint)

    #make min sector
    upsector = [0 for i in range(N)]
    downsector = [0 for i in range(N)]
    for i in range(n_particles):
        upsector[i] = 1
    minstate = state(upsector, downsector)
    #print(minstate.getstate())
    minint = minstate.intequiv()
    #print(minint)

    #make all states and select the appropriate ones
    for i in range(minint, maxint + 1):
        newbasis = makestatefromint(N, i)
        if (newbasis.getSz() == S_z and \
        newbasis.getnumparticles() == n_particles):
            basis.append(newbasis)

    return basis

#Given the number of sites, number of particles and total spin,
#returns the number of elements basis set.

def getbasissize(N, n_particles, S_z = 0):
    basis_size = 0

    #make max sector
    upsector = [0 for i in range(N)]
    downsector = [0 for i in range(N)]
    for i in range(1, n_particles+1):
        downsector[-i] = 1
    maxstate = state(upsector, downsector)
    #print(maxstate.getstate())
    maxint = maxstate.intequiv()
    #print(maxint)

    #make min sector
    upsector = [0 for i in range(N)]
    downsector = [0 for i in range(N)]
    for i in range(n_particles):
        upsector[i] = 1
    minstate = state(upsector, downsector)
    #print(minstate.getstate())
    minint = minstate.intequiv()
    #print(minint)

    #make all states and select the appropriate ones
    for i in range(minint, maxint + 1):
        newbasis = makestatefromint(N, i)
        if (newbasis.getSz() == S_z and \
        newbasis.getnumparticles() == n_particles):
            basis_size += 1

    return basis_size

###########################################################################

#basis = createbasis(4, 4, 0)
#
#for b in basis:
#    print(b.getstate(), end = '\t')
#    b.getoccupation(0, 1)
#    print(b.phase)



a = state([1, 1, 1, 0], [0, 0, 1, 1])
print(a.getstate())
#print(a.binequiv())
#print(a.intequiv())

a.destroy(0, 1)
print(a.getstate())
print(a.phase)
a.create(0, 1)
print(a.getstate())
print(a.phase)
a.destroy(2, -1)
print(a.getstate())
print(a.phase)
a.create(2, -1)
print(a.getstate())
print(a.phase)