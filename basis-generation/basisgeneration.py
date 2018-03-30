# -*- coding: utf-8 -*-
"""
Created on Sat Feb  3 17:10:25 2018

@author: AB Sanyal
"""
#import itertools

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
                p += " " + chr(8593) + " "
            if (self.upconfig[i] == 0 and self.downconfig[i] == 1):
                p += " " + chr(8595) + " "
            if (self.upconfig[i] == 1 and self.downconfig[i] == 1):
                p += " " + chr(8645) + " "
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

    #Definition of the c, c_dagger and n operators

    #Creation operator
    def create(self, site, sigma):
        if (sigma == 1 and self.upconfig[site] == 0):

            n1 = 0
            for i in range(site + 1, self.N):
                n1 += self.upconfig[i]

            n2 = 0
            for i in range(site + 1, self.N):
                n2 += self.downconfig[i]

            self.phase *= pow(-1, n1 + n2)
            self.upconfig[site] = 1

        elif (sigma == -1 and self.downconfig[site] == 0):

            n1 = 0
            for i in range(site, self.N):
                n1 += self.upconfig[i]

            n2 = 0
            for i in range(site + 1, self.N):
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
            for i in range(site + 1, self.N):
                n1 += self.upconfig[i]

            n2 = 0
            for i in range(site + 1, self.N):
                n2 += self.downconfig[i]

            self.phase *= pow(-1, n1 + n2)
            self.upconfig[site] = 0

        elif (sigma == -1 and self.downconfig[site] == 1):

            n1 = 0
            for i in range(site, self.N):
                n1 += self.upconfig[i]

            n2 = 0
            for i in range(site, self.N):
                n2 += self.downconfig[i]

            self.phase *= pow(-1, n1 + n2 - 1)
            self.downconfig[site] = 0

        elif (sigma == 1 and self.upconfig[site] == 0):
            self.upconfig = [0 for i in range(self.N)]
            self.downconfig = [0 for i in range(self.N)]
            self.phase = 0

        elif (sigma == -1 and self.downconfig[site] == 0):
            self.upconfig = [0 for i in range(self.N)]
            self.downconfig = [0 for i in range(self.N)]
            self.phase = 0

    #Move a particle from site to site
    def move(self, i, j, sigma):
        """
        Move a particle of spin sigma = +/- 1 from site i to site j
        """
        if (self.getoccupation(i, sigma) != 0):
            self.destroy(i, sigma)
            self.create(j, sigma)
        else:
            self.upconfig = [0 for i in range(self.N)]
            self.downconfig = [0 for i in range(self.N)]
            self.phase = 0


    #Returns the total number of particles in the state
    def getnumparticles(self):
        n = 0
        for i in self.upconfig:
            n += i
        for j in self.downconfig:
            n += j
        return n

    #Get the number of particles on the left block
    def getleftnum(self):
        n = 0
        for i in range(int(self.N / 2)):
            n += self.upconfig[i] + self.downconfig[i]
        return n

    #Get the total Sz of particles on the left block
    def getleftSz(self):
        s = 0
        for i in range(int(self.N / 2)):
            s += 0.5 * self.upconfig[i] - 0.5 * self.downconfig[i]
        return s

    #Returns the total spin of the state
    def getSz(self):
        s = 0
        for i in self.upconfig:
            s += i * 1/2
        for j in self.downconfig:
            s += j * (-1/2)
        return s

    #Occupation number operator for a given site and sigma
    def getoccupation(self, site, sigma):
        self.destroy(site, sigma)
        self.create(site, sigma)

        return self.phase

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

def createbasis(N, n_particles, S_z = 0):
    basis = []

    filledlist = [0 for i in range(2 * N)]

    for i in range(n_particles):
        filledlist[i] = 1

    filledlist.reverse()

    #make max sector
    upsector = filledlist[:n_particles]
    downsector = filledlist[n_particles:]

    maxstate = state(upsector, downsector)
    #print(maxstate.getstate())
    maxint = maxstate.intequiv()
    #print(maxint)

    #make min sector
    filledlist.reverse()

    upsector = filledlist[:n_particles]
    downsector = filledlist[n_particles:]

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


#Create a sub-basis from a given basis with left-sector specified
def createsubbasis(basis, l_n, l_Sz = 0):
    sbasis = []

    for state in basis:
        if (state.getleftnum() == l_n and state.getleftSz() == l_Sz):
            sbasis.append(state)

    return sbasis

#Create a basis with specified parameters and left block parameters
#lbsbasis = Left Block Specified Basis
def createlbsbasis(N, n_particles, S_z, l_n, l_Sz):

    basis = []

    filledlist = [0 for i in range(2 * N)]

    for i in range(n_particles):
        filledlist[i] = 1

    filledlist.reverse()

    #make max sector
    upsector = filledlist[:n_particles]
    downsector = filledlist[n_particles:]

    maxstate = state(upsector, downsector)
    #print(maxstate.getstate())
    maxint = maxstate.intequiv()
    #print(maxint)

    #make min sector
    filledlist.reverse()

    upsector = filledlist[:n_particles]
    downsector = filledlist[n_particles:]

    minstate = state(upsector, downsector)
    #print(minstate.getstate())
    minint = minstate.intequiv()
    #print(minint)

    #make all states and select the appropriate ones
    for i in range(minint, maxint + 1):
        newbasis = makestatefromint(N, i)
        if (
            newbasis.getSz() == S_z and \
            newbasis.getnumparticles() == n_particles and \
            newbasis.getleftnum() == l_n and \
            newbasis.getleftSz() == l_Sz
        ):
            basis.append(newbasis)

    return basis

#create a Fock Space basis with a specified left block parameters
#lfs = Left Fock Space

def createlfsbasis(N, n_particles, S_z, l_n):
    import numpy as np
    basis = []
    spins = list(0.5 * np.array(range(-int(l_n), int(l_n)+1)))
    for s in spins:
        tempbasis = createlbsbasis(N, n_particles, S_z, l_n, s)
        basis += tempbasis
    
    return basis

#Inner product of two states
def innerproduct(a, b):
    if (a.getstate() == b.getstate()):
        return a.phase * b.phase
    else:
        return 0

#Clone a state
def clonestate(a):
    upc = a.upconfig[:]
    downc = a.downconfig[:]
    return state(upc, downc, a.phase)

###########################################################################
