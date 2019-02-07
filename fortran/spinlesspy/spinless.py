# Thu Feb 7 12:52:23 IST 2019

import numpy as np
import os
import statemanip as sm

os.system('rm *.dat')
os.system('clear')

num_sites = 4
num_particles = 2
smax = num_sites-1
ns = num_particles
t = np.complex(1, 0)

s1 = sm.state(smax, ns, 0, 0)
s2 = sm.state(smax, ns, 1, 0)
s3 = sm.state(smax, ns-1, 0, 1)
s4 = sm.state(smax, ns-1, 1, 1)

# Optional, splitting entire list
# f = open('splitstates.dat', 'w')

# sm.splitstates(f, s1)
# sm.splitstates(f, s2)
# sm.splitstates(f, s3)
# sm.splitstates(f, s4)

# f.close()

# Generate splitting at all levels
level_range = range(3, smax + 1)
# level_range = [3]  # For testing
for level in level_range:
    f = open('splitstatesatlevel' + str(level) + '.dat', 'w')

    sm.splitstatesatlevel(f, s1, level)
    sm.splitstatesatlevel(f, s2, level)
    sm.splitstatesatlevel(f, s3, level)
    sm.splitstatesatlevel(f, s4, level)

    f.close()

    f = open('splitstatesatlevelraw' + str(level) + '.dat', 'w')

    sm.splitstatesatlevelraw(f, s1, level)
    sm.splitstatesatlevelraw(f, s2, level)
    sm.splitstatesatlevelraw(f, s3, level)
    sm.splitstatesatlevelraw(f, s4, level)

    f.close()

    f = open('stateordinatesatlevel' + str(level) + '.dat', 'w')

    sizelist = np.loadtxt('splitstatesatlevelraw' +
                          str(level) + '.dat', usecols=(4,), dtype=int)
    rawlist = np.loadtxt('splitstatesatlevelraw' +
                         str(level) + '.dat', dtype=int)

    # print(sizelist)
    # print(rawlist)
    n = len(sizelist)

    for i in range(n):
        # print(sum(sizelist[:i]) + 1)
        f.write(
            str(rawlist[i, 0]) + '\t' +
            str(rawlist[i, 1]) + '\t' +
            str(rawlist[i, 2]) + '\t' +
            str(rawlist[i, 3]) + '\t' +
            str(rawlist[i, 4]) + '\t' +
            str(sum(sizelist[:i]) + 1) + '\n'
        )

    f.close()

# w_list = np.linspace(-5, 5, 1000)
w_list = [0]  # For testing

for w in w_list:

    # Generate the GF files for S = 3

    statesl3list = np.loadtxt('stateordinatesatlevel3.dat', dtype=int)
    # print(len(statesl3list))
    for line in statesl3list:
        ts = line[0]
        tns = line[1]
        ta = line[2]
        tb = line[3]
        tsize = line[4]
        tpos = line[5]
        # print(ts, tns, ta, tb, tsize, tpos)
        t_state = sm.state(ts, tns, ta, tb)
        fname = 'g_3_'+str(tpos)+'.dat'
        if (tsize == 1):
            h = np.array([[t]], dtype=complex)
        else:
            h = np.array([[0, t], [t, 0]], dtype=complex)
        mat = sm.G(h, w)
        mat = mat.reshape((tsize, tsize))
        # print(tpos)
        # print(mat)
        np.savetxt(fname, mat, fmt='%1.8f')

    ####################################################
