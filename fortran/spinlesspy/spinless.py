# Thu Feb 7 12:52:23 IST 2019

import numpy as np
import os
import statemanip as sm
from numpy.linalg import inv
# from numpy import imag, trace
import matplotlib.pyplot as plt
import time

os.system('rm *.dat')
os.system('clear')

num_sites = 12
num_particles = 6
smax = num_sites-1
ns = num_particles
t = np.complex(1, 0)

print(num_sites, 'sites,', num_particles, 'particles.')

wmax = 5
wmin = -5

A_list = []
w_list = np.linspace(wmin, wmax, 10)
# w_list = [0]  # For testing

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

t_start = time.perf_counter()

# Generate splitting at all levels
level_range = range(1, smax + 1)
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


level = smax
while (level > 3):
    f = 'stateordinatesatlevel' + str(level) + '.dat'
    ordinates = np.loadtxt(f, dtype=int)
    # print(ordinates)

    f = open('statelinkedatlevel' + str(level) + '.dat', 'w')

    for line in ordinates:
        ts = line[0]
        tns = line[1]
        ta = line[2]
        tb = line[3]
        tsize = line[4]
        tpos = line[5]

        s1 = sm.state(ts, tns, ta, tb)

        tp1 = 0
        tp2 = 0

        s1s0 = sm.relegate0(s1)
        if (sm.checkvalidity(s1s0) == 1):
            tp1 = tpos

        s1s1 = sm.relegate1(s1)
        if (sm.checkvalidity(s1s1) == 1):
            tp2 = tpos + sm.getstatesize(s1s0)

        # print(sm.getlstate(s1), tp1, tp2)
        p = str(ts) + '\t' + \
            str(tns) + '\t' + \
            str(ta) + '\t' + \
            str(tb) + '\t' + \
            str(tsize) + '\t' + \
            str(tpos) + '\t' + \
            str(tp1) + '\t' + \
            str(tp2) + '\n'
        f.write(p)

    # print('-'*50)
    f.close()
    level -= 1

# quit()

# print('*'*50)
for w in w_list:

    # if (w % 10 == 0):
    # print(w)

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
            h = np.array([[complex(0.0, 0.0)]], dtype=complex)
        else:
            h = np.array([[0, t], [t, 0]], dtype=complex)
        mat = sm.G(h, w)
        mat = mat.reshape((tsize, tsize))
        # print(tpos)
        # print(mat)
        np.savetxt(fname, mat, fmt='%1.8f')

    ####################################################

    level = 4

    while(level <= smax):
        # while(level <= 3):  # Testing line

        # Grouping
        states = np.loadtxt('statelinkedatlevel' + str(level) + '.dat')
        n = len(states)
        # print(level, n)
        i = 0
        while(i < n):
            line1 = states[i]
            ts1 = int(line1[0])
            tns1 = int(line1[1])
            ta1 = int(line1[2])
            tb1 = int(line1[3])
            tsize1 = int(line1[4])
            tpos1 = int(line1[5])
            tchild1 = int(line1[6])
            tchild2 = int(line1[7])

            s1 = sm.state(ts1, tns1, ta1, tb1)

            if (tchild1 != 0 and tchild2 != 0):

                f = 'stateordinatesatlevel' + str(level-1) + '.dat'
                childordinates = np.loadtxt(f, dtype=int)
                for line in childordinates:
                    ts = int(line[0])
                    tns = int(line[1])
                    ta = int(line[2])
                    tb = int(line[3])
                    tsize = int(line[4])
                    tpos = int(line[5])

                    if (tpos == tchild1):
                        sc1 = sm.state(ts, tns, ta, tb)
                        fc1 = 'g_'+str(level-1)+'_'+str(tpos)+'.dat'

                    if (tpos == tchild2):
                        sc2 = sm.state(ts, tns, ta, tb)
                        fc2 = 'g_'+str(level-1)+'_'+str(tpos)+'.dat'

                # print(sm.getlstate(s1),
                #       sm.getlstate(sc1),
                #       sm.getlstate(sc2))

                g11 = np.loadtxt(fc1, dtype=np.complex)
                if (np.shape(g11) == ()):
                    g11 = np.array([[g11]], dtype=complex)

                g22 = np.loadtxt(fc2, dtype=np.complex)
                if (np.shape(g22) == ()):
                    g22 = np.array([[g22]], dtype=complex)

                u12 = sm.connecting_u(sc1, sc2)
                u21 = np.transpose(u12)

                fg11 = inv(inv(g11) - sm.tmm(u12, g22, u21))
                fg22 = inv(inv(g22) - sm.tmm(u21, g11, u12))
                fg12 = sm.tmm(fg11, u12, g22)
                fg21 = sm.tmm(fg22, u21, g11)

                fg = np.block([[fg11, fg12], [fg21, fg22]])

                f = 'g_'+str(level)+'_'+str(tpos1)+'.dat'
                np.savetxt(f, fg, fmt='%1.8f')

            else:
                tchild1 = int(tchild1 + tchild2)
                f = 'stateordinatesatlevel' + str(level-1) + '.dat'
                childordinates = np.loadtxt(f, dtype=int)
                for line in childordinates:
                    ts = int(line[0])
                    tns = int(line[1])
                    ta = int(line[2])
                    tb = int(line[3])
                    tsize = int(line[4])
                    tpos = int(line[5])

                    if (tpos == tchild1):
                        sc1 = sm.state(ts, tns, ta, tb)
                        fc1 = 'g_'+str(level-1)+'_'+str(tpos)+'.dat'

                    g11 = np.loadtxt(fc1, dtype=np.complex)
                    if (np.shape(g11) == ()):
                        g11 = np.array([[g11]], dtype=complex)

                    f = 'g_'+str(level)+'_'+str(tpos1)+'.dat'
                    np.savetxt(f, g11, fmt='%1.8f')

                # print(sm.getlstate(s1),
                #       sm.getlstate(sc1))

            i += 1

        level += 1
        # print('*'*50)

    # Connect the 4 states in pairs at smax level
    states = np.loadtxt('stateordinatesatlevel' + str(smax) + '.dat')
    # Connect 1 and 2
    line1 = states[0]
    ts1 = int(line1[0])
    tns1 = int(line1[1])
    ta1 = int(line1[2])
    tb1 = int(line1[3])
    tsize1 = int(line1[4])
    tpos1 = int(line1[5])
    s1 = sm.state(ts1, tns1, ta1, tb1)
    fs1 = 'g_'+str(ts1)+'_'+str(tpos1)+'.dat'
    g11 = np.loadtxt(fs1, dtype=np.complex)
    if (np.shape(g11) == ()):
        g11 = np.array([[g11]], dtype=complex)

    # print(np.shape(g11))

    line2 = states[1]
    ts2 = int(line2[0])
    tns2 = int(line2[1])
    ta2 = int(line2[2])
    tb2 = int(line2[3])
    tsize2 = int(line2[4])
    tpos2 = int(line2[5])
    s2 = sm.state(ts2, tns2, ta2, tb2)
    fs2 = 'g_'+str(ts2)+'_'+str(tpos2)+'.dat'
    g22 = np.loadtxt(fs2, dtype=np.complex)
    if (np.shape(g22) == ()):
        g22 = np.array([[g22]], dtype=complex)
    # print(np.shape(g22))
    u12 = sm.connecting_u(s1, s2)
    u21 = np.transpose(u12)

    fg11 = inv(inv(g11) - sm.tmm(u12, g22, u21))
    fg22 = inv(inv(g22) - sm.tmm(u21, g11, u12))
    fg12 = sm.tmm(fg11, u12, g22)
    fg21 = sm.tmm(fg22, u21, g11)

    fg = np.block([[fg11, fg12], [fg21, fg22]])

    # print(np.shape(fg))

    f = 'g_'+str(smax+1)+'_'+str(tpos1)+'.dat'
    np.savetxt(f, fg, fmt='%1.8f')

    # Connect 3 and 4

    line1 = states[2]
    ts1 = int(line1[0])
    tns1 = int(line1[1])
    ta1 = int(line1[2])
    tb1 = int(line1[3])
    tsize1 = int(line1[4])
    tpos1 = int(line1[5])
    s1 = sm.state(ts1, tns1, ta1, tb1)
    fs1 = 'g_'+str(ts1)+'_'+str(tpos1)+'.dat'
    g11 = np.loadtxt(fs1, dtype=np.complex)
    if (np.shape(g11) == ()):
        g11 = np.array([[g11]], dtype=complex)

    line2 = states[3]
    ts2 = int(line2[0])
    tns2 = int(line2[1])
    ta2 = int(line2[2])
    tb2 = int(line2[3])
    tsize2 = int(line2[4])
    tpos2 = int(line2[5])
    s2 = sm.state(ts2, tns2, ta2, tb2)
    fs2 = 'g_'+str(ts2)+'_'+str(tpos2)+'.dat'
    g22 = np.loadtxt(fs2, dtype=np.complex)
    if (np.shape(g22) == ()):
        g22 = np.array([[g22]], dtype=complex)

    u12 = sm.connecting_u(s1, s2)
    u21 = np.transpose(u12)

    fg11 = inv(inv(g11) - sm.tmm(u12, g22, u21))
    fg22 = inv(inv(g22) - sm.tmm(u21, g11, u12))
    fg12 = sm.tmm(fg11, u12, g22)
    fg21 = sm.tmm(fg22, u21, g11)

    fg = np.block([[fg11, fg12], [fg21, fg22]])
    # print(np.shape(fg))

    f = 'g_'+str(smax+1)+'_'+str(tpos1)+'.dat'
    np.savetxt(f, fg, fmt='%1.8f')

    # Connect two large blocks

    s1 = sm.state(smax+1, ns, 0, 0)
    s2 = sm.state(smax+1, ns, 1, 0)

    tsize1 = sm.getstatesize(s1)
    tsize2 = sm.getstatesize(s2)

    tpos1 = 1
    tpos2 = tpos1 + sm.getstatesize(s1)

    f = 'g_'+str(smax+1)+'_'+str(tpos1)+'.dat'
    g11 = np.loadtxt(f, dtype=np.complex)

    f = 'g_'+str(smax+1)+'_'+str(tpos2)+'.dat'
    g22 = np.loadtxt(f, dtype=np.complex)

    u12 = sm.connecting_u(s1, s2)
    u21 = np.transpose(u12)

    fg11 = inv(inv(g11) - sm.tmm(u12, g22, u21))
    fg22 = inv(inv(g22) - sm.tmm(u21, g11, u12))

    A = (-1/np.pi) * 1/(tsize1 + tsize2) * \
        np.imag(np.trace(fg11) + np.trace(fg22))

    A_list.append(A)

    print(w, A)

t_stop = time.perf_counter()

print("Time taken per omega loop is",
      round((t_stop-t_start)/len(w_list), 5), 'seconds.')

# plt.plot(w_list, A_list)
# plt.show()


# print(sm.getlstate(s1))
# print(sm.getlstate(s2))
# print(sm.getlstate(s3))
# print(sm.getlstate(s4))
