# Thu Feb 7 12:52:23 IST 2019

import numpy as np
import os
import statemanip as sm
from numpy.linalg import inv
# from numpy import imag, trace
import matplotlib.pyplot as plt

os.system('rm *.dat')
os.system('clear')

num_sites = 6
num_particles = 3
smax = num_sites-1
ns = num_particles
t = np.complex(1, 0)

wmax = -5
wmin = 5

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
A_list = []
w_list = np.linspace(wmin, wmax, 1000)
# w_list = [0]  # For testing

for w in w_list:

    # if (w % 10 == 0):
    print(w)

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

    level = 3

    while(level <= smax):
        # while(level <= 3):  # Testing line

        # Grouping
        states = np.loadtxt('stateordinatesatlevel' + str(level) + '.dat')
        n = len(states)
        # print(level, n)
        i = 0
        while(i < n):
            line1 = states[i]
            ts1 = line1[0]
            tns1 = line1[1]
            ta1 = line1[2]
            tb1 = line1[3]
            tsize1 = line1[4]
            tpos1 = int(line1[5])

            if (i < n-1):
                line2 = states[i+1]
                ts2 = line2[0]
                tns2 = line2[1]
                ta2 = line2[2]
                tb2 = line2[3]
                tsize2 = line2[4]
                tpos2 = int(line2[5])

            s1 = sm.state(ts1, tns1, ta1, tb1)
            s2 = sm.state(ts2, tns2, ta2, tb2)

            tl1 = sm.getletter(s1)
            tl2 = sm.getletter(s2)

            if (i != n):
                if (tl1 == 1 and tl2 == 3):
                    # print(sm.getlstate(s1), sm.getlstate(s2), sep='\t')

                    f = 'g_'+str(level)+'_'+str(tpos1)+'.dat'
                    # print(f)
                    g11 = np.loadtxt(f, dtype=np.complex)
                    if (tsize1 == 1):
                        g11 = [[g11]]

                    f = 'g_'+str(level)+'_'+str(tpos2)+'.dat'
                    # print(f)
                    g22 = np.loadtxt(f, dtype=np.complex)
                    if (tsize2 == 1):
                        g22 = [[g22]]

                    u12 = sm.connecting_u(s1, s2)
                    u21 = np.transpose(u12)

                    fg11 = inv(inv(g11) - sm.tmm(u12, g22, u21))
                    fg22 = inv(inv(g22) - sm.tmm(u21, g11, u12))
                    fg12 = sm.tmm(fg11, u12, g22)
                    fg21 = sm.tmm(fg22, u21, g11)

                    fg = np.block([[fg11, fg12], [fg21, fg22]])

                    f = 'g_'+str(level+1)+'_'+str(tpos1)+'.dat'
                    np.savetxt(f, fg, fmt='%1.8f')

                    i += 1
                elif (tl1 == 2 and tl2 == 4):
                    # print(sm.getlstate(s1), sm.getlstate(s2), sep='\t')

                    f = 'g_'+str(level)+'_'+str(tpos1)+'.dat'
                    # print(f)
                    g11 = np.loadtxt(f, dtype=np.complex)
                    if (tsize1 == 1):
                        g11 = [[g11]]

                    f = 'g_'+str(level)+'_'+str(tpos2)+'.dat'
                    # print(f)
                    g22 = np.loadtxt(f, dtype=np.complex)
                    if (tsize2 == 1):
                        g22 = [[g22]]

                    u12 = sm.connecting_u(s1, s2)
                    u21 = np.transpose(u12)

                    fg11 = inv(inv(g11) - sm.tmm(u12, g22, u21))
                    fg22 = inv(inv(g22) - sm.tmm(u21, g11, u12))
                    fg12 = sm.tmm(fg11, u12, g22)
                    fg21 = sm.tmm(fg22, u21, g11)

                    fg = np.block([[fg11, fg12], [fg21, fg22]])

                    f = 'g_'+str(level+1)+'_'+str(tpos1)+'.dat'
                    np.savetxt(f, fg, fmt='%1.8f')

                    i += 1
                else:
                    f = 'g_'+str(level)+'_'+str(tpos1)+'.dat'
                    # print(f)
                    g11 = np.loadtxt(f, dtype=np.complex)
                    if (tsize1 == 1):
                        g11 = [[g11]]

                    fg = g11.copy()
                    f = 'g_'+str(level+1)+'_'+str(tpos1)+'.dat'
                    np.savetxt(f, fg, fmt='%1.8f')

                    # print(sm.getlstate(s1))

            i += 1
        # print('*'*50)

        level += 1

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

plt.plot(w_list, A_list)
plt.show()
