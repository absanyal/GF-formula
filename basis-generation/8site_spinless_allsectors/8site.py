import basisgeneration as bg
import progbar as pb
import numpy as np
import matplotlib.pyplot as plt
import time
import os
from numpy.linalg import inv
os.system('clear')

N = 12
n = 11
S_z = 0.5 * n

U = 0
Unn = 8
eta = 0.1
t = 1
tprime = t
I = complex(0, 1)

w_min = -5
w_max = Unn * (n-1) + 5
num_points = 100

w_list = np.linspace(w_min, w_max, num_points)

def getunn(s1):
    unn = 0
    c = s1.upconfig[:]
    for i in range(len(c)-1):
        if (c[i] == 1 and c[i+1] == 1):
            unn += Unn
    return unn


def mel(state1, state2):
    # calculate the hopping to right
    term = 0
    for sigma in [-1, 1]:
        for i in range(N):
            for j in range(N):
                if (i != j and abs(i - j) == 1):
                    s2 = bg.clonestate(state2)
                    s2.move(i, j, sigma)
                    termtemp = bg.innerproduct(state1, s2)
                    # print(termtemp)
                    term += termtemp
#                    print((i, j, sigma), state1.getstate(), state2.getstate(),
#                          s2.getstate(), termtemp, sep = '\t')
#        print("*" * 80)

    return term

def overlap(basis1, basis2):
    H = np.zeros((len(basis1), len(basis2)))

    for bi in range(len(basis1)):
        for bj in range(len(basis2)):

            state1 = basis1[bi]
            state2 = basis2[bj]

            if (state1.getleftnum() == state2.getleftnum()):
                H[bi][bj] = t * mel(state1, state2)
                # H[bj][bi] = H[bi][bj]
            if (state1.getleftnum() != state2.getleftnum()):
                H[bi][bj] = tprime * mel(state1, state2)
                # H[bj][bi] = H[bi][bj]

            if (basis1 == basis2 and bi == bj):
                H[bi][bj] += (getunn(state1))
    
    H = np.array(H)
    return H

def z(omega):
    return omega + I * eta


def G(H, omega):
    return np.linalg.inv(z(omega)
                         * np.eye(np.shape(H)[0], dtype=np.complex) - H)

def tmm(A, B, C):
    return np.dot(A, np.dot(B, C))

# Generate required sectors
sectorlist = []

basis = bg.createbasis(N, n, S_z)

for l_n in range(0, n+1):
    l_Sz = 0.5 * l_n
    sector = bg.createsubbasis(basis, l_n, l_Sz)
    sectorlist.append(sector)

# for sector in sectorlist:
#     for s in sector:
#         print(s.getstate(), '\t', getunn(s))
#     print('*'*50)

b1 = sectorlist[1]

H_list = []

for sector in sectorlist:
    H = overlap(sector, sector)
    H_list.append(H)

# print( overlap( sectorlist[1], sectorlist[0] ) )
# print(np.shape( overlap( sectorlist[1], sectorlist[0] ) ))

# print(overlap(sectorlist[2], sectorlist[3]))
# exit()


A_list = []
for w in w_list:

    g_list = []
    for H in H_list:
        g_list.append(G(H, w))
    
    # print(len(H_list))
    # print(len(g_list))
    
    lg_list = []
    lg_list.append(g_list[0])

    i = 1
    while (i < len(g_list)):
        # print(i)
        g_ii = g_list[i]
        g_im1_im1 = g_list[i-1]
        tau_i_im1 = overlap(sectorlist[i], sectorlist[i-1])
        tau_im1_i = np.transpose(tau_i_im1)
        lg_ii = inv( inv(g_ii) - tmm( tau_i_im1, g_im1_im1, tau_im1_i ) )
        lg_list.append(lg_ii)
        i += 1
    
    rg_list = []
    rg_list.append(g_list[len(g_list)-1])

    i = len(g_list) - 2
    while (i >= 0):
        g_ii = g_list[i]
        g_ip1_ip1 = g_list[i+1]
        tau_i_ip1 = overlap(sectorlist[i], sectorlist[i+1])
        tau_ip1_i = np.transpose(tau_i_ip1)
        rg_ii = inv( inv(g_ii) - tmm( tau_i_ip1, g_ip1_ip1, tau_ip1_i ) )
        rg_list.append(rg_ii)
        i -= 1
    
    rg_list.reverse()

    # print(g_list[0], lg_list[0])
    # exit()

    fG_list = []
    i = 0
    while (i < len(g_list)):
        # print(i)
        term1 = inv(g_list[i])

        if (i > 0):
            lg_im1_im1 = lg_list[i-1]
            tau_i_im1 = overlap(sectorlist[i], sectorlist[i-1])
            tau_im1_i = np.transpose(tau_i_im1)
            term2 = tmm( tau_i_im1,  lg_im1_im1, tau_im1_i)
        else:
            term2 = np.zeros((len(term1), len(term1)))
        
        if (i < len(g_list)-1):
            rg_ip1_ip1 = rg_list[i+1]
            tau_i_ip1 = overlap(sectorlist[i], sectorlist[i+1])
            tau_ip1_i = np.transpose(tau_i_ip1)
            term3 = tmm( tau_i_ip1,  rg_ip1_ip1, tau_ip1_i)
        else:
            term3 = np.zeros((len(term1), len(term1)))

        fG = inv(term1 - term2 - term3)
        fG_list.append(fG)
        

        i += 1
    A = 0
    for fG in fG_list:
        A += - (1/(np.pi) * (1/len(basis)) * np.imag( np.trace (fG)))
    
    A_list.append(A)

    print(w, A)

plt.plot(w_list, A_list)
fname = str(N) + '_' + str(n) + '_' + str(int(Unn * 100))
plt.savefig(fname+'.pdf')
plt.clf()

# plt.show()