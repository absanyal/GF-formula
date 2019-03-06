import numpy as np
import os
import statemanip as sm
from numpy.linalg import inv
# from numpy import imag, trace
import matplotlib.pyplot as plt
import time
from math import factorial

os.system('rm *.dat')
os.system('clear')

N_max = 20

# num_sites = 20
# N = num_sites

for num_sites in range(3, N_max + 1):
    N = num_sites
    filling_list = []
    costlist = []
    dcostlist = []
    for p in range(1, num_sites):
        num_particles = p
        smax = num_sites-1
        ns = num_particles
        t = np.complex(1, 0)
        uint = 8.0

        num_w_points = 1

        print(num_sites, 'sites,', num_particles, 'particles.')

        # dirName = 'costdata'
        # if not os.path.exists(dirName):
        #     os.mkdir(dirName)
        #     print("Directory ", dirName,  " created.")
        # else:
        #     print("Directory ", dirName,  " already exists.")


        # w_list = [0]  # For testing

        s1 = sm.state(smax+1, ns, 0, 0)
        s2 = sm.state(smax+1, ns, 1, 0)


        # Optional, splitting entire list
        f = open('splitstates.dat', 'w')

        sm.splitstates(f, s1)
        sm.splitstates(f, s2)

        f.close()

        sizelist = np.loadtxt('splitstates.dat', usecols = (1,))

        cost = sum(sizelist ** 3)
        # print(cost)

        directcost = (factorial(N) / (factorial(p) * factorial(N-p)))**3

        filling_list.append(p/N)
        costlist.append(cost)
        dcostlist.append(directcost)

    plt.plot(filling_list, costlist, label = 'Algorithm cost')
    plt.plot(filling_list, dcostlist, label = 'Direct inversion cost')
    plt.legend()
    plt.savefig(str(N)+'_site_cost.pdf')
    plt.clf()
    # plt.show()
    print('*'*50)

