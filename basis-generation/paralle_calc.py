# Mon Apr 9 07:39:18 IST 2018

import os
import time
import multiprocessing

import basisgeneration as bg

os.system('cls')
os.system('clear')

N = 8
n = 8
Sz = (0.5 * n) % 1

numlist = list(range(0, n + 1))
numlist.reverse()

if __name__ == '__main__':

    # Unparallelized method
    # tup1 = time.perf_counter()
    # basis = []
    # for l_n in numlist:
    #     print("Generating basis for n_l=", l_n)
    #     temp_basis = bg.createlfsbasis(N, n, Sz, l_n)
    #     basis += temp_basis
    # tup2 = time.perf_counter()

    # print("Time taken without parallelization:", tup2 - tup1)

    # Parallelized method

    tp1 = time.perf_counter()

    pool = multiprocessing.Pool()
    args = [(N, n, Sz, i) for i in numlist]
    #numlist = tuple(numlist)
    listofbasis = pool.starmap(
        bg.createlfsbasis, args)

    # basis = []
    # for sb in listofbasis:
    #     basis += sb

    # i = 0
    # for s in basis:
    #     print(i, s.getstate(), (s.getleftnum(), s.getleftSz()), sep='\t')
    #     i += 1

    tp2 = time.perf_counter()

    print("Time taken with parallelization:", tp2 - tp1)
