import numpy as np
from numpy import dot, transpose, shape, pi, imag, trace, array, block, round
from numpy.linalg import inv
import matplotlib.pyplot as plt

I = np.complex(0, 1)
eta = 0.1

u = 0


def z(w):
    return w + complex(0, 1) * eta


def G(H, w):
    return inv(z(w) * np.eye(len(H)) - H)


def tmm(a, b, c):
    return dot(a, dot(b, c))


w_list = np.linspace(-6, 6, 1000)
# w_list = [0]
A = []

for w in w_list:

    H_4a2 = np.array([
        [u, 1, 0],
        [1, 0, 1],
        [0, 1, u]
    ])

    H_4c2 = np.array([
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, u]
    ])

    U_4a2_4c2 = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0]
    ])

    U_4c2_4a2 = transpose(U_4a2_4c2)

    g_4a2 = G(H_4a2, w)
    g_4c2 = G(H_4c2, w)

    G_4a2_4a2 = inv(inv(g_4a2) - tmm(U_4a2_4c2, g_4c2, U_4c2_4a2))
    G_4c2_4c2 = inv(inv(g_4c2) - tmm(U_4c2_4a2, g_4a2, U_4a2_4c2))
    G_4a2_4c2 = tmm(G_4a2_4a2, U_4a2_4c2, g_4c2)
    G_4c2_4a2 = tmm(G_4c2_4c2, U_4c2_4a2, g_4a2)

    g_5a2 = block([
        [G_4a2_4a2, G_4a2_4c2],
        [G_4c2_4a2, G_4c2_4c2]
    ])

    ######################################################################

    H_4b1 = np.array([
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0]
    ])

    H_4d1 = np.array([
        [u]
    ])

    g_4b1 = G(H_4b1, w)
    g_4d1 = G(H_4d1, w)

    U_4b1_4d1 = np.array([
        [0],
        [0],
        [1]
    ])

    U_4d1_4b1 = transpose(U_4b1_4d1)

    G_4b1_4b1 = inv(inv(g_4b1) - tmm(U_4b1_4d1, g_4d1, U_4d1_4b1))
    G_4d1_4d1 = inv(inv(g_4d1) - tmm(U_4d1_4b1, g_4b1, U_4b1_4d1))
    G_4b1_4d1 = tmm(G_4b1_4b1, U_4b1_4d1, g_4d1)
    G_4d1_4b1 = tmm(G_4d1_4d1, U_4d1_4b1, g_4b1)

    g_5c2 = block([
        [G_4b1_4b1, G_4b1_4d1],
        [G_4d1_4b1, G_4d1_4d1]
    ])

    #####################################################################

    H_4a1 = np.array([
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0]
    ])

    H_4c1 = np.array([
        [0]
    ])

    g_4a1 = G(H_4a1, w)
    g_4c1 = G(H_4c1, w)

    U_4a1_4c1 = np.array([
        [0],
        [0],
        [1]
    ])

    U_4c1_4a1 = transpose(U_4b1_4d1)

    G_4a1_4a1 = inv(inv(g_4a1) - tmm(U_4a1_4c1, g_4c1, U_4c1_4a1))
    G_4c1_4c1 = inv(inv(g_4c1) - tmm(U_4c1_4a1, g_4a1, U_4a1_4c1))
    G_4a1_4c1 = tmm(G_4a1_4a1, U_4a1_4c1, g_4c1)
    G_4c1_4a1 = tmm(G_4c1_4c1, U_4c1_4a1, g_4a1)

    g_5b1 = block([
        [G_4a1_4a1, G_4a1_4c1],
        [G_4c1_4a1, G_4c1_4c1]
    ])

    #####################################################################

    H_5d1 = np.array([
        [u]
    ])

    g_5d1 = G(H_5d1, w)

    #####################################################################

    U_5a2_5c2 = np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0]
    ])

    U_5c2_5a2 = transpose(U_5a2_5c2)

    G_5a2_5a2 = inv(inv(g_5a2) - tmm(U_5a2_5c2, g_5c2, U_5c2_5a2))
    G_5c2_5c2 = inv(inv(g_5c2) - tmm(U_5c2_5a2, g_5a2, U_5a2_5c2))
    G_5a2_5c2 = tmm(G_5a2_5a2, U_5a2_5c2, g_5c2)
    G_5c2_5a2 = tmm(G_5c2_5c2, U_5c2_5a2, g_5a2)

    g_6a2 = block([
        [G_5a2_5a2, G_5a2_5c2],
        [G_5c2_5a2, G_5c2_5c2]
    ])

    #####################################################################

    U_5b1_5d1 = np.array([
        [0],
        [0],
        [0],
        [1],
    ])

    U_5d1_5b1 = transpose(U_5b1_5d1)

    G_5b1_5b1 = inv(inv(g_5b1) - tmm(U_5b1_5d1, g_5d1, U_5d1_5b1))
    G_5d1_5d1 = inv(inv(g_5d1) - tmm(U_5d1_5b1, g_5b1, U_5b1_5d1))
    G_5b1_5d1 = tmm(G_5b1_5b1, U_5b1_5d1, g_5d1)
    G_5d1_5b1 = tmm(G_5d1_5d1, U_5d1_5b1, g_5b1)

    g_6c2 = block([
        [G_5b1_5b1, G_5b1_5d1],
        [G_5d1_5b1, G_5d1_5d1]
    ])

    #####################################################################
    U_6a2_6c2 = np.array([
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0],
        [0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0]
    ])

    U_6c2_6a2 = transpose(U_6a2_6c2)

    G_6a2_6a2 = inv(inv(g_6a2) - tmm(U_6a2_6c2, g_6c2, U_6c2_6a2))
    G_6c2_6c2 = inv(inv(g_6c2) - tmm(U_6c2_6a2, g_6a2, U_6a2_6c2))
    # G_6a2_6c2 = tmm(G_6a2_6a2, U_6a2_6c2, g_6c2)
    # G_6c2_6a2 = tmm(G_6c2_6c2, U_6c2_6a2, g_6a2)

    # g_6a2 = block([
    #     [G_6a2_6a2, G_6a2_6c2],
    #     [G_6c2_6a2, G_6c2_6c2]
    # ])

    A.append(-1/(pi * 15) * imag(trace(G_6a2_6a2) + trace(G_6c2_6c2)))
#########################################################################

plt.plot(w_list, A)
plt.show()

f = open('dosformula.dat', 'w')
for i in range(len(w_list)):
    f.write(str(w_list[i]) + '\t' + str(A[i]) + '\n')
f.close()

# H = block([
#     [H_4a1, U_4a1_4c1],
#     [U_4c1_4a1, H_4c1]
# ])

# # print(H)
# print(round(G(H, w) - g_5b1))
