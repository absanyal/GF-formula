import numpy as np
from numpy import dot, transpose, shape, pi, imag, trace, array, block, round
from numpy.linalg import inv
import matplotlib.pyplot as plt

I = np.complex(0, 1)
eta = 0.1

u = 8


def z(w):
    return w + complex(0, 1) * eta


def G(H, w):
    return inv(z(w) * np.eye(len(H)) - H)


def tmm(a, b, c):
    return dot(a, dot(b, c))


w_list = np.linspace(4, 27, 1000)
# w_list = [0]
A = []

for w in w_list:
    # print(w)
    H_5a4 = np.array([
        [3*u]
    ])
    g_5a4 = G(H_5a4, w)
    ########################################################################

    H_4b3 = np.array([
        [2*u]
    ])

    H_4d3 = np.array([
        [2*u, 1, 0],
        [1, 2*u, 1],
        [0, 1, 3*u]
    ])

    g_4b3 = G(H_4b3, w)
    g_4d3 = G(H_4d3, w)

    U_4b3_4d3 = np.array([
        [1, 0, 0]
    ])

    U_4d3_4b3 = transpose(U_4b3_4d3)

    G_4b3_4b3 = inv(inv(g_4b3) - tmm(U_4b3_4d3, g_4d3, U_4d3_4b3))
    G_4d3_4d3 = inv(inv(g_4d3) - tmm(U_4d3_4b3, g_4b3, U_4b3_4d3))
    G_4b3_4d3 = tmm(G_4b3_4b3, U_4b3_4d3, g_4d3)
    G_4d3_4b3 = tmm(G_4d3_4d3, U_4d3_4b3, g_4b3)

    g_5c4 = block([
        [G_4b3_4b3, G_4b3_4d3],
        [G_4d3_4b3, G_4d3_4d3]
    ])

    #####################################################################

    H_4a3 = np.array([
        [2*u]
    ])

    H_4c3 = np.array([
        [u, 1, 0],
        [1, u, 1],
        [0, 1, 2*u]
    ])

    g_4a3 = G(H_4a3, w)
    g_4c3 = G(H_4c3, w)

    U_4a3_4c3 = np.array([
        [1, 0, 0]
    ])

    U_4c3_4a3 = transpose(U_4a3_4c3)

    G_4a3_4a3 = inv(inv(g_4a3) - tmm(U_4a3_4c3, g_4c3, U_4c3_4a3))
    G_4c3_4c3 = inv(inv(g_4c3) - tmm(U_4c3_4a3, g_4a3, U_4a3_4c3))
    G_4a3_4c3 = tmm(G_4a3_4a3, U_4a3_4c3, g_4c3)
    G_4c3_4a3 = tmm(G_4c3_4c3, U_4c3_4a3, g_4a3)

    g_5b3 = block([
        [G_4a3_4a3, G_4a3_4c3],
        [G_4c3_4a3, G_4c3_4c3]
    ])

    ####################################################################

    H_4b2 = np.array([
        [2*u, 1, 0],
        [1, u, 1],
        [0, 1, 2*u]
    ])

    g_4b2 = G(H_4b2, w)

    H_4d2 = np.array([
        [2*u, 1, 0],
        [1, 2*u, 1],
        [0, 1, 3*u]
    ])

    g_4d2 = G(H_4d2, w)

    U_4b2_4d2 = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0]
    ])

    U_4d2_4b2 = transpose(U_4b2_4d2)

    G_4b2_4b2 = inv(inv(g_4b2) - tmm(U_4b2_4d2, g_4d2, U_4d2_4b2))
    G_4d2_4d2 = inv(inv(g_4d2) - tmm(U_4d2_4b2, g_4b2, U_4b2_4d2))
    G_4b2_4d2 = tmm(G_4b2_4b2, U_4b2_4d2, g_4d2)
    G_4d2_4b2 = tmm(G_4d2_4d2, U_4d2_4b2, g_4b2)

    g_5d3 = block([
        [G_4b2_4b2, G_4b2_4d2],
        [G_4d2_4b2, G_4d2_4d2]
    ])

    #####################################################################

    U_5a4_5c4 = np.array([
        [1, 0, 0, 0]
    ])

    U_5c4_5a4 = transpose(U_5a4_5c4)

    G_5a4_5a4 = inv(inv(g_5a4) - tmm(U_5a4_5c4, g_5c4, U_5c4_5a4))
    G_5c4_5c4 = inv(inv(g_5c4) - tmm(U_5c4_5a4, g_5a4, U_5a4_5c4))
    G_5a4_5c4 = tmm(G_5a4_5a4, U_5a4_5c4, g_5c4)
    G_5c4_5a4 = tmm(G_5c4_5c4, U_5c4_5a4, g_5a4)

    g_6a4 = block([
        [G_5a4_5a4, G_5a4_5c4],
        [G_5c4_5a4, G_5c4_5c4]
    ])

    #####################################################################

    U_5b3_5d3 = np.array([
        [0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0]
    ])

    U_5d3_5b3 = transpose(U_5b3_5d3)

    G_5b3_5b3 = inv(inv(g_5b3) - tmm(U_5b3_5d3, g_5d3, U_5d3_5b3))
    G_5d3_5d3 = inv(inv(g_5d3) - tmm(U_5d3_5b3, g_5b3, U_5b3_5d3))
    G_5b3_5d3 = tmm(G_5b3_5b3, U_5b3_5d3, g_5d3)
    G_5d3_5b3 = tmm(G_5d3_5d3, U_5d3_5b3, g_5b3)

    g_6c4 = block([
        [G_5b3_5b3, G_5b3_5d3],
        [G_5d3_5b3, G_5d3_5d3]
    ])

    #####################################################################
    U_6a4_6c4 = np.array([
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
    ])

    U_6c4_6a4 = transpose(U_6a4_6c4)

    G_6a4_6a4 = inv(inv(g_6a4) - tmm(U_6a4_6c4, g_6c4, U_6c4_6a4))
    G_6c4_6c4 = inv(inv(g_6c4) - tmm(U_6c4_6a4, g_6a4, U_6a4_6c4))
    # G_6a4_6c4 = tmm(G_6a4_6a4, U_6a4_6c4, g_6c4)
    # G_6c4_6a4 = tmm(G_6c4_6c4, U_6c4_6a4, g_6a4)

    A.append(-1/(pi * 15) * imag(trace(G_6a4_6a4) + trace(G_6c4_6c4)))
#########################################################################

plt.plot(w_list, A)
plt.show()

f = open('dosformula.dat', 'w')
for i in range(len(w_list)):
    f.write(str(w_list[i]) + '\t' + str(A[i]) + '\n')
f.close()

# H = block([
#     [H_4b2, U_4b2_4d2],
#     [U_4d2_4b2, H_4d2]
# ])

# # print(H)
# print(round(G(H, w) - g_5d3))
