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


w_list = np.linspace(-5, 19, 5000)
# w_list = [0]
A = []

for w in w_list:
    H_3c3 = np.array([[2*u]])
    g_3c3 = G(H_3c3, w)

    g_4a3 = g_3c3[:]

    H_3b2 = np.array([[u]])
    g_3b2 = G(H_3b2, w)

    H_3d2 = np.array([
        [u, 1],
        [1, 2*u]])
    g_3d2 = G(H_3d2, w)

    U_3b2_3d2 = np.array([[1, 0]])
    U_3d2_3b2 = transpose(U_3b2_3d2)

    G_3b2_3b2 = inv(inv(g_3b2) - tmm(U_3b2_3d2, g_3d2, U_3d2_3b2))
    G_3d2_3d2 = inv(inv(g_3d2) - tmm(U_3d2_3b2, g_3b2, U_3b2_3d2))
    G_3b2_3d2 = tmm(G_3b2_3b2, U_3b2_3d2, g_3d2)
    G_3d2_3b2 = tmm(G_3d2_3d2, U_3d2_3b2, g_3b2)

    g_4c3 = block([[G_3b2_3b2, G_3b2_3d2], [G_3d2_3b2, G_3d2_3d2]])

    U_4a3_4c3 = np.array([[1, 0, 0]])
    U_4c3_4a3 = transpose(U_4a3_4c3)

    G_4a3_4a3 = inv(inv(g_4a3) - tmm(U_4a3_4c3, g_4c3, U_4c3_4a3))
    G_4c3_4c3 = inv(inv(g_4c3) - tmm(U_4c3_4a3, g_4a3, U_4a3_4c3))
    G_4a3_4c3 = tmm(G_4a3_4a3, U_4a3_4c3, g_4c3)
    G_4c3_4a3 = tmm(G_4c3_4c3, U_4c3_4a3, g_4a3)

    g_5a3 = block([
        [G_4a3_4a3, G_4a3_4c3],
        [G_4c3_4a3, G_4c3_4c3]
    ])
    #########################################################################

    H_4b2 = np.array([
        [u, 1, 0],
        [1, 0, 1],
        [0, 1, u]
    ])
    H_4d2 = np.array([
        [u, 1, 0],
        [1, u, 1],
        [0, 1, 2*u]
    ])
    U_4b2_4d2 = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0]
    ])
    U_4d2_4b2 = transpose(U_4b2_4d2)

    g_4b2 = G(H_4b2, w)
    g_4d2 = G(H_4d2, w)

    G_4b2_4b2 = inv(inv(g_4b2) - tmm(U_4b2_4d2, g_4d2, U_4d2_4b2))
    G_4d2_4d2 = inv(inv(g_4d2) - tmm(U_4d2_4b2, g_4b2, U_4b2_4d2))
    G_4b2_4d2 = tmm(G_4b2_4b2, U_4b2_4d2, g_4d2)
    G_4d2_4b2 = tmm(G_4d2_4d2, U_4d2_4b2, g_4b2)

    g_5c3 = np.block([
        [G_4b2_4b2, G_4b2_4d2],
        [G_4d2_4b2, G_4d2_4d2]
    ])

    #########################################################################
    U_5a3_5c3 = block([
        [np.zeros((1, 3)), np.zeros((1, 3))],
        [np.eye(3), np.zeros((3, 3))]
    ])

    U_5c3_5a3 = transpose(U_5a3_5c3)

    # print(U_5a3_5c3)

    G_5a3_5a3 = inv(inv(g_5a3) - tmm(U_5a3_5c3, g_5c3, U_5c3_5a3))
    G_5c3_5c3 = inv(inv(g_5c3) - tmm(U_5c3_5a3, g_5a3, U_5a3_5c3))
    G_5a3_5c3 = tmm(G_5a3_5a3, U_5a3_5c3, g_5c3)
    G_5c3_5a3 = tmm(G_5c3_5c3, U_5c3_5a3, g_5a3)

    g_6a3 = block([
        [G_5a3_5a3, G_5a3_5c3],
        [G_5c3_5a3, G_5c3_5c3]
    ])

    #########################################################################

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

    g_5b2 = np.block([
        [G_4a2_4a2, G_4a2_4c2],
        [G_4c2_4a2, G_4c2_4c2]
    ])

    #########################################################################

    H_4b1 = np.array([
        [u, 1, 0],
        [1, u, 1],
        [0, 1, u]
    ])
    g_4b1 = G(H_4b1, w)

    H_4d1 = np.array([[2*u]])
    g_4d1 = G(H_4d1, w)

    U_4b1_4d1 = np.array([
        [0],
        [0],
        [1]
    ])

    U_4d1_4b1 = transpose(U_4b1_4d1)

    # print(shape(U_4b1_4d1))

    G_4b1_4b1 = inv(inv(g_4b1) - tmm(U_4b1_4d1, g_4d1, U_4d1_4b1))
    G_4d1_4d1 = inv(inv(g_4d1) - tmm(U_4d1_4b1, g_4b1, U_4b1_4d1))
    G_4b1_4d1 = tmm(G_4b1_4b1, U_4b1_4d1, g_4d1)
    G_4d1_4b1 = tmm(G_4d1_4d1, U_4d1_4b1, g_4b1)

    g_5d2 = block([
        [G_4b1_4b1, G_4b1_4d1],
        [G_4d1_4b1, G_4d1_4d1]
    ])

    #########################################################################

    U_5b2_5d2 = block([
        [np.zeros((3, 3)), np.zeros((3, 1))],
        [np.eye(3), np.zeros((3, 1))]
    ])

    U_5d2_5b2 = transpose(U_5b2_5d2)

    # print(shape(U_5b2_5d2))

    G_5b2_5b2 = inv(inv(g_5b2) - tmm(U_5b2_5d2, g_5d2, U_5d2_5b2))
    G_5d2_5d2 = inv(inv(g_5d2) - tmm(U_5d2_5b2, g_5b2, U_5b2_5d2))
    G_5b2_5d2 = tmm(G_5b2_5b2, U_5b2_5d2, g_5d2)
    G_5d2_5b2 = tmm(G_5d2_5d2, U_5d2_5b2, g_5b2)

    g_6c3 = block([
        [G_5b2_5b2, G_5b2_5d2],
        [G_5d2_5b2, G_5d2_5d2]
    ])

    #########################################################################

    U_6a3_6c3 = block([
        [np.zeros((4, 6)), np.zeros((4, 4))],
        [np.eye(6), np.zeros((6, 4))]
    ])

    U_6c3_6a3 = transpose(U_6a3_6c3)

    G_6a3_6a3 = inv(inv(g_6a3) - tmm(U_6a3_6c3, g_6c3, U_6c3_6a3))
    G_6c3_6c3 = inv(inv(g_6c3) - tmm(U_6c3_6a3, g_6a3, U_6a3_6c3))

    A.append(-1/(pi * 20) * imag(trace(G_6a3_6a3) + trace(G_6c3_6c3)))

    #########################################################################

    # print(U_5a3_5c3)

    # H1 = np.array([
    #     [0]
    # ])

    # H1 = np.array([
    #     [0, 1, 0, 0, 0, 0],
    #     [1, 0, 1, 1, 0, 0],
    #     [0, 1, 0, 0, 1, 0],
    #     [0, 1, 0, 0, 1, 0],
    #     [0, 0, 1, 1, 0, 1],
    #     [0, 0, 0, 0, 1, 0]
    # ])

    # H = block([
    #     [H1, U_5b2_5d2],
    #     [U_5d2_5b2, H2]
    # ])

    # print(round(G(H, w) - g_6c3))

    # H = block([
    #     [np.zeros((1, 3)), np.zeros((1, 3))],
    #     [np.eye(3), np.zeros((3, 3))]
    # ])

    # print(U_5a3_5c3 - H)

plt.plot(w_list, A)
plt.show()

f = open('dosformula.dat', 'w')
for i in range(len(w_list)):
    f.write(str(w_list[i]) + '\t' + str(A[i]) + '\n')
f.close()
