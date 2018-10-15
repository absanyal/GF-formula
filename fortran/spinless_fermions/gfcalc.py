import numpy as np
from numpy import dot, transpose, shape, pi, imag, trace, array
from numpy.linalg import inv
import matplotlib.pyplot as plt

I = np.complex(0, 1)
eta = 0.1


def z(w):
    return w + complex(0, 1) * eta


def G(H, w):
    return inv(z(w) * np.eye(len(H)) - H)


def tmm(a, b, c):
    return dot(a, dot(b, c))


w_list = np.linspace(-5, 5, 1000)
# w_list = [0]
A = []

H_3a2 = [0]
H_3c2 = [
    [0, 1],
    [1, 0]]
H_3b1 = [
    [0, 1],
    [1, 0]]
H_3d1 = [0]
tau_3a2_3c2 = [[1, 0]]
tau_3c2_3a2 = transpose(tau_3a2_3c2)
tau_3b1_3d1 = [[1], [0]]
tau_3d1_3b1 = transpose(tau_3b1_3d1)
tau_ac_bd = [
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0]]
tau_bd_ac = transpose(tau_ac_bd)

# print(tau_bd_ac)

# H = [
#     [0, 1, 0, 0, 0, 0],
#     [1, 0, 1, 1, 0, 0],
#     [0, 1, 0, 0, 1, 0],
#     [0, 1, 0, 0, 1, 0],
#     [0, 0, 1, 1, 0, 1],
#     [0, 0, 0, 0, 1, 0]
# ]

H = array([
    [0, 1, 0],
    [1, 0, 1],
    [0, 1, 0]
])

for w in w_list:
    fG = G(H, w)
    snf = 1/len(H)
    A.append((-snf / pi) * np.imag(np.trace(fG)))

plt.plot(w_list, A)

# print(imag(fG))

A = []

# for w in w_list:

#     # print(w)

#     g_3a2 = G(H_3a2, w)
#     g_3c2 = G(H_3c2, w)

#     lG_3a2 = g_3a2
#     lG_3c2 = inv(inv(g_3c2) - tmm(tau_3c2_3a2, lG_3a2, tau_3a2_3c2))

#     rG_3c2 = g_3c2
#     rG_3a2 = inv(inv(g_3a2) - tmm(tau_3a2_3c2, lG_3c2, tau_3c2_3a2))

#     fG_3a2 = inv(inv(g_3a2) - tmm(tau_3a2_3c2, rG_3c2, tau_3c2_3a2))
#     fG_3c2 = inv(inv(g_3c2) - tmm(tau_3c2_3a2, lG_3a2, tau_3a2_3c2))

#     fG_tau_3a2_3c2 = tmm(fG_3a2, tau_3a2_3c2, rG_3c2)
#     fG_tau_3c2_3a2 = tmm(rG_3c2, tau_3c2_3a2, fG_3a2)

#     # print( fG_tau_3a2_3c2 == transpose(fG_tau_3c2_3a2) )

#     fG_ac = np.block(
#         [[fG_3a2, fG_tau_3a2_3c2],
#          [fG_tau_3c2_3a2, fG_3c2]]
#     )

#     g_3b1 = G(H_3b1, w)
#     g_3d1 = G(H_3d1, w)

#     lG_3b1 = g_3b1
#     lG_3d1 = inv(inv(g_3d1) - tmm(tau_3d1_3b1, lG_3b1, tau_3b1_3d1))

#     rG_3d1 = g_3d1
#     rG_3b1 = inv(inv(g_3b1) - tmm(tau_3b1_3d1, lG_3d1, tau_3d1_3b1))

#     fG_3b1 = inv(inv(g_3b1) - tmm(tau_3b1_3d1, rG_3d1, tau_3d1_3b1))
#     fG_3d1 = inv(inv(g_3d1) - tmm(tau_3d1_3b1, lG_3b1, tau_3b1_3d1))

#     fG_tau_3b1_3d1 = tmm(lG_3b1, tau_3b1_3d1, fG_3d1)
#     fG_tau_3d1_3b1 = tmm(fG_3d1, tau_3d1_3b1, lG_3b1)

#     # print( fG_tau_3b1_3d1 == transpose(fG_tau_3d1_3b1) )

#     fG_bd = np.block(
#         [[fG_3b1, fG_tau_3b1_3d1],
#          [fG_tau_3d1_3b1, fG_3d1]]
#     )

#     g_ac = fG_ac
#     g_bd = fG_bd

#     lG_ac = g_ac
#     lG_bd = inv(inv(g_bd) - tmm(tau_bd_ac, lG_ac, tau_ac_bd))

#     rG_bd = g_bd
#     rG_ac = inv(inv(g_ac) - tmm(tau_ac_bd, lG_bd, tau_bd_ac))

#     fG_ac = inv(inv(g_ac) - tmm(tau_ac_bd, rG_bd, tau_bd_ac))
#     fG_bd = inv(inv(g_bd) - tmm(tau_bd_ac, lG_ac, tau_ac_bd))

#     snf = 1 / (len(fG_ac) + len(fG_bd))

#     A.append((-snf / pi) * np.imag(np.trace(fG_ac) + np.trace(fG_bd)))

# f = open('dos_smart.dat', 'w')

# for i in range(len(w_list)):
#     f.write(str(w_list[i]) + '\t' + str(A[i]) + '\n')

# f.close()

# plt.plot(w_list, A)

# # plt.show()

# print(fG_ac[1:3, 1:3])

H1 = array([[0]])
H2 = array([[0]])
H3 = array([[0]])

tau12 = array([[0]])
tau21 = array([[0]])
tau23 = array([[0]])
tau32 = array([[0]])

for w in w_list:

    g1 = G(H1, w)
    g2 = G(H2, w)
    g3 = G(H3, w)

    lg1 = inv(inv(g1))
    lg2 = inv(inv(g2) - tmm(tau21, g1, tau12))
    lg3 = inv(inv(g3) - tmm(tau32, g2, tau23))

    rg3 = inv(inv(g3))
    rg2 = inv(inv(g2) - tmm(tau23, g3, tau32))
    rg1 = inv(inv(g1) - tmm(tau12, g2, tau21))

    fG1 = inv(inv(g1) - tmm(tau12, rg2, tau21))
    fG2 = inv(inv(g2) - tmm(tau21, lg1, tau12) - tmm(tau23, rg3, tau32))
    fG3 = inv(inv(g3) - tmm(tau32, lg2, tau23))

    A.append((1/3) * (-1/pi) * imag(trace(fG1) + trace(fG2) + trace(fG2)))

plt.plot(w_list, A)

plt.show()
