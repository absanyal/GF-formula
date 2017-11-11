# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 10:26:37 2017

@author: amit
"""

# imports
from numpy import exp, pi, imag
import matplotlib.pyplot as plt

# constants
eta = 0.1
e1 = 0
e2 = 0
t = 1
imagI = complex(0, 1)
k_list = [0, pi]

# G_k_w


def G_k_w(k, w):
    z = w + complex(0, eta)

    G11 = pow((z - e1) - pow(t, 2) / (z - e2), -1)
    G12 = - (t / (z - e2)) * pow((z - e1) - pow(t, 2) / (z - e2), -1)

    G22 = pow((z - e2) - pow(t, 2) / (z - e1), -1)
    G21 = - (t / (z - e1)) * pow((z - e2) - pow(t, 2) / (z - e1), -1)

    return exp(imagI * k * 0) * G11 + exp(imagI * k * 1) * G12 + \
        exp(imagI * k * 0) * G22 + exp(imagI * k * 1) * G21


def A_k_w(k, w):
    return -(1 / pi) * imag(G_k_w(k, w))


def DOS(w):
    s = 0
    for k in k_list:
        s += A_k_w(k, w) / len(k_list)
    return s


for k in k_list:
    A_list = []
    w_list = []
    w = -6
    while(w <= 6):
        A_list.append(A_k_w(k, w))
        w_list.append(w)
        w += 0.01
    plt.plot(w_list, A_list)

plt.show()

DOS_list = []
w_list = []
w = -6
while(w <= 6):
    DOS_list.append(DOS(w))
    w_list.append(w)
    w += 0.01

plt.plot(w_list, DOS_list)
plt.show()
