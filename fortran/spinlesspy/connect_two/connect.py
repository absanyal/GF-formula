# Fri Mar 15 20:58:15 IST 2019

# A sample to form the 12 site 11 particle problem from the
# 6s5p and 6s6p problems for U=0.0 at w=0.0

# Use np.kron(a, b) for outer product

import numpy as np
import os
import statemanip as sm
from numpy.linalg import inv
# from numpy import imag, trace
import matplotlib.pyplot as plt
import time

os.system('rm *.dat')
os.system('clear')

