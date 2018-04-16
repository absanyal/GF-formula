# -*- coding: utf-8 -*-
"""
Created on Sun Sep 24 20:17:26 2017

@author: AB Sanyal

A module form of the progressbar that can be reused anywhere.
"""

# imports
import sys
from itertools import cycle

# symbols = [
#         chr(8942),
#         chr(8944),
#         chr(8943),
#         chr(8945)
#         ]

symbols = [
    '(  ',
    ' | ',
    '  )',
    ' | '
]

sym_list = cycle(symbols)

progress_symbol = chr(9608)
emptysymbol = chr(9617)

stepsize = 5


def progressbar(v, v_min, v_max):

    # percentage completed
    pc = round((v - v_min) / (v_max - v_min) * 100, 1)

    counter = int(pc / stepsize)

    symbol = next(sym_list)

    message = "Progress: " + symbol + "  " + progress_symbol * counter \
        + emptysymbol * int((100 / stepsize) - counter) + " " + str(pc) + " %"
    sys.stdout.write('\r' + message)

    if (v >= v_max):
        done = "Progress: Done!"
        sys.stdout.write('\r' + " " * len(message))
        sys.stdout.write('\r' + done)
