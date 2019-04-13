###### Tue Mar 26 17:11:59 IST 2019

import numpy as np

n = int(input("Enter row number:"))

def bin2dec(n):
    p = []
    if (n == 0):
        p.append(0)
    while (n>0):
        p.append(n%2)
        n = int(n/2)
    
    p.reverse()
    return p

for i in range(1,n+1):
    print(i, bin2dec(i-1))