import sys
sys.path.append("..")

from mym import *

U = Matrix([[1,    0,    0], 
            [0, -0.6, -0.8], 
            [0, -0.8,  0.6]])

A = Matrix([[1, 3, 4], 
            [3, 1, 2], 
            [4, 2, 1]])

RES = mul(U, mul(A, U))
print(RES)