import sys
sys.path.append("..")

from mym import *

t = 1.0/3

H1 = Matrix([[-t,  -2*t, -2*t], 
             [-2*t, 2*t,  -t], 
             [-2*t, -t,  2*t]])

A = Matrix([[1, 1, 1], 
            [2, -1, -1], 
            [2, -4, 5]])

A1 = mul(H1, A)

print(A1)

H2 = Matrix([[1, 0, 0], 
             [0, 0, 1], 
             [0, 1, 0]])

A2 = mul(H2, A1)

print(A2)

R = A2
Q = mul(H1, H2)

print(Q)
print(mul(Q, R))