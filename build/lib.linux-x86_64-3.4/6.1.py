from mym import *

"""
A = Matrix([[5, 2, 1], [-1, 4, 2], [2, -3, 10]])
B = Matrix([[-12], [20], [3]])

print("solve by jacobi:")
print(jacobi(A, B, delta=1e-4, step=True))

print("\n")
print(gauss_seidel(A, B, delta=1e-4, step=True))

"""

"""
A = Matrix([[5, 2, 1], [-1, 4, 2], [2, -3, 10]])
D = Matrix([[5, 0, 0], [-0.9, 4, 0], [1.8, -2.7, 10]])

#print(A)
#print(D)
D_ = inv(D.copy())
#print(D_)
W = Matrix([[0.9, 0, 0], [0, 0.9, 0], [0, 0, 0.9]])
print(mul(W, mul(D_, A)))
#print(A)
"""


A = Matrix([[4, -1, 0], [-1, 4, -1], [0, -1, 4]])
B = Matrix([[1], [4], [-3]])

print(sor(A, B, 1.03, 5e-6, True))
print(sor(A, B, 1, 5e-6, True))
print(sor(A, B, 1.1, 5e-6, True))