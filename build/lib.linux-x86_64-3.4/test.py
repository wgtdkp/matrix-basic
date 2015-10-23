# --*-- utf-8 --*--
from mym import *

#m = Matrix([[1, 2], [3, 4]])
#print(m)
#n = Matrix([[1, 2], [3, 4]])
#print(mul(m, n))
#print(det(m))
a = Matrix([[1, 2, 3, 4, 5], 
		    [11, 23, 3, -9, 10],
            [19, 5, 9, 23, 0],
            [7, 43, 90, 21, 33],
            [0, 1, 4, 8, 6]
           ])
b = Matrix([[1],
	        [2],
            [3],
            [4],
            [5]
           ])

print("the A matrix is: ")
print(a)
print("gauss_elim the euqation's solution is: ")
print(gauss_elim(A=a.copy(), B=b.copy(), step=False))

"""
print("the A matrix is: ")
print(a)
print("me_gauss_elim the euqation's solution is: ")
print(me_gauss_elim(A=a.copy(), B=b.copy(), step=False))

print("the A matrix is: ")
print(a)
print("tri_decomp the euqation's solution is: ")
print(tri_decomp(A=a.copy(), B=b.copy(), step=False))

print("the A matrix is: ")
print(a)
print("me_tri_decomp the euqation's solution is: ")
print(me_tri_decomp(A=a.copy(), B=b.copy(), step=False))
"""

"""
L = [[1,  11, 19, 7,  0],
     [11, 23, 5,  43, 1],
     [19, 5,  9,  90, 4],
     [7,  43, 90, 21, 8],
     [0,  1,  4,  8,  6],
    ]
"""

#symmetric positive definite matrix
a = Matrix([[23, 2, 3], 
     [2, 5, 3], 
     [3, 3, 11]
    ])
b = Matrix([[1], [2], [3]])
print("the A matrix is: ")
print(a)
print("cholesky_decomp the euqation's solution is: ")
#print(cholesky_decomp(A=a.copy(), B=b.copy()))
print(en_cholesky_decomp(A=a.copy(), B=b.copy()))

a = Matrix([[10, 2, 0],
            [3,  13, 7 ],
            [0,  -3, 8],
           ])
b = Matrix([[1], [2], [3]])
print("the A matrix is: ")
print(a)
print("chasing_method the euqation's solution is: ")
print(chasing_method(A=a.copy(), B=b.copy()))


a = Matrix([[10, 2, 0],
            [3,  13, 7 ],
            [0,  -3, 8],
           ])
b = Matrix([[1], [2], [3]])
print("the A matrix is: ")
print(a)
print("jacobi the euqation's solution is: ")
print(jacobi(A=a.copy(), B=b.copy(), delta=1e-6))

a = Matrix([[10, 2, 0],
            [3,  13, 7 ],
            [0,  -3, 8],
           ])
b = Matrix([[1], [2], [3]])
print("the A matrix is: ")
print(a)
print("Gauss-Seidel the euqation's solution is: ")
print(gauss_seidel(A=a.copy(), B=b.copy(), delta=1e-6))

a = Matrix([[10, 2, 0],
            [3,  13, 7 ],
            [0,  -3, 8],
           ])
b = Matrix([[1], [2], [3]])
print("the A matrix is: ")
print(a)
print("SOR the euqation's solution is: ")
print(sor(A=a.copy(), B=b.copy(), w=1.1, delta=1e-6))

print("the A matrix's norm: ")
print(norm(a, 1))
print(norm(a, 2))
print(norm(a, (1<<31) - 1))