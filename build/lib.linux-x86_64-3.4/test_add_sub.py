from mym import *

A = Matrix([[1, 2], 
            [3, 4]
           ])

B = Matrix([[5, 6], 
            [7, 8]
           ])

print("A + B : ")
print(add(A, B))

print("A - B : ")
print(sub(A, B))

print("A - 5 : ")
print(sub(A, mul_cons(eye(A.size()[0]), 5)))