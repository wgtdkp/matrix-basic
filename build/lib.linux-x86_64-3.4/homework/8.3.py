
A = [[7, 3, -2], [3, 4, -1], [-2, -1, 3]]
V = [1, 1, 1]
U = [1, 1, 1]
while True:
    for i in range(0, 3):
        sum = 0
        for j in range(0, 3):
            sum += A[i][j] * V[j]
        U[i] = sum
    print(U)
    max = 0
    for i in range(0, 3):
        if abs(U[i]) > max:
            max = abs(U[i])
    for i in range(0, 3):
        U[i] /= max
    print(U)
    print("\n")
    cnt = 0
    for i in range(0, 3):
        if abs(U[i] - V[i]) < 1e-3:
            cnt = cnt + 1
    if cnt == 3:
        break
    V = U[:]
