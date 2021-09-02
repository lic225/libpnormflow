def parallel(n):
    m = 2 * (n - 2)
    E = []

    for i in range(1, n-1):
        E.append((0, i))
        E.append((i, n-1))
    return (n, m, E)

def path(n):
    m = n-1
    E = []
    for i in range(n-1):
        E.append((i, i+1))
    return (n, m, E)

def grid(k): # k*k grid
    n = k*k

    def id(i, j):
        return i*k + j
    E = []
    for i in range(k):
        for j in range(k):
            if j + 1 < k:
                E.append((id(i, j), id(i, j+1)))
            if i + 1 < k:
                E.append((id(i, j), id(i+1, j)))
    return (n, len(E), E)

def l2parallel(k):
    n = 2 + k * (k - 1)
    E = [(0, n - 1)]
    for j in range(1, k - 1):
        for i in range(k):
            a = (j - 1) * k + i + 1
            b = j * k + i + 1
            E.append((a, b))

    for i in range(k):
        E.append((0, i + 1))
        E.append((n - i - 1, n - 1))

    return (n, len(E), E)

# (n, m, E) = parallel(100000)
# (n, m, E) = grid(10)
# (n, m, E) = path(1000)
(n, m, E) = l2parallel(100)
print(n, n, m)
for (u, v) in E:
    print(u, v)
