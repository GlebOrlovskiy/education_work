n, k = map(int, input().split())

def fact(n):
    if n == 0:
        return 1
    result = 1
    for i in range(1, n + 1):
        result = result*i
    return result

def C(n, k):
    f_k = fact(k)
    f_n = fact(n)
    f_n_minus_k = fact(n - k)
    return int(f_n/(f_k*f_n_minus_k))

print(C(n + k - 1, k))
