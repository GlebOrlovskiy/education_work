n, k = map(int, input().split())


def fact(n):
    if n == 0:
        return 1
    result = 1
    for i in range(2, n + 1):
        result = result*i
    return result

f_k = fact(k)
f_n = fact(n)
f_n_minus_k = fact(n - k)

print(int(f_n/(f_k*f_n_minus_k)))
