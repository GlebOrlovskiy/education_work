def fast_rocks(n, m):
    if n % 2 == 0 and m % 2 == 0:
        return False
    return True

n, m = map(int, input().split(' '))

if fast_rocks(n, m):
    print('Win')
else:
    print('Lose')
