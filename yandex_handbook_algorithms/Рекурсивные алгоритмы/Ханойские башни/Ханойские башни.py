def HanoiTowers(n, steps, fromPeg, toPeg):
    if n == 1:
        steps.append((fromPeg, toPeg))
        return
    unusedPeg = 6 - fromPeg - toPeg
    HanoiTowers(n - 1, steps, fromPeg, unusedPeg)
    steps.append((fromPeg, toPeg))
    HanoiTowers(n - 1, steps, unusedPeg, toPeg)

n = int(input())
steps = []
HanoiTowers(n, steps, 1, 3)

print(len(steps))
for step in steps:
    print(str(step[0]) + ' ' + str(step[1]))
