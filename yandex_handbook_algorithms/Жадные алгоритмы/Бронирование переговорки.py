n = int(input())

segments = []

for i in range(n):
    l, r = map(int, input().split())
    segments.append((l, r))

segments = sorted(segments, key=lambda seg: seg[1])
new_segments = []

def search(segments, leader = 0):
    if leader >= len(segments) - 1:
        return
      
    r_ = segments[leader][1]
    old_segments = segments[:]

    for i in range(leader + 1, len(old_segments)):
        if old_segments[i][0] <= r_:
            segments.remove(old_segments[i])

    search(segments, leader + 1)

search(segments)
print(len(segments))
