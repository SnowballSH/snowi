from collections import Counter

from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("geometry/closest-pair")


def rnd_pt(c):
  return (rnd(-c, c), rnd(-c, c))


def d2(a, b):
  return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2


def brute(v):
  return min(d2(v[i], v[j]) for i in range(len(v)) for j in range(i))


def check_pair(v, want):
  """Two distinct input elements at the minimum distance."""
  a, b = lib.closest_pair(v)
  check_eq(d2(a, b), want, str(v))
  count = Counter(v)
  check(count[a] >= 1 + (a == b) and count[b] >= 1, str(v))


for it in range(iters(3000)):
  n = rnd(2, 12)
  C = 10**9 if it % 4 == 0 else rnd(1, 6)
  v = [rnd_pt(C) for _ in range(n)]
  if it % 5 == 1:  # all on a vertical or horizontal line
    v = [(0, y) if it % 2 else (x, 0) for x, y in v]
  if it % 7 == 2:
    v[-1] = v[0]
  check_pair(v, brute(v))
for _ in range(iters(200)):  # larger inputs, tiny grid: many ties
  C = rnd(10, 60)
  v = [rnd_pt(C) for _ in range(rnd(20, 200))]
  check_pair(v, brute(v))

C = 10**9
check_pair([(-C, -C), (C, C)], 8 * C * C)
check_pair([(3, 4), (3, 4)], 0)
check_pair([(-C, -C), (C, C), (C, -C)], 4 * C * C)
# Extra fields ride along: indices of the pair.
a, b = lib.closest_pair([(0, 0, 0), (9, 9, 1), (1, 2, 2), (5, 5, 3)])
check_eq(sorted([a[2], b[2]]), [0, 2])
# Floats work unchanged.
a, b = lib.closest_pair([(0.5, 0.25), (3.0, 1.0), (0.75, 0.0), (2.5, 2.5)])
check(abs(d2(a, b) - 0.125) < 1e-12)

n = 200_000
v = [rnd_pt(C) for _ in range(n)]
v[-1] = (v[0][0] + 1, v[0][1])
with Timer("closest pair of 2e5 random points"):
  a, b = lib.closest_pair(v)
check_eq(d2(a, b), 1)
v = [(7, 3 * i - C) for i in range(n)]
v[n // 2] = (7, v[n // 2][1] + 1)
with Timer("closest pair of 2e5 points on a vertical line"):
  a, b = lib.closest_pair(v)
check_eq(d2(a, b), 4)
v = [(i % 3, i // 3) for i in range(n)]
with Timer("closest pair of a 3 x 66667 grid"):
  a, b = lib.closest_pair(v)
check_eq(d2(a, b), 1)
print("OK geometry/closest-pair")
