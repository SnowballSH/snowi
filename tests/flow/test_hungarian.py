import random
from itertools import permutations

from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("flow/hungarian")
hungarian = lib.hungarian


def brute(a, n, m):
  """Minimum over all injections rows -> columns."""
  return min(
    sum(a[i][c[i]] for i in range(n)) for c in permutations(range(m), n)
  )


def check_match(a, n, m, cost, match):
  check_eq(len(match), n)
  check(all(0 <= j < m for j in match))
  check_eq(len(set(match)), n)
  check_eq(sum(a[i][match[i]] for i in range(n)), cost)


def planted(n, m):
  """Duals u, v (v <= 0, 0 on free columns) and a tight injection."""
  perm = list(range(m))
  random.Random(rnd(0, 10**9)).shuffle(perm)
  perm = perm[:n]
  u = [rnd(-(10**6), 10**6) for _ in range(n)]
  v = [0] * m
  for j in perm:
    v[j] = rnd(-(10**6), 0)
  a = [[u[i] + v[j] + rnd(1, 10**6) for j in range(m)] for i in range(n)]
  for i in range(n):
    a[i][perm[i]] = u[i] + v[perm[i]]
  return a, perm, sum(u) + sum(v)


for _ in range(iters(2000)):
  m = rnd(0, 6)
  n = rnd(0, m)
  kind = rnd(0, 2)
  A = [5, 10**9, (1 << 60) - 1][kind]
  a = [
    [
      A * rnd(-1, 1) if kind == 2 and rnd(0, 1) else rnd(-A, A)
      for _ in range(m)
    ]
    for _ in range(n)
  ]
  cost, match = hungarian(a)
  check_eq(cost, brute(a, n, m))
  check_match(a, n, m, cost, match)

for _ in range(iters(1000)):  # Tweak: float costs work as they are
  m = rnd(0, 6)
  n = rnd(0, m)
  a = [[rnd(-(10**6), 10**6) / 997 for _ in range(m)] for _ in range(n)]
  cost, match = hungarian(a)
  check(abs(cost - brute(a, n, m)) < 1e-6)
  check_eq(len(set(match)), n)
  check(abs(sum(a[i][match[i]] for i in range(n)) - cost) < 1e-6)

for _ in range(iters(200)):
  m = rnd(1, 40)
  n = rnd(1, m)
  a, perm, opt = planted(n, m)
  check_eq(hungarian(a), (opt, perm))

a, perm, opt = planted(200, 300)
with Timer("planted 200 x 300"):
  res = hungarian(a)
check_eq(res, (opt, perm))

n = 300  # a[i][j] = -i j: identity is optimal (rearrangement)
a = [[-i * j for j in range(n)] for i in range(n)]
with Timer("300 x 300, a[i][j] = -i j (many Dijkstra steps)"):
  cost, match = hungarian(a)
check_eq(cost, -sum(i * i for i in range(n)))
check_match(a, n, n, cost, match)
print("OK flow/hungarian")
