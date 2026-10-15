import math
from math import gcd

from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("geometry/minkowski", "geometry/convex-hull")


def rnd_pt(c):
  return (rnd(-c, c), rnd(-c, c))


def cr(o, a, b):
  return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def on_seg(a, b, p):
  return (
    cr(a, b, p) == 0
    and min(a[0], b[0]) <= p[0] <= max(a[0], b[0])
    and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])
  )


def check_sum(a, b):
  """The sum is hull(S), S = {a + b}: vertices in S, strictly convex,
  counter-clockwise from its smallest point, containing all of S."""
  s = lib.minkowski(a, b)
  S = {(p[0] + q[0], p[1] + q[1]) for p in a for q in b}
  n = len(s)
  check(1 <= n <= len(a) + len(b), str((a, b)))
  check(all(p in S for p in s))
  check_eq(s[0], min(S))
  if n == 1:
    check_eq(len(S), 1)
  if n == 2:
    check(s[0] != s[1] and all(on_seg(s[0], s[1], p) for p in S))
  if n >= 3:
    check(all(cr(s[i - 2], s[i - 1], s[i]) > 0 for i in range(n)))
    check(all(cr(s[0], s[i], s[i + 1]) > 0 for i in range(1, n - 1)))
    check(all(cr(s[i - 1], s[i], p) >= 0 for p in S for i in range(n)))
  return s


def check_support(a, b, s, d):
  """The extreme value in direction d adds up."""

  def best(v):
    return max(p[0] * d[0] + p[1] * d[1] for p in v)

  check_eq(best(s), best(a) + best(b))


def rnd_convex(n, c):
  return lib.convex_hull([rnd_pt(c) for _ in range(n)])


for it in range(iters(2000)):
  C = 999_999_997 if it % 4 == 0 else rnd(1, 6)
  a = rnd_convex(rnd(1, 12), C)
  b = rnd_convex(rnd(1, 12), C)
  if it % 5 == 1:  # a parallel copy: every edge direction shared
    k = rnd(1, 2) if C < 10 else 1
    b = lib.convex_hull(
      [(k * x + rnd(-3, 3), k * y + rnd(-3, 3)) for x, y in a]
    )
  r = rnd(0, len(a) - 1)
  a = a[r:] + a[:r]
  s = check_sum(a, b)
  for _ in range(10):
    check_support(a, b, s, rnd_pt(10**9))

tri, sq = [(0, 0), (2, 0), (0, 1)], [(0, 0), (1, 0), (1, 1), (0, 1)]
check_eq(lib.minkowski(tri, sq), [(0, 0), (3, 0), (3, 1), (1, 2), (0, 2)])
check_eq(lib.minkowski([(1, 2)], [(3, 4)]), [(4, 6)])
check_eq(lib.minkowski([(1, 1)], tri), [(1, 1), (3, 1), (1, 2)])
s1, s2, s3 = [(0, 0), (2, 1)], [(5, 3), (1, 1)], [(0, 3), (0, 0)]
check_eq(lib.minkowski(s1, s2), [(1, 1), (7, 4)])
check_eq(lib.minkowski(s1, s3), [(0, 0), (2, 1), (2, 4), (0, 3)])
C = 10**9
big = [(-C, -C), (C, -C), (C, C), (-C, C)]
check_eq(lib.minkowski(big, big), [(2 * x, 2 * y) for x, y in big])

# Collision of A and B + t is t in A - B; distance between disjoint
# A and B + t is the distance from t to A - B.
for it in range(iters(300)):
  D = 5 * 10**8 if it % 2 else 10
  A, B = rnd_convex(rnd(1, 8), D), rnd_convex(rnd(1, 8), D)
  t = rnd_pt(D)
  d = lib.minkowski(A, [(-x, -y) for x, y in B])
  meet = lib.in_hull(d, t) >= 0
  Bt = [(x + t[0], y + t[1]) for x, y in B]
  want = any(lib.in_hull(Bt, p) >= 0 for p in A) or any(
    lib.in_hull(A, p) >= 0 for p in Bt
  )
  for i in range(len(A)):
    for j in range(len(Bt)):
      want = want or lib.seg_inter(A[i - 1], A[i], Bt[j - 1], Bt[j]) != []
  check_eq(meet, want)
  if meet or it % 2:
    continue
  f = [(x - t[0], y - t[1]) for x, y in d]
  dist = min(lib.seg_dist(f[i - 1], f[i], (0, 0)) for i in range(len(f)))
  brute = min(
    lib.seg_dist(X[i - 1], X[i], p)
    for X, Y in ((A, Bt), (Bt, A))
    for i in range(len(X))
    for p in Y
  )
  check(abs(dist - brute) < 1e-9 * max(1.0, brute))

# Large: 2e5 + 2e5 vertices near circles of radius 1e9.
n = 200_000
pa = lib.convex_hull(
  [
    (
      round(C * math.cos(2 * math.pi * i / n)),
      round(C * math.sin(2 * math.pi * i / n)),
    )
    for i in range(n)
  ]
)
pb = lib.convex_hull(
  [
    (
      round(C / 2 * math.cos(2 * math.pi * (i + 0.5 * (i % 2)) / n)),
      round(C / 2 * math.sin(2 * math.pi * (i + 0.5 * (i % 2)) / n)),
    )
    for i in range(n)
  ]
)
with Timer(f"sum of {len(pa)} + {len(pb)} vertices"):
  s = lib.minkowski(pa, pb)
dirs = set()
for v in (pa, pb):
  for i in range(len(v)):
    ex, ey = v[i][0] - v[i - 1][0], v[i][1] - v[i - 1][1]
    g = gcd(ex, ey)
    dirs.add((ex // g, ey // g))
check_eq(len(s), len(dirs))
check(all(cr(s[i - 2], s[i - 1], s[i]) > 0 for i in range(len(s))))
for _ in range(50):
  check_support(pa, pb, s, rnd_pt(C))
print("OK geometry/minkowski")
