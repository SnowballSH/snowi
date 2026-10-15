import math

from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("geometry/convex-hull")


def rnd_pt(c):
  return (rnd(-c, c), rnd(-c, c))


def cr(o, a, b):
  return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def brute_on(a, b, p):
  return (
    cr(a, b, p) == 0
    and min(a[0], b[0]) <= p[0] <= max(a[0], b[0])
    and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])
  )


def in_triangle(a, b, c, p):
  if cr(a, b, c) == 0:  # degenerate: the segments cover it
    return False
  x, y, z = cr(a, b, p), cr(b, c, p), cr(c, a, p)
  return (x >= 0 and y >= 0 and z >= 0) or (x <= 0 and y <= 0 and z <= 0)


def brute_vertices(pts):
  """Points not in a closed triangle or segment of the others."""
  pts = sorted(set(pts))
  n, res = len(pts), set()
  for i in range(n):
    others = pts[:i] + pts[i + 1 :]
    m, p = n - 1, pts[i]
    covered = any(
      brute_on(others[a], others[b], p)
      for a in range(m)
      for b in range(a + 1, m)
    ) or any(
      in_triangle(others[a], others[b], others[c], p)
      for a in range(m)
      for b in range(a + 1, m)
      for c in range(b + 1, m)
    )
    if not covered:
      res.add(p)
  return res


def brute_in_hull(h, p):
  n = len(h)
  if n <= 2:
    return 0 if n and brute_on(h[0], h[-1], p) else -1
  res = 1
  for i in range(n):
    c = cr(h[i - 1], h[i], p)
    if c < 0:
      return -1
    if c == 0:
      res = 0
  return res


def check_hull(pts):
  h = lib.convex_hull(pts)
  want = brute_vertices(pts)
  check_eq(len(h), len(want), str(pts))
  check_eq(set(h), want)
  check_eq(h[0], min(want))
  n = len(h)
  if n >= 3:
    for i in range(n):
      check(cr(h[i - 2], h[i - 1], h[i]) > 0)
      if 1 <= i < n - 1:
        check(cr(h[0], h[i], h[i + 1]) > 0)
  best = max(lib.dist2(a, b) for a in pts for b in pts)
  a, b = lib.hull_diameter(h)
  check_eq(lib.dist2(a, b), best)
  check(a in want and b in want)
  return h


for it in range(iters(1200)):
  n = rnd(1, 9)
  C = rnd(1, 4) if it % 3 else 10**9
  pts = [rnd_pt(C) for _ in range(n)]
  if it % 4 == 0:  # collinear runs
    k = [rnd(-3, 3) for _ in range(n)]
    pts = [(pts[0][0] + t, pts[0][1] + 2 * t) for t in k]
  if it % 5 == 0:
    pts.append(pts[0])
  h = check_hull(pts)
  if C < 10:
    for x in range(-C - 1, C + 2):
      for y in range(-C - 1, C + 2):
        check_eq(lib.in_hull(h, (x, y)), brute_in_hull(h, (x, y)))
  for q in range(20):
    p = rnd_pt(C)
    if q % 3 == 0:
      p = h[rnd(0, len(h) - 1)]
    if q % 3 == 1:
      a, b = h[rnd(0, len(h) - 1)], h[rnd(0, len(h) - 1)]
      p = ((a[0] + b[0]) // 2, (a[1] + b[1]) // 2)
    check_eq(lib.in_hull(h, p), brute_in_hull(h, p))

check_eq(lib.convex_hull([]), [])
check_eq(lib.convex_hull([(2, 2), (2, 2)]), [(2, 2)])
sq = [(0, 0), (2, 0), (2, 2), (0, 2), (1, 1), (1, 0), (2, 1)]
h = lib.convex_hull(sq)
check_eq(h, [(0, 0), (2, 0), (2, 2), (0, 2)])
check_eq(
  [lib.in_hull(h, p) for p in ((1, 1), (0, 0), (0, 1), (3, 1))], [1, 0, 0, -1]
)
check_eq(lib.in_hull([], (0, 0)), -1)
# Sorted by x only (the sort by x, then y, is still needed).
xs = [(0, 2), (0, 3), (0, 0), (1, 0), (2, 3), (3, 0)]
check_eq(lib.convex_hull(xs), [(0, 0), (3, 0), (2, 3), (0, 3)])
seg = lib.convex_hull([(0, 0), (3, 3), (1, 1)])
check_eq(seg, [(0, 0), (3, 3)])
check_eq((lib.in_hull(seg, (2, 2)), lib.in_hull(seg, (4, 4))), (0, -1))

C, n = 10**9, 200_000
pts = [rnd_pt(C) for _ in range(n)]
pts[:4] = [(-C, -C), (C, -C), (-C, C), (C, C)]
with Timer("hull of 2e5 points + 2e5 in_hull + diameter"):
  h = lib.convex_hull(pts)
  ok = all(lib.in_hull(h, p) >= 0 for p in pts)
  a, b = lib.hull_diameter(h)
check(ok and h[0] == (-C, -C) and len(h) == 4)
check_eq(lib.dist2(a, b), 8 * C * C)
circle = [
  (round(C * math.cos(t)), round(C * math.sin(t)))
  for t in (2 * math.pi * i / n for i in range(n))
]
with Timer("hull of 2e5 points in convex position + diameter"):
  h = lib.convex_hull(circle)
  a, b = lib.hull_diameter(h)
check(len(h) > n // 2 and lib.dist2(a, b) >= 4 * C * C - 4 * C)
check(all(cr(h[i - 2], h[i - 1], h[i]) > 0 for i in range(len(h))))
print("OK geometry/convex-hull")
