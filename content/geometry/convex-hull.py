"""
Title: Convex hull
Description: `convex_hull` returns the vertices of the convex hull
  counter-clockwise, starting from the smallest point (by `x`, then
  `y`), without collinear points. `in_hull` returns 1 if `p` is
  strictly inside such a hull, 0 if on its boundary and -1 if
  outside. `hull_diameter` returns the farthest pair of points.
Precondition: `in_hull` and `hull_diameter` take a hull as returned
  by `convex_hull` (non-empty for `hull_diameter`). Exact for
  integers.
Time: $O(n log n)$ for the hull, $O(log n)$ per `in_hull` query,
  $O(n)$ for `hull_diameter`.
Usage:
  h = convex_hull(pts)  # 1 or 2 points if degenerate
  where = in_hull(h, p)  # 1 in, 0 edge, -1 out
  a, b = hull_diameter(h); d2 = dist2(a, b)
Tweak: To keep points in the middle of hull edges, pop on `< 0`
  instead of `<= 0`; if all points are collinear the middle ones
  then appear twice, so return the sorted `pts` in that case.
  `in_hull` and `hull_diameter` need the hull without them.
Example:
  >>> pts = [(0, 0), (1, 2), (2, 1), (3, 3), (4, 0), (5, 2)]
  >>> h = convex_hull(pts); h
  [(0, 0), (4, 0), (5, 2), (3, 3), (1, 2)]
  >>> in_hull(h, (2, 1)), in_hull(h, (2, 0)), in_hull(h, (5, 0))
  (1, 0, -1)
  >>> a, b = hull_diameter(h); a, b, dist2(a, b)
  ((5, 2), (0, 0), 29)
  >>> d = [math.sqrt(dist2(h[i - 1], h[i])) for i in range(5)]
  >>> round(sum(d), 3)  # perimeter
  12.944
Notes: The answer only ever involves hull vertices when the
  objective is a maximum of a linear function (max of $a x + b y$
  over points), a farthest pair, a smallest enclosing shape or
  "mixtures" (any weighted average of options). If every element
  gets $x$ added, each candidate (a subarray, a window) is a line
  in $x$ with its length as slope, and the upper envelope of these
  lines is the upper hull of the points (slope, intercept).
  Decide whether collinear points belong to the hull and whether
  touching counts before coding. Points that arrive sorted by `x`,
  then `y`, without duplicates, can skip the sort (each chain is
  $O(n)$); sorted by `x` alone is not enough.
Uses: geometry/lines
Source: kactl ConvexHull, PointInsideHull, HullDiameter (CC0);
  UCup S4-S5 Training Packet 20.
Status: stress-tested against brute force
"""


def convex_hull(pts):
  pts = sorted(set(pts))
  if len(pts) <= 2:
    return pts
  h = []
  for chain in (pts, pts[::-1]):
    s = len(h)
    for p in chain:
      while len(h) >= s + 2 and cross3(h[-2], h[-1], p) <= 0:
        h.pop()
      h.append(p)
    h.pop()
  return h


def in_hull(h, p):
  n = len(h)
  if n < 3:
    return 0 if n and on_segment(h[0], h[-1], p) else -1
  a, b = cross3(h[0], h[1], p), cross3(h[0], h[-1], p)
  if a < 0 or b > 0:
    return -1
  lo, hi = 1, n - 1
  while hi - lo > 1:
    mid = (lo + hi) // 2
    if cross3(h[0], h[mid], p) >= 0:
      lo = mid
    else:
      hi = mid
  c = cross3(h[lo], h[hi], p)
  if c < 0:
    return -1
  return 0 if c == 0 or a == 0 or b == 0 else 1


def hull_diameter(h):
  n, j = len(h), 1 % len(h)
  best, res = -1, None
  for i in range(n):
    a, b = h[i], h[(i + 1) % n]
    ex, ey = b[0] - a[0], b[1] - a[1]
    while True:
      c, d = h[j], h[(j + 1) % n]
      if ex * (d[1] - c[1]) - ey * (d[0] - c[0]) <= 0:
        break
      j = (j + 1) % n
    for p in (a, b):
      d2 = dist2(p, h[j])
      if d2 > best:
        best, res = d2, (p, h[j])
  return res
