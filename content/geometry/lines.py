"""
Title: Lines and segments
Description: `side_of(s, e, p)` is 1 if `p` is left of the
  directed line `s -> e`, 0 on it and -1 right of it.
  `on_segment` tests the closed segment. `line_inter` returns
  `(1, point)`, or `(0, None)` for parallel lines and `(-1, None)`
  for the same line. `seg_inter` intersects two closed segments: a
  list of no point, one point, or the two ends of a collinear
  overlap, sorted. `line_dist` is the signed distance from `p` to
  line `ab` (positive on the left), `seg_dist` the distance to
  segment `ab`. `line_proj` projects `p` onto line `ab`, or
  reflects it across the line when `refl` is true.
Precondition: The two points defining a line differ. With integer
  input every decision is exact; computed points and distances are
  floats (the points of `line_inter` and `seg_inter` are the exact
  ones, correctly rounded).
Time: $O(1)$.
Usage:
  side_of(a, b, p), on_segment(a, b, p)
  kind, x = line_inter(a, b, c, d)
  v = seg_inter(a, b, c, d)  # len(v) <= 2
  seg_dist(a, b, p), line_dist(a, b, p)
  q, r = line_proj(a, b, p), line_proj(a, b, p, True)
Tweak: Exact computed points: pass `Fraction` coordinates (or
  replace `/` by `Fraction(num, den)`). For non-integer floats, the
  `== 0` in `on_segment` should become a tolerance:
  `seg_dist(s, e, p) < eps`.
Example:
  >>> seg_inter((0, 0), (2, 2), (2, 0), (0, 2))
  [(1.0, 1.0)]
  >>> seg_inter((0, 0), (4, 0), (3, 0), (1, 0))
  [(1, 0), (3, 0)]
  >>> line_inter((0, 0), (1, 0), (0, 1), (1, 1))
  (0, None)
  >>> line_proj((0, 0), (1, 1), (2, 0), True)
  (0.0, 2.0)
Notes: Two segments cross properly iff the ends of each lie
  strictly on opposite sides of the other (`side_of` products are
  -1); everything else is a touching or collinear case, which
  `on_segment` decides. Read whether touching counts. "Mirror",
  "reflection", "unfold the path" call for `line_proj(.., True)`;
  "distance to a wall or road" is `seg_dist`. Float coordinates
  near $10^9$ round the cross products, so computed points may be
  off by about $10^(-5)$; integer input avoids that.
Uses: geometry/point
Source: kactl SegmentIntersection, LineIntersection,
  SegmentDistance, LineProjectionReflection (CC0).
Status: stress-tested against brute force
"""

import math


def side_of(s, e, p):
  c = cross3(s, e, p)
  return (c > 0) - (c < 0)


def on_segment(s, e, p):
  return cross3(p, s, e) == 0 and dot(sub(s, p), sub(e, p)) <= 0


def line_inter(s1, e1, s2, e2):
  d = cross(sub(e1, s1), sub(e2, s2))
  if d == 0:
    return (-1 if cross3(s1, e1, s2) == 0 else 0), None
  q = cross3(s2, e2, s1)
  x = (s1[0] * d + (e1[0] - s1[0]) * q) / d
  return 1, (x, (s1[1] * d + (e1[1] - s1[1]) * q) / d)


def seg_inter(a, b, c, d):
  oa, ob = cross3(c, d, a), cross3(c, d, b)
  oc, od = cross3(a, b, c), cross3(a, b, d)
  if oa * ob < 0 and oc * od < 0:
    t = ob - oa
    x = (a[0] * ob - b[0] * oa) / t
    return [(x, (a[1] * ob - b[1] * oa) / t)]
  s = set()
  for p, u, v in ((a, c, d), (b, c, d), (c, a, b), (d, a, b)):
    if on_segment(u, v, p):
      s.add(p)
  return sorted(s)


def line_dist(a, b, p):
  return cross3(a, b, p) / math.hypot(b[0] - a[0], b[1] - a[1])


def seg_dist(a, b, p):
  if dot(sub(p, a), sub(b, a)) <= 0:
    return math.hypot(p[0] - a[0], p[1] - a[1])
  if dot(sub(p, b), sub(a, b)) <= 0:
    return math.hypot(p[0] - b[0], p[1] - b[1])
  return abs(line_dist(a, b, p))


def line_proj(a, b, p, refl=False):
  vx, vy = b[0] - a[0], b[1] - a[1]
  t = (1 + refl) * cross3(a, b, p) / (vx * vx + vy * vy)
  return (p[0] + vy * t, p[1] - vx * t)
