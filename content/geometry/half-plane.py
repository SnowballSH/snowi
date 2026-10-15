"""
Title: Half-plane intersection
Description: Intersects the half-planes left of the directed lines
  `(s, e)` inside the box $[-"big", "big"]^2$ and returns the
  vertices of the result counter-clockwise. Sorts the lines by
  angle and keeps the useful ones in a deque.
Precondition: Points are tuples, `s != e` in each half-plane;
  results are floats. Points within `eps` $= 10^(-9)$ of a line,
  plus $10^(-15) (|x| + |y|)$ for rounding far out, count as on
  it: the result is empty when the intersection is empty or has no
  interior (a point or a segment), and has no repeated or collinear
  vertices. Directions closer than `eps` count as parallel. Fine
  for coordinates up to about $10^4$.
Time: $O(n log n)$; about a second for $2 dot 10^5$ half-planes.
Usage:
  hs = [(s, e), ...]  # keep left of s -> e
  poly = half_plane(hs)  # [] if empty
Tweak: Scale `eps` with larger coordinates. The default box only
  keeps unbounded results finite; pass a smaller `big` when the
  problem bounds the plane.
Example:
  >>> sq = [(0, 0), (4, 0), (4, 4), (0, 4)]
  >>> hs = [(sq[i - 1], sq[i]) for i in range(4)]
  >>> half_plane(hs + [((0, 1), (4, 3))])  # also above a line
  [(0.0, 1.0), (4.0, 3.0), (4.0, 4.0), (0.0, 4.0)]
  >>> half_plane(hs + [((0, 5), (1, 5))])  # also y >= 5: none
  []
  >>> # Largest circle inside: shift edges left by r, bisect
  >>> def shift(s, e, r):  # the line s -> e moved left by r
  ...   k = r / math.hypot(e[0] - s[0], e[1] - s[1])
  ...   ox, oy = -(e[1] - s[1]) * k, (e[0] - s[0]) * k
  ...   return (s[0] + ox, s[1] + oy), (e[0] + ox, e[1] + oy)
  >>> tri = [(0, 0), (4, 0), (0, 3)]  # inradius 1
  >>> ts = [(tri[i - 1], tri[i]) for i in range(3)]
  >>> lo, hi = 0.0, 4.0
  >>> for _ in range(50):
  ...   mid = (lo + hi) / 2
  ...   ok = half_plane([shift(s, e, mid) for s, e in ts])
  ...   lo, hi = (mid, hi) if ok else (lo, mid)
  >>> round(lo, 6)
  1.0
Notes: Decide whether touching counts. This code drops results
  without interior, so to find a point on the closed side of
  every line (a single point may be the answer), first move each
  line outward by a little more than `eps`. A vertex on the box
  means the true intersection is unbounded. Doubles near the box
  at $10^9$ are only good to about $10^(-7)$: without the
  rounding term in `out`, a line through a box corner would
  split it into two vertices. Linear programs in two variables
  are half-plane intersections: maximise the objective over the
  vertices. "A point that sees every wall" (polygon kernel),
  "the region closer to site $a$ than to all others" (a Voronoi
  cell: keep the $a$ side of each perpendicular bisector) and
  "the largest circle or offset that fits" are the usual
  disguises. Integer input with exact answers needs a rational
  version: every test is the sign of a degree-4 polynomial in
  the coordinates (exact with Python ints).
Source: Zhu Zeyuan's sort-and-deque algorithm; cp-algorithms
  Half-plane intersection (idea).
Status: stress-tested against an exact brute force
"""

import math
from collections import deque


def half_plane(hs, big=1e9):
  eps = 1e-9
  c = [(-big, -big), (big, -big), (big, big), (-big, big)]
  sx, sy, dx, dy, ln = [], [], [], [], []
  for s, e in list(hs) + [(c[i - 1], c[i]) for i in range(4)]:
    sx.append(s[0])
    sy.append(s[1])
    dx.append(e[0] - s[0])
    dy.append(e[1] - s[1])
    ln.append(math.hypot(dx[-1], dy[-1]))
  # + 0.0 turns -0.0 into 0.0
  ang = [math.atan2(y + 0.0, x) for x, y in zip(dx, dy)]

  def at(i, j):  # as line_inter, exact on integer corners
    d = dx[i] * dy[j] - dy[i] * dx[j]
    t = dx[j] * (sy[i] - sy[j]) - dy[j] * (sx[i] - sx[j])
    x = (sx[i] * d + dx[i] * t) / d
    return x, (sy[i] * d + dy[i] * t) / d

  def out(h, p):  # p not more than eps left of line h
    x, y = p[0] - sx[h], p[1] - sy[h]
    tol = eps + 1e-15 * (abs(p[0]) + abs(p[1]))  # + rounding
    return dx[h] * y - dy[h] * x < tol * ln[h]

  dq = deque()
  for h in sorted(range(len(ang)), key=ang.__getitem__):
    while len(dq) > 1 and out(h, at(dq[-2], dq[-1])):
      dq.pop()
    while len(dq) > 1 and out(h, at(dq[0], dq[1])):
      dq.popleft()
    if dq:
      b = dq[-1]
      turn = dx[h] * dy[b] - dy[h] * dx[b]
      if abs(turn) < eps * ln[h] * ln[b]:  # parallel
        if dx[h] * dx[b] + dy[h] * dy[b] < 0:
          return []
        if not out(h, (sx[b], sy[b])):
          continue
        dq.pop()
    dq.append(h)
  while len(dq) > 2 and out(dq[0], at(dq[-2], dq[-1])):
    dq.pop()
  while len(dq) > 2 and out(dq[-1], at(dq[0], dq[1])):
    dq.popleft()
  q = list(dq)
  if len(q) < 3:
    return []
  return [at(q[i], q[(i + 1) % len(q)]) for i in range(len(q))]
