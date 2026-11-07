"""
Title: Circles
Description: `circle_line` meets circle `(c, r)` with line `ab`: no
  point, one (tangent) or two, in the direction `a -> b`.
  `circle_inter` meets two circles: no point, one (touching) or two,
  the first left of the line from `a` to `b`. `tangents` gives the
  common tangents as pairs of touch points (on circle 1, on circle
  2): the outer ones, the inner ones if `r2` is negated, and the
  tangents from point `c2` if `r2 = 0`. `circumcircle` returns the
  center and radius of the circle through three points, `mec` those
  of the smallest circle containing all points (Welzl, randomised).
  `circle_poly` is the signed area of the circle's intersection with
  a simple polygon (positive if counter-clockwise).
Precondition: Results are floats. Lengths within `eps` $= 10^(-9)$
  count as equal (fine for coordinates up to about $10^4$).
  Concentric circles give no point and no tangent. `circumcircle`:
  points not collinear. `mec`: at least one point. `circle_poly`: no
  two consecutive vertices equal.
Time: $O(1)$; `mec` expected $O(n)$, `circle_poly` $O(n)$.
Usage:
  v = circle_line(c, r, a, b)
  w = circle_inter(c1, r1, c2, r2)
  t = tangents(c1, r1, c2, r2)  # t[i][0] on circle 1
  o, rad = circumcircle(a, b, c)
  center, radius = mec(pts)
  area = circle_poly(c, r, poly)
Tweak: For larger coordinates scale `eps` (to about
  $10^(-13)$ times the largest coordinate); `mec` has no `eps`.
Example:
  >>> circle_line((0, 0), 5, (3, -10), (3, 10))
  [(3.0, -4.0), (3.0, 4.0)]
  >>> circle_inter((0, 0), 5, (8, 6), 5)
  [(4.0, 3.0)]
  >>> [t for t, _ in tangents((0, 0), 5, (0, 10), 5)]
  [(5.0, 0.0), (-5.0, 0.0)]
  >>> circumcircle((0, 0), (4, 0), (0, 3))
  ((2.0, 1.5), 2.5)
  >>> mec([(0, 0), (2, 0), (1, 1), (1, -1), (1, 0)])
  ((1.0, 0.0), 1.0)
  >>> sq = [(0, 0), (2, 0), (2, 2), (0, 2)]
  >>> round(circle_poly((0, 0), 1, sq), 9)  # a quarter circle
  0.785398163
Notes: `mec` works relative to its first point, so its precision
  follows the spread of the points, not their distance from the
  origin; it counts points up to `r * 1e-9` outside as inside. With
  integer input, decide the case exactly in integers (circles touch
  iff $d^2 = (r_1 plus.minus r_2)^2$; line `ab` is tangent iff
  $((b - a) times (c - a))^2 = r^2 |b - a|^2$) and use these
  functions only for the points. The minimum enclosing circle hides
  behind "place one station to minimise the farthest distance". For
  the largest circle inside a convex polygon, binary search the
  radius, push every edge inward by it and test whether the
  half-plane intersection is empty. Circumradius from side lengths
  $x, y, z$ and area $S$: $R = x y z / (4 S)$.
Source: kactl CircleLine, CircleIntersection, CircleTangents,
  Circumcircle, MinimumEnclosingCircle, CirclePolygonIntersection
  (CC0).
Status: stress-tested against brute force
"""

import math
import random


def circle_line(c, r, a, b):
  eps = 1e-9
  ux, uy = b[0] - a[0], b[1] - a[1]
  n2 = ux * ux + uy * uy
  t = ((c[0] - a[0]) * ux + (c[1] - a[1]) * uy) / n2
  x, y = a[0] + ux * t, a[1] + uy * t
  d = math.hypot(c[0] - x, c[1] - y)
  if d > r + eps:
    return []
  if d > r - eps:
    return [(x, y)]
  k = math.sqrt((r * r - d * d) / n2)
  return [(x - ux * k, y - uy * k), (x + ux * k, y + uy * k)]


def circle_inter(a, r1, b, r2):
  eps = 1e-9
  vx, vy = b[0] - a[0], b[1] - a[1]
  d, hi, lo = math.hypot(vx, vy), r1 + r2, abs(r1 - r2)
  if d < eps or d > hi + eps or d < lo - eps:
    return []
  t = (d * d + r1 * r1 - r2 * r2) / (2 * d * d)
  x, y = a[0] + vx * t, a[1] + vy * t
  if d > hi - eps or d < lo + eps:
    return [(x, y)]
  h = math.sqrt((d - lo) * (d + lo) * (hi - d) * (hi + d))
  h /= 2 * d * d
  return [(x - vy * h, y + vx * h), (x + vy * h, y - vx * h)]


def tangents(c1, r1, c2, r2):
  eps = 1e-9
  dx, dy = c2[0] - c1[0], c2[1] - c1[1]
  dr, d2 = r1 - r2, dx * dx + dy * dy
  d = math.sqrt(d2)
  if d < eps or d < abs(dr) - eps:
    return []
  h = 0.0 if d < abs(dr) + eps else math.sqrt(d2 - dr * dr)
  out = []
  for s in (-1, 1) if h else (-1,):
    vx = (dx * dr - dy * h * s) / d2
    vy = (dy * dr + dx * h * s) / d2
    p1 = (c1[0] + vx * r1, c1[1] + vy * r1)
    out.append((p1, (c2[0] + vx * r2, c2[1] + vy * r2)))
  return out


def circumcircle(a, b, c):
  ux, uy = c[0] - a[0], c[1] - a[1]
  vx, vy = b[0] - a[0], b[1] - a[1]
  u2, v2 = ux * ux + uy * uy, vx * vx + vy * vy
  w = 2 * (ux * vy - uy * vx)
  ox = a[0] - (uy * v2 - vy * u2) / w
  oy = a[1] + (ux * v2 - vx * u2) / w
  return (ox, oy), math.hypot(ox - a[0], oy - a[1])


def mec(ps):
  hypot, slack = math.hypot, 1 + 1e-9
  ps = random.sample(ps, len(ps))
  bx, by = ps[0][0], ps[0][1]
  ps = [(p[0] - bx, p[1] - by) for p in ps]
  ox, oy, r = 0.0, 0.0, 0.0
  for i in range(len(ps)):
    pi = ps[i]
    if hypot(pi[0] - ox, pi[1] - oy) <= r * slack:
      continue
    ox, oy, r = pi[0], pi[1], 0.0
    for j in range(i):
      pj = ps[j]
      if hypot(pj[0] - ox, pj[1] - oy) <= r * slack:
        continue
      ox, oy = (pi[0] + pj[0]) / 2, (pi[1] + pj[1]) / 2
      r = hypot(pi[0] - ox, pi[1] - oy)
      for k in range(j):
        pk = ps[k]
        if hypot(pk[0] - ox, pk[1] - oy) > r * slack:
          (ox, oy), r = circumcircle(pi, pj, pk)
  return (bx + ox, by + oy), r


def circle_poly(c, r, poly):
  k, total = r * r / 2, 0.0

  def arc(ax, ay, bx, by):  # sector of the turn a -> b
    return math.atan2(ax * by - ay * bx, ax * bx + ay * by) * k

  for i in range(len(poly)):
    px, py = poly[i - 1][0] - c[0], poly[i - 1][1] - c[1]
    qx, qy = poly[i][0] - c[0], poly[i][1] - c[1]
    dx, dy = qx - px, qy - py
    d2 = dx * dx + dy * dy
    a, s, t = (dx * px + dy * py) / d2, 1, 0
    h = a * a - (px * px + py * py - r * r) / d2
    if h > 0:  # the chord p + d [s, t] lies inside
      h = math.sqrt(h)
      s, t = max(0.0, -a - h), min(1.0, -a + h)
    if t < 0 or s >= 1:
      total += arc(px, py, qx, qy)
    else:
      ux, uy = px + dx * s, py + dy * s
      wx, wy = px + dx * t, py + dy * t
      total += arc(px, py, ux, uy) + arc(wx, wy, qx, qy)
      total += (ux * wy - uy * wx) / 2
  return total
