"""
Title: Point and vector
Description: Points and vectors are plain tuples `(x, y)`.
  `cross3(o, a, b) > 0` iff `o -> a -> b` turns left. `angle_sort`:
  polar angle order from the positive x-axis, counter-clockwise;
  `angle_cmp` is the same order as a three-way comparison. Both are
  exact for integers.
Precondition: Extra tuple fields after `x, y` are ignored, so
  `(x, y, i)` carries an index along. Angle order: no zero vectors;
  equal directions tie.
Time: $O(1)$ each; `angle_sort` $O(n log n)$.
Usage:
  o, a, b = (0, 0), (3, 1), (x, y)
  c, turn = cross(a, b), cross3(o, a, b)
  d2 = dist2(a, b)  # squared distance
  angle_sort(vecs)  # in place
Tweak: Real distance: `math.sqrt(dist2(a, b))` (`math.dist` would
  also count extra tuple fields).
Example:
  >>> v = [(1, 0), (0, -2), (-1, 0), (1, 1), (2, 2)]
  >>> angle_sort(v); v
  [(1, 0), (1, 1), (2, 2), (-1, 0), (0, -2)]
  >>> cross3((0, 0), (1, 0), (1, 1)), dist2((1, 1), (4, 5))
  (1, 25)
  >>> sorted([(0, 1), (1, 0)], key=cmp_to_key(angle_cmp))
  [(1, 0), (0, 1)]
Notes: The sign of a cross product answers "left or right turn?",
  "which side of a line?" and "clockwise or not?"; with integer
  input it never needs an epsilon. Compare lengths with `dist2` to
  stay in integers. Sorting by `atan2` alone is wrong for
  coordinates above about $10^7$ (distinct directions round to the
  same float); `angle_sort` sorts by `atan2`, then sorts each run of
  keys within $10^(-12)$ of each other with `angle_cmp`. That is
  about 5 times faster than `cmp_to_key` alone, and no worse than it
  when every key ties. To sort around a center `c`, sort
  `(x - cx, y - cy, i)`. In hot loops inline the arithmetic on local
  variables instead of calling these.
Source: kactl Point (CC0); half-plane angle order is folklore.
Status: stress-tested against brute force
"""

import math
from functools import cmp_to_key


def sub(a, b):
  return (a[0] - b[0], a[1] - b[1])


def dot(a, b):
  return a[0] * b[0] + a[1] * b[1]


def cross(a, b):
  return a[0] * b[1] - a[1] * b[0]


def cross3(o, a, b):
  ox, oy = o[0], o[1]
  return (a[0] - ox) * (b[1] - oy) - (a[1] - oy) * (b[0] - ox)


def dist2(a, b):
  dx, dy = a[0] - b[0], a[1] - b[1]
  return dx * dx + dy * dy


def half(a):  # angle in [pi, 2 pi)
  return a[1] < 0 or (a[1] == 0 and a[0] < 0)


def angle_cmp(a, b):
  ha, hb = half(a), half(b)
  if ha != hb:
    return 1 if ha else -1
  c = cross(a, b)
  return (c < 0) - (c > 0)


def angle_sort(v):
  tau, atan2 = 2 * math.pi, math.atan2
  v.sort(key=lambda p: atan2(p[1], p[0]) % tau)
  k = [atan2(p[1], p[0]) % tau for p in v]
  i, n = 0, len(v)
  while i < n:  # sort each run of float near-ties exactly
    j = i + 1
    while j < n and k[j] - k[j - 1] < 1e-12:
      j += 1
    if j - i > 1:
      v[i:j] = sorted(v[i:j], key=cmp_to_key(angle_cmp))
    i = j
