"""
Title: Minkowski sum
Description: Returns the Minkowski sum $A + B = {a + b}$ of two
  convex polygons, counter-clockwise from its smallest point (by
  `x`, then `y`).
Precondition: `a`, `b` non-empty, convex, counter-clockwise and
  without collinear vertices, as from `convex_hull` (one or two
  points are fine). The sum is then the same kind of polygon, with
  at most $n + m$ vertices. Exact for integers.
Time: $O(n + m)$.
Usage:
  s = minkowski(a, b)
Tweak: $A - B$: negate the points of `b` first (a half turn keeps
  the order counter-clockwise). Sum of $k$ polygons: sum them in
  pairs, or merge all edges by angle at once.
Example:
  >>> tri = [(0, 0), (2, 0), (0, 1)]
  >>> sq = [(0, 0), (1, 0), (1, 1), (0, 1)]
  >>> minkowski(tri, sq)
  [(0, 0), (3, 0), (3, 1), (1, 2), (0, 2)]
  >>> # Do A and B + t meet (touching counts)? Iff t is in A - B.
  >>> d = minkowski(sq, [(-x, -y) for x, y in tri])
  >>> def meet(t):
  ...   return all(cross3(d[i - 1], d[i], t) >= 0
  ...              for i in range(len(d)))
  >>> meet((0, 0)), meet((1, 1)), meet((2, 0))
  (True, True, False)
Notes: The extreme point of $A + B$ in a direction is the sum of
  the extreme points of $A$ and $B$, so "maximise a linear function
  of $a + b$" is one query on the sum. The translations $t$ for
  which $A$ and $B + t$ overlap form $A - B$ (collision checks,
  free space of a moving robot); the distance between convex
  polygons is that from the origin to $A - B$. A sum has one vertex
  per distinct edge direction of the summands. Window problems
  ("best sum + length times $x$ over all subarrays") become upper
  hulls of sums $(i, S_i) + (-j, -S_j)$, built by divide and
  conquer on the index with merges like this one. Reverse
  clockwise input first; the merge itself starts both polygons at
  their smallest point.
Uses: geometry/point
Source: UCup S4-S5 Training Packet 20; cp-algorithms Minkowski
  sum (idea).
Status: stress-tested against brute force
"""


def minkowski(a, b):
  k = a.index(min(a))
  a = a[k:] + a[: k + 1]
  k = b.index(min(b))
  b = b[k:] + b[: k + 1]
  n, m, i, j = len(a) - 1, len(b) - 1, 0, 0
  res = []
  while i < n or j < m:
    res.append((a[i][0] + b[j][0], a[i][1] + b[j][1]))
    c = 0  # > 0: a's edge turns less, take it first
    if i < n and j < m:
      c = cross(sub(a[i + 1], a[i]), sub(b[j + 1], b[j]))
    if i < n and c >= 0:
      i += 1
    if j < m and c <= 0:
      j += 1
  return res
