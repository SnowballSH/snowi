"""
Title: Closest pair of points
Description: Returns two points of `pts` at the smallest distance.
  Divide and conquer, bottom-up: with the points sorted by `x`,
  merge sort them by `y` in blocks of doubling size; each merge
  compares the points near the middle `x` with their next
  neighbours by `y`.
Precondition: At least two points, tuples `(x, y, ...)`; extra
  fields ride along, so `(x, y, i)` gives the indices. Equal points
  are allowed (distance 0). Exact for integers.
Time: $O(n log n)$; under a second for $2 dot 10^5$ points.
Usage:
  a, b = closest_pair(pts); d2 = dist2(a, b)
Tweak: Floats work as is. For the indices pass `(x, y, i)`.
Example:
  >>> pts = [(0, 0), (5, 1), (2, 3), (6, 3), (0, 5)]
  >>> a, b = closest_pair(pts); sorted([a, b]), dist2(a, b)
  ([(5, 1), (6, 3)], 5)
  >>> # Tricky function: min of (i - j)^2 + (S_i - S_j)^2
  >>> a = [1, -3, 4, -1, 2]
  >>> s = [sum(a[:i + 1]) for i in range(len(a))]
  >>> p, q = closest_pair([(i, x) for i, x in enumerate(s)])
  >>> dist2(p, q)
  2
Notes: "Minimum distance between any two" of $10^5$ or more
  points, or "minimise $(i - j)^2 + g(i, j)^2$" (a distance in
  disguise), call for this. Manhattan distance becomes the $L_oo$
  distance $max(|Delta x|, |Delta y|)$ of the points
  $(x + y, x - y)$: measure that instead, and compare `best` with
  $|Delta x|$ and $|Delta y|$ instead of their squares. The closest
  red-blue pair cannot use this (points of one color may be close
  together, so the strips are not small); it needs a Voronoi
  diagram of one color. Key sorts are slow in PyPy, so the merges
  are written out; only one sort by `x` remains.
Uses: geometry/point
Source: Shamos and Hoey's divide and conquer, made iterative.
Status: stress-tested against brute force
"""


def closest_pair(pts):
  n = len(pts)
  xs = [p[0] for p in pts]
  order = sorted(range(n), key=xs.__getitem__)
  xs = [xs[i] for i in order]
  ys = [pts[i][1] for i in order]
  best, bi, bj = float("inf"), 0, 1
  idx = list(range(n))  # each aligned block sorted by y
  w = 1
  while w < n:  # merge blocks [lo, m) and [m, hi) by y
    for lo in range(0, n - w, 2 * w):
      m, hi = lo + w, min(lo + 2 * w, n)
      blk, i, j = [], lo, m
      while i < m or j < hi:
        if j == hi or (i < m and ys[idx[i]] <= ys[idx[j]]):
          blk.append(idx[i])
          i += 1
        else:
          blk.append(idx[j])
          j += 1
      idx[lo:hi] = blk
      xm = xs[m]
      strip = [t for t in blk if (xs[t] - xm) ** 2 < best]
      for a in range(len(strip)):
        i = strip[a]
        for b in range(a + 1, len(strip)):
          j = strip[b]
          dy = ys[j] - ys[i]
          if dy * dy >= best:
            break
          dx = xs[j] - xs[i]
          if dx * dx + dy * dy < best:
            best, bi, bj = dx * dx + dy * dy, i, j
    w *= 2
  return pts[order[bi]], pts[order[bj]]
