"""
Title: Ternary search
Description: `ternary(lo, hi, f)` returns the smallest `x` in
  `[lo, hi]` that minimises an integer function `f`, by binary
  search for the first `x` with `f(x) <= f(x + 1)`.
  `golden(a, b, f)` returns a minimiser of a real function on
  `[a, b]` by golden-section search.
Precondition: `ternary`: `lo <= hi`; `f` strictly decreases up to
  its first minimum and never decreases after it. Every convex
  `f` qualifies (a flat bottom is fine); a flat stretch before
  the minimum breaks it. `golden`: `f` strictly decreases, then
  strictly increases, except on a flat bottom.
Time: `ternary`: $2 ceil(log_2 n)$ calls of `f` for $n$ values.
  `golden`: 102 calls; each shrinks the interval by 0.618.
Usage:
  x = ternary(0, n - 1, cost)
  t = golden(0, 1e9, g)
Tweak: `f` may return anything comparable (ints, floats, tuples).
Example:
  >>> ternary(-10, 10, lambda x: (x - 3) ** 2)
  3
  >>> ternary(0, 10, lambda x: abs(x - 4) + abs(x - 7))
  4
  >>> # Points move along a line, x_i(t) = p[i] + v[i] t. The
  >>> # spread max - min is convex in t: when are they closest?
  >>> p, v = [0, 10, 3], [2, -1, 0]
  >>> def spread(t):
  ...   xs = [a + b * t for a, b in zip(p, v)]
  ...   return max(xs) - min(xs)
  >>> t = golden(0, 100, spread)
  >>> round(t, 6), round(spread(t), 6)
  (3.333333, 3.666667)
Notes: Unimodality does not survive sums but convexity does:
  sums and maxima of convex functions are convex ($sum_i |x -
  a_i|$, the farthest of several points, the spread of moving
  points). If $F(x, y)$ is convex in both together,
  $min_y F(x, y)$ is convex in $x$: nest two searches. A product
  of monotone functions such as $i dot f(i)$ is often not
  unimodal: print `f` on small cases before trusting it. Near a
  smooth minimum `f` is flat, so `golden` finds the minimiser
  only to about $10^(-8)$ relative error, but the minimum value
  to full precision: print `f(x)`, not a formula in `x`.
  On integers prefer `ternary`; rounding `golden`'s answer can
  miss the best integer. Golden section costs one call per step
  (factor 0.618), plain ternary search two (factor 0.667 per
  step, 0.816 per call). Each call of `f` is a Python call: keep
  `f` a plain function over flat lists. Cues: "choose a time,
  position or parameter minimising the maximum (or the total)
  cost", "the answer first improves, then gets worse".
Source: kactl (CC0).
Status: stress-tested against brute force
"""


def ternary(lo, hi, f):
  while lo < hi:
    mid = (lo + hi) // 2
    if f(mid) <= f(mid + 1):
      hi = mid
    else:
      lo = mid + 1
  return lo


def golden(a, b, f):
  r = (5**0.5 - 1) / 2
  x1, x2 = b - r * (b - a), a + r * (b - a)
  f1, f2 = f(x1), f(x2)
  for _ in range(100):
    if f1 < f2:
      b, x2, f2 = x2, x1, f1
      x1 = b - r * (b - a)
      f1 = f(x1)
    else:
      a, x1, f1 = x1, x2, f2
      x2 = a + r * (b - a)
      f2 = f(x2)
  return (a + b) / 2
