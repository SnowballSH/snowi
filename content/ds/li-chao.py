"""
Title: Li Chao tree
Description: Minimum of lines $k x + m$ at the points of a fixed
  sorted list `xs`. `add_segment` adds a line that exists only
  for `x` in `[l, r)`. `query(x)` returns `LiChao.INF` if no line
  covers `x`. The convex hull trick for PyPy.
Precondition: `xs` is sorted and non-empty; `query(x)` needs `x`
  in `xs`; `l` and `r` may be any values. Values $k x + m$ stay
  below `INF` $= 2^63 - 1$ (raise `INF` for larger ones).
Time: $O(log n)$ per line or query, $O(log^2 n)$ per segment.
Usage:
  lc = LiChao(xs)  # the query points
  lc.add_line(k, m); lc.add_segment(k, m, l, r)
  y = lc.query(x)
Tweak: Maximum: add `(-k, -m)` and negate the answer. Online `x`
  in a large range: use `xs = range(lo, hi)` if it fits in memory,
  else read the queries first and compress.
Example:
  >>> lc = LiChao(list(range(8)))
  >>> lc.add_line(1, 0); lc.add_line(-1, 5)
  >>> [lc.query(x) for x in range(8)]
  [0, 1, 2, 2, 1, 0, -1, -2]
  >>> lc.add_segment(0, -5, 2, 4)  # y = -5 on [2, 4)
  >>> [lc.query(x) for x in range(8)]
  [0, 1, -5, -5, 1, 0, -1, -2]
  >>> h, C = [1, 3, 4, 7], 2  # frog: (h[j] - h[i])^2 + C
  >>> lc, dp = LiChao(sorted(h)), [0]
  >>> lc.add_line(-2 * h[0], h[0] ** 2)
  >>> for j in range(1, 4):
  ...   dp.append(lc.query(h[j]) + h[j] ** 2 + C)
  ...   lc.add_line(-2 * h[j], dp[j] + h[j] ** 2)
  >>> dp
  [0, 6, 9, 20]
Notes: Think of it when a DP transition expands to
  $f(i) g(j) + a(i) + b(j)$: each `i` is a line with slope $f(i)$
  and intercept $a(i)$, queried at $x = g(j)$. Signals: squared
  differences $(h_j - h_i)^2$, products of prefix sums, "for
  every query `c`, the best value minus `c` times size". Slopes
  and queries may come in any order, and segments handle lines
  valid on a range of `x` only. Add line `j` after computing
  `dp[j]`, never before. Python cannot overflow, but values
  beyond $2^63$ become big ints and are several times slower.
  Each node keeps the line that wins at its middle point; the
  loser can win on one side only and moves down into that child.
  A query takes the minimum along the leaf-to-root path.
Source: Li Chao tree as in cp-algorithms; Packet 06 (UCup).
Status: stress-tested against brute force
"""

from bisect import bisect_left


class LiChao:
  INF = (1 << 63) - 1

  def __init__(self, xs):
    size = 1
    while size < len(xs):
      size *= 2
    self.size, self.lg = size, size.bit_length() - 1
    self.xs = list(xs) + [xs[-1]] * (size - len(xs))
    self.K = [0] * (2 * size)  # node v holds K[v] x + M[v]
    self.M = [self.INF] * (2 * size)

  def _push(self, k, m, v, h):  # node v has height h
    xs, K, M = self.xs, self.K, self.M
    lo = (v << h) - self.size
    while True:
      x = xs[lo]
      a = k * x + m < K[v] * x + M[v]
      x = xs[lo + (1 << h) - 1]
      if a == (k * x + m < K[v] * x + M[v]):
        if a:  # one line wins on the whole node
          K[v], M[v] = k, m
        return
      h -= 1
      mid = lo + (1 << h)
      x = xs[mid]
      if k * x + m < K[v] * x + M[v]:
        k, K[v] = K[v], k
        m, M[v] = M[v], m
        a = not a
      v *= 2
      if not a:  # the loser can only win on the right
        v += 1
        lo = mid

  def add_line(self, k, m):
    self._push(k, m, 1, self.lg)

  def add_segment(self, k, m, l, r):  # x in [l, r)
    size, xs = self.size, self.xs
    l = bisect_left(xs, l) + size
    r = bisect_left(xs, r) + size
    h = 0
    while l < r:
      if l & 1:
        self._push(k, m, l, h)
        l += 1
      if r & 1:
        r -= 1
        self._push(k, m, r, h)
      l >>= 1
      r >>= 1
      h += 1

  def query(self, x):
    K, M = self.K, self.M
    v = bisect_left(self.xs, x) + self.size
    y = self.INF
    while v:
      w = K[v] * x + M[v]
      if w < y:
        y = w
      v >>= 1
    return y
