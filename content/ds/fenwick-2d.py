"""
Title: 2D Fenwick tree
Description: Point add and rectangle sums on a dense `n` x `m`
  grid. `query(x, y)` is the sum over `[0, x) x [0, y)`;
  `rect(x1, y1, x2, y2)` the sum over `[x1, x2) x [y1, y2)`.
Precondition: 0-indexed, `0 <= x <= n`, `0 <= y <= m`; `add` at
  `x = n` or `y = m` does nothing. Memory is $n m$ values in one
  flat list, so $n m <= 10^7$ or so.
Time: $O(log n log m)$ per operation.
Usage:
  fw = Fenwick2D(n, m)
  fw.add(x, y, 5)
  s = fw.rect(x1, y1, x2, y2)
Tweak: Rectangle add with point query: add `v` at `(x1, y1)` and
  `(x2, y2)`, `-v` at `(x1, y2)` and `(x2, y1)`; the value at
  `(x, y)` is `query(x + 1, y + 1)`.
Example:
  >>> fw = Fenwick2D(3, 4)
  >>> fw.add(0, 0, 1); fw.add(1, 2, 5); fw.add(2, 3, 7)
  >>> fw.query(2, 3), fw.rect(1, 2, 3, 4), fw.rect(0, 1, 3, 2)
  (6, 12, 0)
Notes: For sparse points with big coordinates, go offline: sort
  the distinct `x` of all updates, and for each Fenwick node `i`
  collect (sorted, unique) the `y` of updates whose `x` passes
  through `i`. Each node then holds a 1D Fenwick tree over its
  own `y` list, searched with `bisect`: $O(log^2 n)$ per operation
  and $O(k log k)$ memory for `k` updates. If there are no updates
  between queries, a sweep over `x` with a 1D Fenwick tree over
  `y` is simpler and faster. Cues: "count points in a rectangle",
  "sum of a submatrix with point changes".
Source: kactl FenwickTree2d (CC0) for the offline variant.
Status: stress-tested against brute force
"""


class Fenwick2D:
  def __init__(self, n, m):
    self.n, self.m = n, m
    self.t = [0] * (n * m)

  def add(self, x, y, delta):
    t, n, m = self.t, self.n, self.m
    while x < n:
      j, row = y, x * m
      while j < m:
        t[row + j] += delta
        j |= j + 1
      x |= x + 1

  def query(self, x, y):  # sum of [0, x) x [0, y)
    t, m, s = self.t, self.m, 0
    while x > 0:
      j, row = y, (x - 1) * m - 1
      while j > 0:
        s += t[row + j]
        j &= j - 1
      x &= x - 1
    return s

  def rect(self, x1, y1, x2, y2):
    q = self.query
    return q(x2, y2) - q(x1, y2) - q(x2, y1) + q(x1, y1)
