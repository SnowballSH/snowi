"""
Title: Segment tree
Description: Bottom-up point assign, range query over a monoid.
  The largest `r` with `pred(query(l, r))` true is
  `max_right(l, pred)`, and the smallest such `l` is
  `min_left(r, pred)`.
Precondition: 0-indexed, half-open ranges. `op` need not be
  commutative (left operand is the left part). `pred(e)` is true
  and `pred` is monotone: once false for a range, false for
  every range containing it.
Time: $O(n)$ build, $O(log n)$ per operation.
Usage:
  st = SegTree(a, add, 0)  # from operator import add
  st.set(i, 5); st.get(i); st.query(l, r); st.query_all()
  # position of the k-th one (k >= 1) in a 0/1 list, or n:
  i = st.max_right(0, lambda s: s < k)
Tweak: Max: `SegTree(a, max, -INF)`; then
  `max_right(l, lambda v: v < x)` is the first index `>= l` with
  `a[i] >= x`. Min: `SegTree(a, min, INF)`. Speed: inline `op` into
  `set`/`query` (30%); pack pairs as `(x << 20) + y` below $2^63$
  (tuple nodes are 3 times slower).
Example:
  >>> st = SegTree([3, 1, 4, 1, 5], max, -1)
  >>> st.query(1, 4), st.query_all()
  (4, 5)
  >>> st.max_right(1, lambda v: v < 4)  # first a[i] >= 4, i >= 1
  2
  >>> st.set(2, 0); st.min_left(4, lambda v: v < 3)
  1
  >>> MOD = 998244353
  >>> def op(f, g):  # maps x -> b x + c: "f, then g"
  ...   return g[0] * f[0] % MOD, (g[0] * f[1] + g[1]) % MOD
  >>> fs = SegTree([(2, 1), (3, 0), (1, 5)], op, (1, 0))
  >>> b, c = fs.query(0, 3); b * 10 + c  # (2*10 + 1) * 3 + 5
  68
Notes: Think of it when a problem has point updates and range
  questions about any associative summary, or asks for "the first
  position after `l` where ...": `max_right` walks down the tree
  in $O(log n)$ instead of binary searching with queries. On sums,
  `s <= x` is monotone only when all values are `>= 0`. The query
  keeps a left and a right accumulator, so order is preserved:
  swapping them is invisible for max or sum but wrong for maps.
  Leaves past `n` hold `e`, so `e` must be a true identity. Use an
  int `INF` such as `1 << 62`, not `math.inf`: a float among the
  ints made a max tree 1.6 times slower.
Source: AtCoder Library segtree (CC0); kactl (CC0).
Status: stress-tested against brute force
"""


class SegTree:
  def __init__(self, values, op, e):
    n, size = len(values), 1
    while size < n:
      size *= 2
    t = [e] * (2 * size)
    t[size : size + n] = values
    for i in range(size - 1, 0, -1):
      t[i] = op(t[2 * i], t[2 * i + 1])
    self.n, self.size, self.t = n, size, t
    self.op, self.e = op, e

  def set(self, i, v):
    t, op = self.t, self.op
    i += self.size
    t[i] = v
    i >>= 1
    while i:
      t[i] = op(t[2 * i], t[2 * i + 1])
      i >>= 1

  def get(self, i):
    return self.t[i + self.size]

  def query_all(self):
    return self.t[1]

  def query(self, l, r):
    t, op = self.t, self.op
    a = b = self.e
    l += self.size
    r += self.size
    while l < r:
      if l & 1:
        a = op(a, t[l])
        l += 1
      if r & 1:
        r -= 1
        b = op(t[r], b)
      l >>= 1
      r >>= 1
    return op(a, b)

  def max_right(self, l, pred):
    if l == self.n:
      return l
    t, op, size, acc = self.t, self.op, self.size, self.e
    l += size
    while True:
      while l % 2 == 0:
        l >>= 1
      if not pred(op(acc, t[l])):
        while l < size:
          l *= 2
          if pred(op(acc, t[l])):
            acc = op(acc, t[l])
            l += 1
        return l - size
      acc = op(acc, t[l])
      l += 1
      if (l & -l) == l:
        return self.n

  def min_left(self, r, pred):
    if r == 0:
      return 0
    t, op, size, acc = self.t, self.op, self.size, self.e
    r += size
    while True:
      r -= 1
      while r > 1 and r % 2:
        r >>= 1
      if not pred(op(t[r], acc)):
        while r < size:
          r = 2 * r + 1
          if pred(op(t[r], acc)):
            acc = op(t[r], acc)
            r -= 1
        return r + 1 - size
      acc = op(t[r], acc)
      if (r & -r) == r:
        return 0
