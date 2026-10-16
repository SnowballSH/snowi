"""
Title: Lazy segment tree
Description: Range apply and range fold. Nodes form a monoid under
  `op` with identity `e`; tags act on nodes through `mapping(f, s)`
  and compose by `composition(f, g)` = "`g`, then `f`", with
  identity `id`.
Precondition: 0-indexed, half-open ranges. `mapping(f, op(a, b))
  = op(mapping(f, a), mapping(f, b))` and `mapping(id, s) = s`.
  `pred(e)` is true and `pred` is monotone.
Time: $O(n)$ build, $O(log n)$ per operation.
Usage:
  from operator import add  # range add, range min:
  st = LazySegTree(a, min, INF, add, add, 0)
  st.apply(l, r, 5); st.query(l, r); st.query_all()
  st.set(i, 7); st.get(i)
  r = st.max_right(l, lambda v: v > x)
Tweak: Range sum needs the length: nodes `(sum, len)`, leaves
  `(a_i, 1)`, `mapping` returns `(s[0] + f * s[1], s[1])`. Assign:
  `id = None`, `mapping` returns `s if f is None else f`,
  `composition` `g if f is None else f`. Pack `(sum, len)` as
  `(sum << 20) + len` while $|"sum"| < 2^43$: `op = add`,
  `mapping(f, s) = s + (f * (s & M) << 20)`, `M = 2**20 - 1`.
Example:
  >>> MOD = 998244353  # range affine a_i = b a_i + c, range sum
  >>> def op(x, y):  # node: (sum, len)
  ...   return (x[0] + y[0]) % MOD, x[1] + y[1]
  >>> def mapping(f, s):  # tag: (b, c)
  ...   return (f[0] * s[0] + f[1] * s[1]) % MOD, s[1]
  >>> def comp(f, g):  # g, then f
  ...   return f[0] * g[0] % MOD, (f[0] * g[1] + f[1]) % MOD
  >>> a = [(x, 1) for x in [1, 2, 3, 4]]
  >>> st = LazySegTree(a, op, (0, 0), mapping, comp, (1, 0))
  >>> st.apply(1, 4, (2, 1)); st.apply(0, 2, (1, 10))
  >>> [st.get(i)[0] for i in range(4)], st.query(1, 3)
  ([11, 15, 7, 9], (22, 2))
Notes: Check the three laws on paper before coding: tags compose,
  act on a whole node in $O(1)$, and distribute over `op`. Sums
  need `len`; products of tags must stay in a small family. The
  new tag goes on the outside: `composition(new, old)`. If a tag
  cannot update a node's summary without looking inside (chmin
  with sum, `%=` with sum), use segment tree beats. Point apply is
  `apply(i, i + 1, f)`. `max_right` finds e.g. the first index
  after `l` with value `<= x` under range adds: `pred = v > x` on
  the min.
Source: AtCoder Library lazy_segtree (CC0).
Status: stress-tested against brute force
"""


class LazySegTree:
  def __init__(self, values, op, e, mapping, composition, id):
    n = len(values)
    lg = max(n - 1, 0).bit_length()
    size = 1 << lg
    d = [e] * (2 * size)
    d[size : size + n] = values
    for i in range(size - 1, 0, -1):
      d[i] = op(d[2 * i], d[2 * i + 1])
    self.n, self.lg, self.size, self.d = n, lg, size, d
    self.lz = [id] * size
    self.op, self.e, self.mapping = op, e, mapping
    self.composition, self.id = composition, id

  def _put(self, k, f):
    self.d[k] = self.mapping(f, self.d[k])
    if k < self.size:
      self.lz[k] = self.composition(f, self.lz[k])

  def _push(self, k):
    f = self.lz[k]
    if f != self.id:
      self._put(2 * k, f)
      self._put(2 * k + 1, f)
      self.lz[k] = self.id

  def _pull(self, k):
    d = self.d
    d[k] = self.op(d[2 * k], d[2 * k + 1])

  def _push_ends(self, l, r):  # ancestors not inside [l, r)
    for i in range(self.lg, 0, -1):
      if l >> i << i != l:
        self._push(l >> i)
      if r >> i << i != r:
        self._push((r - 1) >> i)

  def _pull_ends(self, l, r):
    for i in range(1, self.lg + 1):
      if l >> i << i != l:
        self._pull(l >> i)
      if r >> i << i != r:
        self._pull((r - 1) >> i)

  def set(self, i, v):
    i += self.size
    self._push_ends(i, i + 1)
    self.d[i] = v
    self._pull_ends(i, i + 1)

  def get(self, i):
    i += self.size
    self._push_ends(i, i + 1)
    return self.d[i]

  def query_all(self):
    return self.d[1]

  def query(self, l, r):
    l += self.size
    r += self.size
    self._push_ends(l, r)
    d, op = self.d, self.op
    a = b = self.e
    while l < r:
      if l & 1:
        a = op(a, d[l])
        l += 1
      if r & 1:
        r -= 1
        b = op(d[r], b)
      l >>= 1
      r >>= 1
    return op(a, b)

  def apply(self, l, r, f):
    l += self.size
    r += self.size
    self._push_ends(l, r)
    a, b = l, r
    while a < b:
      if a & 1:
        self._put(a, f)
        a += 1
      if b & 1:
        b -= 1
        self._put(b, f)
      a >>= 1
      b >>= 1
    self._pull_ends(l, r)

  def max_right(self, l, pred):
    if l == self.n:
      return l
    d, op, size, acc = self.d, self.op, self.size, self.e
    l += size
    self._push_ends(l, l + 1)
    while True:
      while l % 2 == 0:
        l >>= 1
      if not pred(op(acc, d[l])):
        while l < size:
          self._push(l)
          l *= 2
          if pred(op(acc, d[l])):
            acc = op(acc, d[l])
            l += 1
        return l - size
      acc = op(acc, d[l])
      l += 1
      if (l & -l) == l:
        return self.n

  def min_left(self, r, pred):
    if r == 0:
      return 0
    d, op, size, acc = self.d, self.op, self.size, self.e
    r += size
    self._push_ends(r - 1, r)
    while True:
      r -= 1
      while r > 1 and r % 2:
        r >>= 1
      if not pred(op(d[r], acc)):
        while r < size:
          self._push(r)
          r = 2 * r + 1
          if pred(op(d[r], acc)):
            acc = op(d[r], acc)
            r -= 1
        return r + 1 - size
      acc = op(d[r], acc)
      if (r & -r) == r:
        return 0
