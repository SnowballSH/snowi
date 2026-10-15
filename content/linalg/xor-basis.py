"""
Title: XOR basis
Description: Linear basis over $FF_2$ of the inserted numbers: the
  set of all subset XORs (the span) is exactly the set of XORs of
  basis vectors. `insert` returns whether `x` was independent,
  `max_xor(x)` / `min_xor(x)` give the largest / smallest `x ^ v`
  over `v` in the span, and `kth(k)` the `k`-th smallest span
  element. `PrefixXorBasis` answers "max subset XOR of
  `a[l..r]`" while `a` is inserted left to right.
Precondition: $0 <= x < 2^"bits"$. `kth`: $0 <= k < 2^"rank"$
  (`k = 0` gives 0). `PrefixXorBasis`: insert `a[r]` with `i = r`
  for `r = 0, 1, ...`; `max_xor(l)` then covers `a[l..r]`,
  inclusive.
Time: $O("bits")$ per operation; $2 dot 10^5$ inserts and queries
  of 60-bit values take about 0.5 s in PyPy.
Usage:
  xb = XorBasis()
  indep = xb.insert(x)
  best, low = xb.max_xor(), xb.min_xor(y)
  inside = xb.contains(y)
  v = xb.kth(k)  # xb.rank: basis size
  pb = PrefixXorBasis()
  for r in range(n): pb.insert(a[r], r)  # then pb.max_xor(l)
Tweak: `XorBasis(30)` for values below $2^30$ is faster. More than
  63 bits works (big ints) but is slower.
Example:
  >>> xb = XorBasis(4)
  >>> [xb.insert(x) for x in (6, 5, 3, 9)]  # 3 = 6 ^ 5
  [True, True, False, True]
  >>> xb.rank, xb.max_xor(), xb.min_xor(4), xb.contains(4)
  (3, 15, 1, False)
  >>> [xb.kth(k) for k in range(8)]  # the whole span, sorted
  [0, 3, 5, 6, 9, 10, 12, 15]
  >>> 2 ** (4 - xb.rank)  # subsets of (6, 5, 3, 9) with XOR 15
  2
  >>> pb = PrefixXorBasis(3)
  >>> for r, x in enumerate([5, 6, 3, 4]): pb.insert(x, r)
  >>> [pb.max_xor(l) for l in range(4)]  # a[l..3]
  [7, 7, 7, 4]
Notes: `n` numbers of rank $r$: $2^r$ distinct subset XORs, and
  each one (0 included, via the empty set) is the XOR of exactly
  $2^(n - r)$ subsets. Non-empty subsets with XOR 0 exist iff
  $n > r$; the smallest non-zero span element is `kth(1)`. The
  nonzero slots `p[b]` with `pos[b] >= l` of a `PrefixXorBasis`
  form a basis of `a[l..r]` (empty slots have `pos[b] = 0` too),
  so their count is its rank. For any other query `insert` them
  into an `XorBasis` (copying `p` would leave `rank` at 0). As `l`
  decreases the span grows only at the `pos[b]`.
  Cues: "maximum XOR of a subset", "can x be formed", "number of
  distinct XORs", "XOR of a path where edges may repeat" (insert
  `d[u] ^ d[v] ^ w` for non-tree edges, answer
  `min_xor(d[s] ^ d[t])`), "XOR of chosen segments" (basis of
  prefix XORs). See the linear algebra facts.
Source: Folklore linear basis; prefix (timestamp) basis from
  Codeforces 1100F.
Status: stress-tested against span enumeration
"""


class XorBasis:
  def __init__(self, bits=60):
    self.p = [0] * bits  # p[b]: vector with top bit b, or 0
    self.rank = 0

  def insert(self, x):
    p = self.p
    while x:
      b = x.bit_length() - 1
      if not p[b]:
        p[b] = x
        self.rank += 1
        return True
      x ^= p[b]
    return False

  def max_xor(self, x=0):
    for v in reversed(self.p):
      if x ^ v > x:
        x ^= v
    return x

  def min_xor(self, x):
    for v in reversed(self.p):
      if x ^ v < x:
        x ^= v
    return x

  def contains(self, x):
    return self.min_xor(x) == 0

  def kth(self, k):  # pivot bits of the answer spell out k
    p, x, i = self.p, 0, self.rank
    for b in range(len(p) - 1, -1, -1):
      if p[b]:
        i -= 1
        if (x >> b ^ k >> i) & 1:
          x ^= p[b]
    return x


class PrefixXorBasis:
  def __init__(self, bits=60):
    self.p = [0] * bits
    self.pos = [0] * bits  # each slot keeps the newest vector

  def insert(self, x, i):
    p, pos = self.p, self.pos
    while x:
      b = x.bit_length() - 1
      if not p[b] or pos[b] < i:
        p[b], x = x, p[b]
        pos[b], i = i, pos[b]
        if not x:  # took an empty slot
          return
      x ^= p[b]

  def max_xor(self, l):
    p, pos, x = self.p, self.pos, 0
    for b in range(len(p) - 1, -1, -1):
      if pos[b] >= l and x ^ p[b] > x:
        x ^= p[b]
    return x
