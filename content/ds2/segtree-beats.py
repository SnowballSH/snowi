"""
Title: Segment tree beats
Description: Range chmin ($a_i = min(a_i, x)$), range chmax
  ($a_i = max(a_i, x)$) and range add, with range sum, max and
  min. A node keeps its max, strict second max and count of the
  max, and the same for the min. Chmin stops at a node where only
  the max class changes and fixes the sum from the count.
Precondition: 0-indexed, half-open ranges. Values, `x` and range
  sums stay below $4 dot 10^18$ in absolute value.
Time: $O((n + q) log^2 n)$ amortized with range add, and
  $O((n + q) log n)$ with only chmin, chmax and queries. In PyPy
  (fresh process) $10^5$ random chmin and sum operations at
  $n = 10^5$ take 1.4 s; random mixes of all four kinds take 2.1 s
  for $3 dot 10^4$, 2.8 s for $5 dot 10^4$ and 9 s for
  $2 dot 10^5$ (C++: 0.06 s to 0.6 s). The JIT setting in Tweak
  cuts the mixes to 1.1 s, 1.8 s and 7.8 s; beyond $5 dot 10^4$
  mixed operations use C++.
Usage:
  st = Beats(a)
  st.chmin(l, r, x); st.chmax(l, r, x); st.add(l, r, x)
  s, mx, mn = st.query(l, r)  # sum, max, min of [l, r)
  # an empty range gives (0, -Beats.INF, Beats.INF)
Tweak: Assign `x` on `[l, r)`: `chmin(l, r, x)`, `chmax(l, r,
  x)`. For only chmin and sum, drop the min lists and `_max`.
  More node data goes into new lists, `_pull`, and every `_add`,
  `_min`, `_max` that changes it, in $O(1)$ per node. Speed:
  `import pypyjit` and
  `pypyjit.set_param("max_unroll_recursion=-1")` at the start of
  the program shorten the JIT warm-up.
Example:
  >>> st = Beats([5, 1, 4, 1, 3])
  >>> st.chmin(0, 5, 3); st.add(1, 4, 2); st.chmax(0, 2, 4)
  >>> st.query(0, 5), st.query(2, 3)  # 4 4 5 3 3
  ((19, 5, 3), (5, 5, 5))
  >>> st.chmin(2, 4, 0); st.chmax(2, 4, 0)  # assign 0
  >>> st.query(0, 5)  # 4 4 0 0 3
  (11, 4, 0)
Notes: Think of it when a range update clamps values ("set all
  above `x` to `x`", "everybody gets at least `x`") and sums or
  counts are asked, which a lazy tag cannot do. Without adds, the
  bound charges each extra recursion to the number of distinct
  values in a node: it drops when the max class merges with the
  second, and a chmin or chmax raises it by at most one in each of
  the $O(log n)$ partly covered nodes. One range add can raise it
  by half a node's length, so with adds Ji's proof needs a finer
  potential (over "tag classes") and gets $O(log^2 n)$; random
  tests visit about $4 log_2 n$ nodes per operation, like a plain
  segment tree. The same "recurse while something
  changes" walk handles `a_i %= x`, `a_i /= x`,
  `a_i = floor(sqrt(a_i))` and `a_i &= x` on values $>= 0$: stop
  where nothing changes (`mx < x` for `%=`) and update leaves;
  each value shrinks $O(log V)$ times. Pitfalls: the second max is
  strict and $-infinity$ at a leaf, adds must not shift the
  sentinels, and a node with one or two distinct values links the
  max and min fields (`_min`, `_max`). The recursion runs on an
  explicit stack; nodes to recompute are pulled in reverse visit
  order, which puts every child before its parent.
Source: Ji Ruyi, segment tree beats (2016).
Status: stress-tested against brute force
"""


class Beats:
  INF = (1 << 63) - 1  # sentinel, never shifted

  def __init__(self, a):
    n, INF = len(a), self.INF
    size = 1
    while size < n:
      size *= 2
    m = 2 * size  # top-down tree on [0, n): children 2k, 2k + 1
    self.n = n
    self.S, self.LZ, self.LEN = [0] * m, [0] * m, [0] * m
    self.MX, self.MX2 = [-INF] * m, [-INF] * m
    self.MN, self.MN2 = [INF] * m, [INF] * m
    self.CMX, self.CMN = [0] * m, [0] * m
    order, stack = [], [1, 0, n] if n else []
    while stack:
      hi, lo, k = stack.pop(), stack.pop(), stack.pop()
      if hi - lo == 1:
        self.S[k] = self.MX[k] = self.MN[k] = a[lo]
        self.CMX[k] = self.CMN[k] = self.LEN[k] = 1
        continue
      order.append(k)
      mid = (lo + hi) >> 1
      stack += (2 * k, lo, mid, 2 * k + 1, mid, hi)
    for k in reversed(order):
      self._pull(k)

  def _pull(self, k):
    a, b = 2 * k, 2 * k + 1
    S, MX, MX2, CMX = self.S, self.MX, self.MX2, self.CMX
    MN, MN2, CMN = self.MN, self.MN2, self.CMN
    S[k] = S[a] + S[b]
    self.LEN[k] = self.LEN[a] + self.LEN[b]
    x, y = MX[a], MX[b]
    if x == y:
      MX[k], CMX[k] = x, CMX[a] + CMX[b]
      MX2[k] = max(MX2[a], MX2[b])
    elif x > y:
      MX[k], CMX[k], MX2[k] = x, CMX[a], max(MX2[a], y)
    else:
      MX[k], CMX[k], MX2[k] = y, CMX[b], max(x, MX2[b])
    x, y = MN[a], MN[b]
    if x == y:
      MN[k], CMN[k] = x, CMN[a] + CMN[b]
      MN2[k] = min(MN2[a], MN2[b])
    elif x < y:
      MN[k], CMN[k], MN2[k] = x, CMN[a], min(MN2[a], y)
    else:
      MN[k], CMN[k], MN2[k] = y, CMN[b], min(x, MN2[b])

  def _add(self, k, x):
    INF = self.INF
    self.S[k] += x * self.LEN[k]
    self.MX[k] += x
    self.MN[k] += x
    self.LZ[k] += x
    if self.MX2[k] != -INF:
      self.MX2[k] += x
    if self.MN2[k] != INF:
      self.MN2[k] += x

  def _min(self, k, x):  # needs MX2[k] < x < MX[k]
    MX, MN = self.MX, self.MN
    self.S[k] += (x - MX[k]) * self.CMX[k]
    if MN[k] == MX[k]:
      MN[k] = x
    elif self.MN2[k] == MX[k]:
      self.MN2[k] = x
    MX[k] = x

  def _max(self, k, x):  # needs MN[k] < x < MN2[k]
    MX, MN = self.MX, self.MN
    self.S[k] += (x - MN[k]) * self.CMN[k]
    if MX[k] == MN[k]:
      MX[k] = x
    elif self.MX2[k] == MN[k]:
      self.MX2[k] = x
    MN[k] = x

  def _push(self, k):
    MX, MN, LZ = self.MX, self.MN, self.LZ
    for c in (2 * k, 2 * k + 1):
      if LZ[k]:
        self._add(c, LZ[k])
      if MX[c] > MX[k]:
        self._min(c, MX[k])
      if MN[c] < MN[k]:
        self._max(c, MN[k])
    LZ[k] = 0

  def _update(self, l, r, x, kind):  # 0 chmin, 1 chmax, 2 add
    MX, MX2, MN, MN2 = self.MX, self.MX2, self.MN, self.MN2
    stack, pulls = [1, 0, self.n], []
    while stack:
      hi, lo, k = stack.pop(), stack.pop(), stack.pop()
      if r <= lo or hi <= l:
        continue
      inside = l <= lo and hi <= r
      if kind == 0:
        if MX[k] <= x:
          continue
        if inside and MX2[k] < x:
          self._min(k, x)
          continue
      elif kind == 1:
        if MN[k] >= x:
          continue
        if inside and MN2[k] > x:
          self._max(k, x)
          continue
      elif inside:
        self._add(k, x)
        continue
      self._push(k)
      pulls.append(k)
      mid = (lo + hi) >> 1
      stack += (2 * k, lo, mid, 2 * k + 1, mid, hi)
    for k in reversed(pulls):
      self._pull(k)

  def chmin(self, l, r, x):
    self._update(l, r, x, 0)

  def chmax(self, l, r, x):
    self._update(l, r, x, 1)

  def add(self, l, r, x):
    self._update(l, r, x, 2)

  def query(self, l, r):  # (sum, max, min) of [l, r)
    S, MX, MN = self.S, self.MX, self.MN
    s, mx, mn = 0, -self.INF, self.INF
    stack = [1, 0, self.n]
    while stack:
      hi, lo, k = stack.pop(), stack.pop(), stack.pop()
      if r <= lo or hi <= l:
        continue
      if l <= lo and hi <= r:
        s += S[k]
        mx = max(mx, MX[k])
        mn = min(mn, MN[k])
        continue
      self._push(k)
      mid = (lo + hi) >> 1
      stack += (2 * k, lo, mid, 2 * k + 1, mid, hi)
    return s, mx, mn
