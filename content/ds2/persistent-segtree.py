"""
Title: Persistent segment tree
Description: Point add and range sum over `[0, n)`; every `add`
  returns the root of a new version. Root 0 is the all-zero array.
  `kth(a, b, k)` reads version `b` minus version `a` as counts and
  returns the smallest `i` whose count over `[0, i]` exceeds `k`:
  the `k`-th smallest (0-indexed) value of a subarray when version
  `i` counts the compressed values of `a[:i]`.
Precondition: 0-indexed, half-open ranges, `n >= 1`. For `kth`,
  every count of `b` minus `a` is $>= 0$ and `0 <= k <` their
  total.
Time: $O(log n)$ per operation; $2 dot 10^5$ adds plus
  $2 dot 10^5$ `kth` at $n = 2 dot 10^5$ take about 1.5 s.
  Memory: each `add` creates at most $ceil(log_2 n) + 1$ nodes, 24
  bytes each in three flat lists.
Usage:
  st = PersistentSegTree(n)
  v1 = st.add(0, i, 5); v2 = st.add(v1, j, 3)
  s = st.query(v1, l, r)  # sum of [l, r) in version v1
  idx = st.kth(v1, v2, k)  # k-th unit of v2 - v1
Tweak: Point set: `add(v, i, x - query(v, i, i + 1))`. Max or
  min: set the new leaf, recompute the copied path from its
  children, and query by combining whole nodes (`add`'s delta and
  `query`'s prefix difference are only for sums).
Example:
  >>> from bisect import bisect_left
  >>> a = [5, 1, 4, 1, 3, 9]
  >>> vals = sorted(set(a))  # [1, 3, 4, 5, 9]
  >>> st = PersistentSegTree(len(vals))
  >>> root = [0]  # root[i]: counts of the values in a[:i]
  >>> for x in a:
  ...   root.append(st.add(root[-1], bisect_left(vals, x), 1))
  >>> l, r = 1, 5  # a[1:5] = [1, 4, 1, 3]
  >>> [vals[st.kth(root[l], root[r], k)] for k in range(4)]
  [1, 1, 3, 4]
  >>> c = bisect_left(vals, 4)  # how many x < 4 in a[1:5]
  >>> st.query(root[r], 0, c) - st.query(root[l], 0, c)
  3
  >>> v = st.add(root[6], 0, -2)  # versions stay intact
  >>> st.query(v, 0, 5), st.query(root[6], 0, 5)
  (4, 6)
Notes: Think of it for "k-th smallest in a range", "count values
  in `[x, y)` in a subarray" online, "answer for the array as it
  was after update `t`", or when a sweep must be queried at any
  earlier moment (versions over time or over a sorted key). The
  trick is prefix versions over the value axis: version `r` minus
  version `l` is the histogram of `a[l:r]`, and one walk down both
  trees finds the k-th unit. On a tree, a version per vertex built
  from its parent gives path histograms as
  $v + u - "lca" - "parent"("lca")$; walk four roots at once.
  Memory is the usual limit: $2 dot 10^5$ updates on
  $n = 2 dot 10^5$ make about $3.8 dot 10^6$ nodes. On a static
  array a wavelet matrix answers the same k-th smallest and count
  queries online with about half the memory; offline, parallel
  binary search also works.
Source: Folklore ("chairman tree"); path copying.
Status: stress-tested against brute force over all versions
"""


class PersistentSegTree:
  def __init__(self, n):
    self.n = n
    self.L, self.R, self.S = [0], [0], [0]  # node 0: empty tree

  def add(self, x, i, delta):  # returns the new root
    L, R, S = self.L, self.R, self.S
    root, lo, hi = len(S), 0, self.n
    while True:
      y = len(S)
      L.append(L[x])
      R.append(R[x])
      S.append(S[x] + delta)
      if hi - lo == 1:
        return root
      mid = (lo + hi) >> 1
      if i < mid:
        L[y], x, hi = y + 1, L[x], mid
      else:
        R[y], x, lo = y + 1, R[x], mid

  def prefix(self, x, r):  # sum of [0, r) in version x
    L, R, S = self.L, self.R, self.S
    lo, hi, s = 0, self.n, 0
    while x and lo < r:
      if hi <= r:
        return s + S[x]
      mid = (lo + hi) >> 1
      if r <= mid:
        x, hi = L[x], mid
      else:
        s += S[L[x]]
        x, lo = R[x], mid
    return s

  def query(self, x, l, r):
    return self.prefix(x, r) - self.prefix(x, l)

  def kth(self, a, b, k):  # k-th unit of b minus a
    L, R, S = self.L, self.R, self.S
    lo, hi = 0, self.n
    while hi - lo > 1:
      mid = (lo + hi) >> 1
      c = S[L[b]] - S[L[a]]
      if k < c:
        a, b, hi = L[a], L[b], mid
      else:
        k -= c
        a, b, lo = R[a], R[b], mid
    return lo
