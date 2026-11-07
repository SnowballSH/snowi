"""
Title: Sparse table
Description: Static range minimum.
Precondition: 0-indexed, `l < r`. The operation is associative
  and idempotent, `op(x, x) = x`: min, max, gcd, bitwise and/or.
Time: $O(n log n)$ time and memory to build, $O(1)$ per query.
Usage:
  st = SparseTable(a)
  m = st.query(l, r)  # min of [l, r)
Tweak: Index of the minimum: build from `(a[i] << 20) + i` and
  take `& (2**20 - 1)`; ties go to the smallest index. Needs
  `n <= 2^20` and $|a_i| < 2^43$; otherwise use tuples `(a[i], i)`
  (packed big ints run twice as slow).
Example:
  >>> st = SparseTable([5, 2, 4, 2, 7, 1])
  >>> st.query(0, 3), st.query(2, 5), st.query(4, 5)
  (2, 2, 7)
  >>> a = [5, 2, 4, 2, 7, 1]  # leftmost min of [1, 5) is at 1
  >>> st = SparseTable([(x << 20) + i for i, x in enumerate(a)])
  >>> st.query(1, 5) & (2**20 - 1)
  1
Notes: Use it when the array never changes and there are many
  queries: $10^6$ range minimum queries, LCA as the shallowest
  vertex between two Euler tour positions, or the LCP of two
  suffixes from the LCP array. A non-idempotent operation (sum,
  product) double counts the overlap; use prefix sums instead.
  In PyPy, $n = 10^6$ builds in about 0.3 s, answers $10^6$
  queries in about 0.2 s and takes about 150 MB ($n log_2 n$ ints);
  under a tight memory limit use a segment tree instead.
Source: kactl RMQ (CC0).
Status: stress-tested against brute force
"""


class SparseTable:
  def __init__(self, values):
    t = [list(values)]
    k = 1
    while 2 * k <= len(values):
      prev = t[-1]
      m = len(prev) - k  # indexing: no slice copy, less memory
      t.append([min(prev[i], prev[i + k]) for i in range(m)])
      k *= 2
    self.t = t

  def query(self, l, r):
    j = (r - l).bit_length() - 1
    row = self.t[j]
    return min(row[l], row[r - (1 << j)])
