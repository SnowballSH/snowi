"""
Title: Fenwick tree
Description: Point add and prefix sums. `lower_bound(s)` returns
  the smallest `i` whose prefix sum over `[0, i]` is at least `s`:
  `n` if none is, `-1` if `s <= 0`.
Precondition: 0-indexed. `lower_bound` needs all values `>= 0`.
Time: $O(log n)$ per operation.
Usage:
  fw = Fenwick(n)
  fw.add(i, 5)
  s = fw.query(r) - fw.query(l)  # sum of [l, r)
Example:
  >>> fw = Fenwick(4)
  >>> for i, x in enumerate([3, 0, 2, 5]): fw.add(i, x)
  >>> fw.query(3), fw.lower_bound(4), fw.lower_bound(11)
  (5, 2, 4)
Notes: Build in $O(n)$: set `t[i] = a[i]`, then for each `i` add
  `t[i]` to `t[i | (i + 1)]` when that index is `< n`.
Source: kactl FenwickTree (CC0).
Status: stress-tested against a list
"""


class Fenwick:
  def __init__(self, n):
    self.n = n
    self.t = [0] * n

  def add(self, i, delta):
    t, n = self.t, self.n
    while i < n:
      t[i] += delta
      i |= i + 1

  def query(self, r):  # sum of [0, r)
    t, s = self.t, 0
    while r > 0:
      s += t[r - 1]
      r &= r - 1
    return s

  def lower_bound(self, s):
    if s <= 0:
      return -1
    t, pos, pw = self.t, 0, 1 << self.n.bit_length()
    while pw:
      if pos + pw <= self.n and t[pos + pw - 1] < s:
        pos += pw
        s -= t[pos - 1]
      pw >>= 1
    return pos
