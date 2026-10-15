"""
Title: Wavelet matrix
Description: Static array of non-negative integers answering order
  statistics on any subarray: `kth(l, r, k)` is the `k`-th
  smallest value (0-indexed) of `a[l:r]`;
  `count_less(l, r, x)` counts its values `< x`,
  `count(l, r, lo, hi)` its values in `[lo, hi)`, and
  `rank(x, r)` the copies of `x` in `a[:r]`. Level `h` stores a
  prefix count of the ones in bit `h` of every value, then moves
  the values whose bit `h` is 0 to the front (stably): that is the
  order level `h - 1` sees.
Precondition: 0-indexed, half-open ranges, `a[i] >= 0`. `kth`
  needs `0 <= k < r - l`; any integer `x` is fine elsewhere.
Time: $O(n log V)$ build and $O(log V)$ per query, where $V$ is
  the largest value plus one. At $n = 2 dot 10^5$, $V = 10^9$:
  build 0.3 s, $2 dot 10^5$ `kth` about 1 s (random memory access
  dominates). Memory: $n log V$ list slots of 8 bytes.
Usage:
  wm = WaveletMatrix(a)
  v = wm.kth(l, r, k)  # k-th smallest of a[l:r]
  c = wm.count_less(l, r, x); d = wm.count(l, r, lo, hi)
  copies = wm.rank(x, r)  # of x in a[:r]
Tweak: Negative or huge values: compress to indices in `vals`
  first (this also makes $log V = log n$ and saves memory), then
  map `kth` back through `vals`. Sum of the `k` smallest in
  `a[l:r]`: per level `h`, also keep prefix sums of the values in
  their order after level `h`'s partition (zeros first). In `kth`,
  add the sum over `[l - l1, r - r1)` whenever the walk goes
  right; at the end add `k * res`, the copies of the result among
  the `k` smallest. A prefix count at every position is faster in
  PyPy than packed words with a popcount table (no
  `int.bit_count` before Python 3.10), at 64 times the memory of
  the C++ bit vectors.
Example:
  >>> wm = WaveletMatrix([5, 1, 4, 1, 3, 9])
  >>> wm.kth(1, 5, 2)  # [1, 4, 1, 3] sorted is 1 1 3 4
  3
  >>> wm.count_less(1, 5, 4), wm.count(0, 6, 3, 6), wm.rank(1, 4)
  (3, 3, 2)
  >>> l, r, x = 0, 6, 5  # largest value < x in a[l:r], if c > 0
  >>> c = wm.count_less(l, r, x); wm.kth(l, r, c - 1)
  4
  >>> wm.kth(l, r, c)  # smallest value >= x, if c < r - l
  5
  >>> from bisect import bisect_left
  >>> a = [-7, 10**18, 3, -7]
  >>> vals = sorted(set(a))  # compressed: [-7, 3, 10**18]
  >>> wc = WaveletMatrix([bisect_left(vals, x) for x in a])
  >>> z = bisect_left(vals, 0)  # values < 0 have indices < z
  >>> vals[wc.kth(0, 4, 3)], wc.count_less(0, 4, z)
  (1000000000000000000, 2)
Notes: Think of it for static "k-th smallest in a range", "median
  of a subarray", "how many values in `[x, y)` lie in `a[l:r]`"
  and "closest value to `x` in a range", online and without
  updates. Against a persistent segment tree in PyPy it needs
  about half the memory and builds 3 to 4 times faster; `kth`
  takes about as long, or 1.5 times less on compressed values.
  Level `h` sends the values whose bit `h` is 0 to the front and
  the others behind all `zeros[h]` zeros, so a range maps to one
  range on either side by counting ones, and a query walks one
  root-to-leaf path. The order at level `h` is not sorted by the
  higher bits as a number (bit `h + 1` is the main key), so data
  kept per level, like the sums in Tweak, must follow the
  partitions. On a tree, build it on the Euler tour for subtree
  queries. With point updates use offline methods (parallel
  binary search, sorted sweeps) or sqrt decomposition instead.
Source: Claude, Navarro, Ordonez, "The wavelet matrix" (2012).
Status: stress-tested against brute force
"""


class WaveletMatrix:
  def __init__(self, a):
    n = self.n = len(a)
    lg = self.lg = max(a, default=0).bit_length()
    # ones[h * (n + 1) + i]: ones in the first i bits of level h
    self.ones = ones = [0] * (lg * (n + 1))
    self.zeros = [0] * lg
    for h in range(lg - 1, -1, -1):
      base, c = h * (n + 1), 0
      for i, x in enumerate(a):
        ones[base + i] = c
        c += x >> h & 1
      ones[base + n] = c
      self.zeros[h] = n - c
      low = [x for x in a if not x >> h & 1]
      a = low + [x for x in a if x >> h & 1]

  def kth(self, l, r, k):
    ones, zeros, m = self.ones, self.zeros, self.n + 1
    res = 0
    for h in range(self.lg - 1, -1, -1):
      l1, r1 = ones[h * m + l], ones[h * m + r]
      z = r - l - (r1 - l1)
      if k < z:
        l -= l1
        r -= r1
      else:
        k -= z
        res |= 1 << h
        l, r = zeros[h] + l1, zeros[h] + r1
    return res

  def count_less(self, l, r, x):
    if x >> self.lg:
      return 0 if x < 0 else r - l
    ones, zeros, m = self.ones, self.zeros, self.n + 1
    res = 0
    for h in range(self.lg - 1, -1, -1):
      l1, r1 = ones[h * m + l], ones[h * m + r]
      if x >> h & 1:
        res += r - l - (r1 - l1)
        l, r = zeros[h] + l1, zeros[h] + r1
      else:
        l -= l1
        r -= r1
    return res

  def count(self, l, r, lo, hi):  # 0 if lo >= hi
    c = self.count_less(l, r, hi) - self.count_less(l, r, lo)
    return max(c, 0)

  def rank(self, x, r):
    return self.count(0, r, x, x + 1)
