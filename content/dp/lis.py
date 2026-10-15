"""
Title: Longest increasing subsequence
Description: Indices of one longest strictly increasing
  subsequence of `a`, in increasing order. While it runs,
  `k + 1` is the length of the longest one ending at `i`.
Precondition: Values only need `<`. `a` may be empty.
Time: $O(n log n)$.
Usage:
  idx = lis(a)  # a[idx[0]] < a[idx[1]] < ...
  length = len(idx)
Tweak: Longest non-decreasing: `bisect_right` instead of
  `bisect_left`. Decreasing: negate the values.
Example:
  >>> idx = lis([3, 1, 4, 1, 5, 9, 2, 6])
  >>> idx, [[3, 1, 4, 1, 5, 9, 2, 6][i] for i in idx]
  ([1, 2, 4, 7], [1, 4, 5, 6])
  >>> lis([2, 2, 2]), lis([])
  ([2], [])
  >>> # Envelopes nested strictly in both width and height:
  >>> # equal widths sorted by height down never chain.
  >>> env = [(5, 4), (6, 4), (6, 7), (2, 3)]
  >>> env.sort(key=lambda e: (e[0], -e[1]))
  >>> len(lis([h for w, h in env]))
  3
Notes: `tail[k]` is the smallest last value of an increasing
  subsequence of length `k + 1` seen so far; it stays sorted.
  Cues: fewest deletions to make an array sorted is
  `n - LIS`; fewest changes to make it strictly increasing with
  any integers is `n` minus the longest non-decreasing
  subsequence of `a[i] - i` (if values must stay `>= 1`, first
  drop every `i` with `a[i] <= i`). Dilworth: the fewest
  non-increasing subsequences covering `a` equals the LIS length.
  Pairs ordered in two coordinates: sort by one, LIS on the
  other, with the tie order chosen as in the Example. `i` lies in
  some LIS iff $"end"[i] + "start"[i] - 1 = L$ (`start` from the
  reversed, negated list); in every LIS iff also no other such
  `j` has the same `end[j]`. Counting LIS: a Fenwick tree over
  values storing (length, count). Tuples work as values, but
  ints are several times faster in PyPy.
Source: kactl LIS (CC0).
Status: stress-tested against an $O(n^2)$ DP
"""

from bisect import bisect_left


def lis(a):
  tail, at = [], []  # at[k]: index of tail[k]
  prev = [-1] * len(a)
  for i, x in enumerate(a):
    k = bisect_left(tail, x)
    if k == len(tail):
      tail.append(x)
      at.append(i)
    else:
      tail[k] = x
      at[k] = i
    if k:
      prev[i] = at[k - 1]
  res, i = [], at[-1] if at else -1
  while i >= 0:
    res.append(i)
    i = prev[i]
  return res[::-1]
