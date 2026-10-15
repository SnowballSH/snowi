"""
Title: Divide and conquer optimisation
Description: One layer of a DP,
  $"new"[i] = min_(j < i) "old"[j] + "cost"(j, i)$ for every
  `i`, when the smallest optimal `j` is non-decreasing in `i`.
  It solves the middle `i` by scanning, then each half with the
  `j` range cut at that optimum (explicit stack, no recursion).
Precondition: `cost` satisfies the quadrangle inequality
  $c(a, c) + c(b, d) <= c(a, d) + c(b, c)$ for
  $a <= b <= c <= d$, or the smallest optimal `j` is known to
  be non-decreasing in `i` (a monotone choice among tied `j` is
  not enough). Unreachable states are `1 << 62`, also used for
  `new[0]` and as a cap. Keep values below $2^63$: big ints are
  slow.
Time: $O(n log n)$ calls to `cost` per layer; about
  $10^7$ calls per second in PyPy.
Usage:
  # dp[i]: best split of the first i elements; cost(j, i)
  # is the cost of one group [j, i).
  dp = [1 << 62] * (n + 1)
  dp[0] = 0
  for _ in range(K):
    dp = dc_layer(dp, cost)  # dp[n]: K groups
Tweak: Maximum: negate `old` and `cost` (then the reverse
  inequality must hold). Only some `i` needed, or `j <= i`
  allowed: change the first stack entry and the scan bound `m`.
Example:
  >>> a = [4, 1, 3, 2, 5]  # K groups, sum of squared group sums
  >>> pre = [0]
  >>> for x in a: pre.append(pre[-1] + x)
  >>> cost = lambda j, i: (pre[i] - pre[j]) ** 2
  >>> dp = [1 << 62] * 6
  >>> dp[0] = 0
  >>> for K in range(1, 4):
  ...   dp = dc_layer(dp, cost)
  ...   print(K, dp[5])
  1 225
  2 113
  3 75
Notes: The quadrangle inequality says crossing intervals cost
  at most nested ones. It holds for $f(S_i - S_j)$ with $f$
  convex and $S$ non-decreasing (squared sums), for "pairs
  inside the group" counts such as inversions or equal pairs,
  and for sums of such costs. When unsure, brute-force `opt` on
  small cases and check that it never decreases. Costs that are
  not $O(1)$ (equal pairs in `[j, i)` for a minimum, distinct
  values for a maximum): keep a window with add/remove like Mo's
  algorithm and move it to `[j, i)` inside `cost`; the scans move
  it $O(n log n)$ per layer. If the cost has the quadrangle
  inequality, the answer for exactly `K` groups is convex in
  `K`, so the Aliens trick (`dp/wqs`) with a fast inner DP avoids
  the factor `K`. This does not apply to
  $"dp"[i] = min_(j < i) "dp"[j] + "cost"(j, i)$ on a single
  array, where `dp[j]` is still unknown: use CDQ divide and
  conquer on `i`, or a deque of candidates with binary search.
Source: Classic; cp-algorithms "Divide and Conquer DP".
Status: stress-tested against an $O(n^2)$ DP
"""


def dc_layer(old, cost):
  n = len(old)
  res = [1 << 62] * n
  # Solve i in [l, r) knowing the optimal j lies in [lo, hi].
  stack = [1, n, 0, n - 1]
  while stack:
    hi = stack.pop()
    lo = stack.pop()
    r = stack.pop()
    l = stack.pop()
    if l >= r:
      continue
    m = (l + r) >> 1
    best, val = lo, old[lo] + cost(lo, m)
    for j in range(lo + 1, min(m, hi + 1)):
      v = old[j] + cost(j, m)
      if v < val:
        best, val = j, v
    if val < res[m]:
      res[m] = val
    stack += (l, m, lo, best, m + 1, r, best, hi)
  return res
