"""
Title: Knuth optimisation
Description: The interval DP
  $"dp"[l][r] = w(l, r) + min_(l < k < r) "dp"[l][k] + "dp"[k][r]$
  for all $0 <= l < r <= n$, with $"dp"[i][i + 1] = 0$, in
  $O(n^2)$ instead of $O(n^3)$: the smallest optimal split
  satisfies $"opt"[l][r - 1] <= "opt"[l][r] <= "opt"[l + 1][r]$.
Precondition: `w` satisfies the quadrangle inequality (DP
  techniques) and is monotone, $w(b, c) <= w(a, d)$ for
  $a <= b <= c <= d$. Sums of values over the
  interval, with values $>= 0$, satisfy both.
Time: $O(n^2)$ time and memory; $n = 2000$ takes about a second
  in PyPy.
Usage:
  # merge adjacent piles; a merge costs the sum of both
  dp = knuth(n, lambda l, r: pre[r] - pre[l])
  ans = dp[0][n]
Tweak: Leaf values: set `dp[i][i + 1]` before the loops; they
  do not change `opt`.
Example:
  >>> from itertools import accumulate
  >>> pre = [0, *accumulate([3, 1, 4, 1, 5])]  # adjacent piles
  >>> knuth(5, lambda l, r: pre[r] - pre[l])[0][5]
  32
  >>> # Optimal BST on keys searched f[i] times, cost = sum of
  >>> # f[i] * depth (root depth 1). Root k splits [l, k) and
  >>> # [k + 1, r): dp[l][r] covers keys l..r-2 (dp[l][l + 1]
  >>> # is the empty tree), so call knuth(n + 1, ...).
  >>> pre = [0, *accumulate([4, 2, 6, 3])]  # f
  >>> knuth(5, lambda l, r: pre[r - 1] - pre[l])[0][5]
  26
Notes: Signals: merging neighbours where a merge costs the
  total size (slimes, file merging, optimal alphabetic code),
  cutting a stick or string where a cut costs the current
  length (cut positions as boundaries), optimal BSTs. When
  unsure, test both conditions on all quadruples of small
  random inputs. If any two piles may merge (not only
  neighbours), it is Huffman instead: a greedy with a heap.
  Merging `K > 2` piles at once breaks the form; use an
  $O(n^3 / K)$ DP over `(l, r)` and the pile count modulo
  `K - 1`. The same `opt` bounds speed up layer DPs
  (`dp/dc-opt`).
Source: D. Knuth (1971), F. Yao (1980).
Status: stress-tested against the $O(n^3)$ DP
"""


def knuth(n, w):
  dp = [[0] * (n + 1) for _ in range(n + 1)]
  opt = [[0] * (n + 1) for _ in range(n + 1)]
  for i in range(n):
    opt[i][i + 1] = i
  for l in range(n - 2, -1, -1):
    dl, ol, hi = dp[l], opt[l], opt[l + 1]
    for r in range(l + 2, n + 1):
      lo = max(l + 1, ol[r - 1])
      best, bk = dl[lo] + dp[lo][r], lo
      for k in range(lo + 1, hi[r] + 1):
        v = dl[k] + dp[k][r]
        if v < best:
          best, bk = v, k
      ol[r] = bk
      dl[r] = best + w(l, r)
  return dp
