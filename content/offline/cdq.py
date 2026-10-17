"""
Title: CDQ divide and conquer
Description: 3D dominance counting: `cnt[i]` is the number of
  points `j != i` with $a_j <= a_i$, $b_j <= b_i$ and
  $c_j <= c_i$; equal points count each other.
Precondition: Points are tuples (or lists) of three comparable
  values; `c` is compressed inside.
Time: $O(n log^2 n)$, recursion depth $log_2 n$. In PyPy,
  $n = 10^5$ takes about 1.3 s and $2 dot 10^5$ about 2.8 s.
Usage:
  cnt = dominance(pts)  # pts: list of (a, b, c)
Tweak: Strict `<` in `b`: use `<` in the sweep; in `c`:
  `query(c[x])`. Then equal points never count each other, so drop
  the final loop. Strict in `a` only: make it `b` or `c`. Weights:
  add `wt[y]` and undo `wt[idx[j]]` instead of 1; copies of a point
  get `cnt[last] + wt[last] - wt[i]`. DP over the order (`j` must be
  final before it feeds `i`): recurse left, do the cross step, then
  recurse right, sorting copies of both halves by `b` instead of
  merging.
Example:
  >>> dominance([(1, 3, 2), (2, 1, 1), (3, 2, 3), (4, 4, 2)])
  [0, 0, 1, 2]
  >>> dominance([(1, 1, 1), (2, 2, 2), (1, 1, 1), (0, 5, 0)])
  [1, 2, 1, 0]

  # CDQ on time: op t inserts (X[t], Y[t]) or, if ask[t], counts
  # inserted points <= (X[t], Y[t]). Keep idx in time order (time
  # is a), merge by Y, Fenwick over X (compressed). Only inserts
  # of [l, m) add (and are undone), only questions of [m, r)
  # read; no copy fix.
  for i in range(m, r):
    t = idx[i]
    while k < m and Y[idx[k]] <= Y[t]:
      if not ask[idx[k]]:
        add(X[idx[k]], 1)
      merged.append(idx[k])
      k += 1
    if ask[t]:
      res[t] += query(X[t] + 1)
    merged.append(t)
Notes: Think of it for "count pairs `j < i` with two more
  conditions", a 3D partial order, interleaved insertions and
  rectangle queries given offline (split a rectangle into four
  prefix corners with signs $plus.minus 1$), or
  `dp[i] = 1 + max dp[j]` over `j < i` with $x_j <= x_i$ and
  $y_j <= y_i$ (prefix-max Fenwick, reset touched cells by
  assignment). Undo the Fenwick inserts; clearing it per call costs
  $O(n)$ each time. Two conditions only: sort and sweep, no CDQ.
  Four: CDQ inside CDQ, $O(n log^3 n)$, or bitsets ($n^2 / 64$).
  Other cross steps work too: a merge of two convex hulls, a
  two-pointer count, a convolution. A static rectangle count needs
  no CDQ: sweep $x$ with a Fenwick tree over $y$.
Uses: ds/fenwick
Source: Chen Danqi (CDQ); UCup S4-S5 Training Packet 14.
Status: stress-tested against brute force
"""


def dominance(p):
  n = len(p)
  order = sorted(range(n), key=p.__getitem__)
  rank = {v: i for i, v in enumerate(sorted({x[2] for x in p}))}
  b, c = [x[1] for x in p], [rank[x[2]] for x in p]
  cnt, idx, fw = [0] * n, order[:], Fenwick(n)
  add, query = fw.add, fw.query

  def rec(l, r):
    if r - l < 2:
      return
    m = (l + r) // 2
    rec(l, m)
    rec(m, r)  # halves now sorted by b
    k, merged = l, []
    for i in range(m, r):  # cross step: [l, m) into [m, r)
      x = idx[i]
      while k < m and b[idx[k]] <= b[x]:
        y = idx[k]
        add(c[y], 1)
        merged.append(y)
        k += 1
      cnt[x] += query(c[x] + 1)
      merged.append(x)
    for j in range(l, k):
      add(c[idx[j]], -1)
    idx[l:r] = merged + idx[k:m]

  rec(0, n)
  for i in range(n - 2, -1, -1):  # copies count as the last
    if p[order[i]] == p[order[i + 1]]:
      cnt[order[i]] = cnt[order[i + 1]]
  return cnt
