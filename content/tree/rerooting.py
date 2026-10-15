"""
Title: Rerooting
Description: A tree DP for every root at once. Rooted anywhere, the
  value of the subtree of `v` is `fin(acc, v)`, where `acc` merges
  `lift(x, c, v)` over the children `c` of `v` (`x` is the value of
  the subtree of `c`), starting from `e`. `dp[v]` is the value of
  the whole tree rooted at `v`. The values on directed edges come
  too: with `par` from root 0 and `p = par[v]`, `dn[v]` is
  `lift(x, v, p)` for `v`'s side of the edge and `up[v]` is
  `lift(y, p, v)` for `p`'s side (`e` at vertex 0).
Precondition: `adj` holds both directions of the $n - 1$ edges of
  a tree on vertices `0..n-1`, $n >= 1$. `merge` is associative
  and commutative with identity `e`; it needs no inverse.
Time: $O(n)$ calls of `merge`, `lift` and `fin`; iterative. PyPy,
  $n = 2 dot 10^5$: about 0.3 s with int values, 0.5 s with tuples.
Usage:
  # eccentricity: the largest distance from each vertex
  rr = Reroot(adj, 0, max, lambda x, u, v: x + 1,
              lambda acc, v: acc)
  rr.dp[v]  # also rr.par, rr.dn, rr.up
Tweak: `lift(x, u, v)` carries `u`'s side across the edge to `v`
  and may depend on the direction (edges that must be reversed).
  Edge weights: let `adj` hold `(v, w)` pairs, keep `pw[v] = w`
  (the edge to `par[v]`) in the BFS, and pass `pw[v]` and `pw[c]`
  to the two `lift` calls. For speed, inline the three functions.
Example:
  >>> adj = [[1], [0, 2, 3], [1, 4], [1], [2]]  # 0-1 1-2 1-3 2-4
  >>> rr = Reroot(adj, 0, max, lambda x, u, v: x + 1,
  ...             lambda acc, v: acc)
  >>> rr.dp, rr.par  # eccentricities; parents from root 0
  ([3, 2, 2, 3, 3], [-1, 0, 1, 1, 2])
  >>> # Sum of distances from every vertex; values (size, sum).
  >>> merge = lambda a, b: (a[0] + b[0], a[1] + b[1])
  >>> lift = lambda x, u, v: (x[0], x[1] + x[0])  # one edge on
  >>> fin = lambda x, v: (x[0] + 1, x[1])
  >>> rr = Reroot(adj, (0, 0), merge, lift, fin)
  >>> [s for _, s in rr.dp]
  [8, 5, 6, 8, 9]
  >>> rr.dn[2], rr.up[2]  # edge 1-2: {2, 4} and {1, 0, 3} sides
  ((2, 3), (3, 5))
  >>> # Connected vertex sets containing v, mod a prime: product
  >>> # over neighbours of (1 + count on their side); no division.
  >>> MOD = 998244353
  >>> cs = Reroot(adj, 1, lambda a, b: a * b % MOD,
  ...             lambda x, u, v: (x + 1) % MOD, lambda x, v: x)
  >>> cs.dp
  [7, 12, 10, 7, 6]
Notes: Cues: "for every vertex as the root", "sum of distances (or
  farthest vertex) from each vertex", "the best vertex to root the
  tree at", answers on both sides of every edge (`dn`, `up`: a
  value per directed edge). Write the rooted DP as `fin` of merged
  `lift`s first and check it against a single root; rerooting then
  needs no new ideas. When `merge` has an inverse (sums), the
  direct formula is shorter: sum of distances is `ans[c] = ans[p]
  + n - 2 size[c]`. Each vertex costs $O(deg)$ merges, so a heavy
  value (a list of size $s$) costs $O(n s)$; keep it small or
  shift indices instead of copying.
Source: Folklore; prefix and suffix merges as in maspy's library
  (CC0).
Status: stress-tested against a DP from every root
"""


class Reroot:
  def __init__(self, adj, e, merge, lift, fin):
    n = len(adj)
    par, q, lo = [-1] * n, [0], [n] * (n + 1)
    for i in range(n):  # BFS; kids of q[i]: q[lo[i]:lo[i + 1]]
      u = q[i]
      lo[i] = len(q)
      for v in adj[u]:
        if v != par[u]:
          par[v] = u
          q.append(v)
    dp, dn, up = [e] * n, [e] * n, [e] * n
    for i in range(n - 1, 0, -1):
      v, acc = q[i], e
      for j in range(lo[i], lo[i + 1]):
        acc = merge(acc, dn[q[j]])
      dn[v] = lift(fin(acc, v), v, par[v])
    for i in range(n):
      v = q[i]
      acc = up[v]  # prefix merges, then suffix merges
      for j in range(lo[i], lo[i + 1]):
        c = q[j]
        up[c] = acc
        acc = merge(acc, dn[c])
      dp[v], acc = fin(acc, v), e
      for j in range(lo[i + 1] - 1, lo[i] - 1, -1):
        c = q[j]
        up[c] = lift(fin(merge(up[c], acc), v), v, c)
        acc = merge(dn[c], acc)
    self.par, self.dp, self.dn, self.up = par, dp, dn, up
