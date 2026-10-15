"""
Title: Hopcroft-Karp matching
Description: Maximum matching in a bipartite graph with `n` left
  and `m` right vertices. `ml[u]` is the right partner of left `u`
  and `mr[v]` the left partner of right `v` (-1 if free). `cover()`
  then gives a minimum vertex cover (Kőnig): left and right
  vertices touching every edge, as many as the matching.
Precondition: Edges go from left `u` in `[0, n)` to right `v` in
  `[0, m)`; both sides are numbered from 0 separately.
Time: $O(m' sqrt(n + m))$ for $m'$ edges.
Usage:
  hk = HopcroftKarp(n, m)
  hk.add_edge(u, v)
  k = hk.max_matching()  # then hk.ml, hk.mr
  cl, cr = hk.cover()  # after max_matching
Tweak: Maximum independent set: the vertices not in `cover()`,
  $n + m - k$ of them. Several units per left vertex (b-matching)
  or weights: use Dinic or min-cost flow.
Example:
  >>> hk = HopcroftKarp(3, 3)
  >>> for u, v in [(0, 0), (0, 1), (1, 0), (2, 0)]:
  ...   hk.add_edge(u, v)
  >>> hk.max_matching(), hk.ml, hk.mr
  (2, [1, 0, -1], [1, 0, -1])
  >>> hk.cover()  # left 0 and right 0 touch all four edges
  ([0], [0])
  >>> # Minimum path cover of the DAG 0->1, 1->2, 0->3:
  >>> # n minus a matching of out-copies to in-copies.
  >>> pc = HopcroftKarp(4, 4)
  >>> for u, v in [(0, 1), (1, 2), (0, 3)]: pc.add_edge(u, v)
  >>> 4 - pc.max_matching()
  2
Notes: Matching facts (bipartite only): minimum vertex cover =
  maximum matching (Kőnig); maximum independent set and minimum
  edge cover (no isolated vertices) are both $n + m - k$. Hall: all
  left vertices can be matched iff every left set $X$ has
  $|N(X)| >= |X|$, and in general $k = n - max_X (|X| - |N(X)|)$.
  Cues: "pair up", "assign each to one", "maximum set with no two
  conflicting" when conflicts only join two classes (grid cells by
  color parity, even and odd numbers), "fewest rows and columns
  covering all marks" (Kőnig). In a path cover, paths start at the
  `v` with `mr[v] < 0` and continue `v -> ml[v]`. The cover is read
  off the last search: left vertices not reached by alternating
  paths from free left vertices, plus right vertices reached.
Source: kactl hopcroftKarp (CC0), restructured.
Status: stress-tested against bitmask DP and Kőnig certificates
"""


class HopcroftKarp:
  def __init__(self, n, m):
    self.n, self.m = n, m
    self.g = [[] for _ in range(n)]
    self.ml, self.mr = [-1] * n, [-1] * m

  def add_edge(self, u, v):
    self.g[u].append(v)

  def max_matching(self):
    g, n, ml, mr = self.g, self.n, self.ml, self.mr
    while True:
      lvl = self.lvl = [-1] * n
      q = [u for u in range(n) if ml[u] < 0]
      for u in q:
        lvl[u] = 0
      lim = n  # layer of the shortest augmenting paths
      for u in q:
        lu = lvl[u]
        if lu >= lim:
          break
        for v in g[u]:
          w = mr[v]
          if w < 0:
            lim = lu
          elif lvl[w] < 0:
            lvl[w] = lu + 1
            q.append(w)
      if lim == n:
        return n - ml.count(-1)
      it = [0] * n
      for root in range(n):
        if ml[root] >= 0:
          continue
        stack = [root]  # stack[i] is on layer i
        while stack:
          u, d = stack[-1], len(stack)
          a, i = g[u], it[u]
          while i < len(a):
            w = mr[a[i]]
            if w < 0 and d - 1 == lim or w >= 0 and lvl[w] == d:
              break
            i += 1
          it[u] = i
          if i == len(a):
            stack.pop()
            if stack:
              it[stack[-1]] += 1
          elif mr[a[i]] < 0:
            for x in stack:
              v = g[x][it[x]]
              ml[x], mr[v] = v, x
            break
          else:
            lvl[mr[a[i]]] = -1  # visit once per phase
            stack.append(mr[a[i]])

  def cover(self):
    lvl, mr = self.lvl, self.mr
    cl = [u for u in range(self.n) if lvl[u] < 0]
    cr = [v for v, u in enumerate(mr) if u >= 0 and lvl[u] >= 0]
    return cl, cr
