"""
Title: General matching
Description: Maximum matching in any undirected graph (Edmonds'
  blossom algorithm). Returns the size `k` and `mate`, where
  `mate[v]` is the partner of `v`, or -1 if `v` is unmatched.
Precondition: 0-indexed adjacency lists, each edge in both `g[u]`
  and `g[v]`. Self-loops and parallel edges are fine.
Time: $O(n^3)$, at most one search per vertex; far faster in
  practice: about a second for a random graph with $n = 10^5$,
  $m = 2 dot 10^5$, about 0.2 s for paths and random trees.
Usage:
  g = [[] for _ in range(n)]
  g[u].append(v); g[v].append(u)
  k, mate = general_matching(g)
Tweak: Bipartite graphs: `flow/hopcroft-karp` is faster. Weighted
  matching in a general graph needs the much longer weighted
  blossom algorithm; for $n <= 20$ use a bitmask DP.
Example:
  >>> tri = [[1, 2], [0, 2], [0, 1, 3], [2, 4], [3, 5], [4]]
  >>> general_matching(tri)  # triangle 0 1 2, path 2 3 4 5
  (3, [1, 0, 3, 2, 5, 4])
  >>> # Players alternately move a token to an unvisited
  >>> # neighbour; whoever cannot move loses. The first player
  >>> # wins from s iff every maximum matching covers s.
  >>> g = [[1], [0, 2], [1]]  # path 0 - 1 - 2
  >>> k = general_matching(g)[0]
  >>> for s in range(3):
  ...   h = [[v for v in a if v != s] for a in g]
  ...   h[s] = []
  ...   print(s, general_matching(h)[0] < k)
  0 False
  1 True
  2 False
Notes: Do not make a general graph bipartite by doubling the
  vertices: a triangle then gets 3 edges instead of 1 pair.
  Tutte-Berge: $2k = min_U (n + |U| - "odd"(G - U))$, where
  $"odd"(G - U)$ counts the components of odd size after deleting
  $U$; a perfect matching exists iff $"odd"(G - U) <= |U|$ for all
  $U$. The smallest edge cover (no isolated vertices) has $n - k$
  edges, but Kőnig fails: the largest independent set is NP-hard
  here. Cues: "pair up" when any two items may pair (not two
  classes), dominoes on a graph with odd cycles, the token game
  above. Each vertex starts at most one search: a vertex with no
  augmenting path now never gets one later, and the whole tree of
  a failed search can be dropped for good (Edmonds). Hand-written
  versions that rebuild it are quadratic on a long path with many
  leaves at one end.
Source: Edmonds (1965); union-find blossom bases as in the common
  BFS formulation (Gabow 1976).
Status: stress-tested against bitmask DP and Tutte matrix rank
"""


def general_matching(g):
  n = len(g)
  mate, pre, dsu = [-1] * n, [0] * n, list(range(n))
  lab = [0] * n  # 0 unseen, 1 even, 2 odd, 3 dead
  seen, q, stamp = [0] * n, [], 0

  def base(x):
    while dsu[x] != x:
      dsu[x] = dsu[dsu[x]]
      x = dsu[x]
    return x

  def lca(x, y):
    nonlocal stamp
    stamp += 1
    while True:
      if x >= 0:
        x = base(x)
        if seen[x] == stamp:
          return x
        seen[x] = stamp
        x = pre[mate[x]] if mate[x] >= 0 else -1
      x, y = y, x

  def shrink(x, y, b):
    while base(x) != b:
      pre[x] = y
      y = mate[x]
      if lab[y] == 2:
        lab[y] = 1
        q.append(y)
      if dsu[x] == x:
        dsu[x] = b
      if dsu[y] == y:
        dsu[y] = b
      x = pre[y]

  def augment(root):
    for x in q:  # undo the last search
      lab[x] = 0
      dsu[x] = x
    q[:] = [root]
    lab[root] = 1
    for u in q:
      if lab[u] != 1:
        continue
      for v in g[u]:
        if lab[v] == 0:
          lab[v], pre[v] = 2, u
          q.append(v)
          w = mate[v]
          if w >= 0:
            lab[w] = 1
            q.append(w)
            continue
          while v >= 0:
            u = pre[v]
            w = mate[u]
            mate[v], mate[u], v = u, v, w
          return 1
        if lab[v] == 1 and base(u) != base(v):
          b = lca(u, v)
          shrink(u, v, b)
          shrink(v, u, b)
    for x in q:  # no augmenting path through them, ever
      lab[x] = 3
    q[:] = []
    return 0

  k = 0
  for v in range(n):
    if mate[v] < 0:
      k += augment(v)
  return k, mate
