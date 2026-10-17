"""
Title: Dominator tree
Description: Lengauer-Tarjan. Returns `idom`: the immediate
  dominator of each vertex (the closest strict dominator), with
  `idom[s] = s` and `idom[v] = -1` if `v` is unreachable from `s`.
  The edges `idom[v] -> v` form a tree rooted at `s`; the dominators
  of `v` are its ancestors.
Precondition: `adj[u]` lists the heads of the edges leaving `u`;
  vertices 0-indexed. Self-loops and parallel edges are fine.
Time: $O((n + m) log n)$, no recursion; PyPy needs about 1 s for
  $n = 2 dot 10^5$, $m = 6 dot 10^5$.
Usage:
  idom = dominator_tree(adj, s)
Tweak: Paths into a sink `t` instead (post-dominators): run it on
  the reversed graph from `t`. A DAG needs no code: in topological
  order, `idom[v]` is the LCA of the predecessors of `v` that are
  reachable from `s` (binary lifting).
Example:
  >>> adj = [[1, 2], [3], [3], [4, 1], [], [4]]
  >>> idom = dominator_tree(adj, 0); idom
  [0, 0, 0, 0, 3, -1]
  >>> def critical(t):  # vertices on every path 0 -> t
  ...   crit = []
  ...   while idom[t] >= 0:  # stays [] if t is unreachable
  ...     crit.append(t)
  ...     if t == idom[t]: break
  ...     t = idom[t]
  ...   return crit
  >>> critical(4), critical(5)
  ([4, 3, 0], [])
  >>> # vertices cut off from 0 when v fails: subtree sizes
  >>> ch = [[] for _ in idom]
  >>> for v, d in enumerate(idom):
  ...   if d >= 0 and d != v: ch[d].append(v)
  >>> order = [0]  # BFS of the dominator tree
  >>> for u in order: order += ch[u]
  >>> lost = [int(d >= 0) for d in idom]
  >>> for v in reversed(order[1:]): lost[idom[v]] += lost[v]
  >>> lost
  [5, 1, 1, 2, 1, 0]
Notes: Think of it for "which vertices (or edges) must every route
  from `s` use", "removing which node disconnects `v` from the
  capital", or "for each vertex, how many become unreachable if it
  fails": the subtree sizes of the dominator tree. On an undirected
  graph use articulation points (`graph/bcc`) instead. Edges:
  subdivide edge `u -> v` with a new vertex `x` (`u -> x -> v`); the
  edge is critical for `v` iff `x` dominates `v`.
Source: Lengauer and Tarjan (1979), the simple $O(m log n)$ version
  with path compression.
Status: stress-tested against removal-based brute force
"""


def dominator_tree(adj, s):
  n = len(adj)
  tin, it = [-1] * n, [0] * n
  tin[s] = 0  # DFS numbers; below, everything is in them
  order, par, st = [s], [-1], [s]
  while st:
    u = st[-1]
    out, i = adj[u], it[u]
    if i == len(out):
      st.pop()
      continue
    it[u] = i + 1
    v = out[i]
    if tin[v] < 0:
      tin[v] = len(order)
      order.append(v)
      par.append(tin[u])
      st.append(v)
  k = len(order)
  start = [0] * (k + 1)  # predecessors of DFS number w, flat
  for u in order:
    for v in adj[u]:
      if tin[v] >= 0:
        start[tin[v] + 1] += 1
  for i in range(k):
    start[i + 1] += start[i]
  pred, fill = [0] * start[k], start[:]
  for u in order:
    for v in adj[u]:
      t = tin[v]
      if t >= 0:
        pred[fill[t]] = tin[u]
        fill[t] += 1
  semi = list(range(k))
  anc, best = semi[:], semi[:]
  idom, head, nxt, path = [0] * k, [-1] * k, [-1] * k, []

  def eval_(v):  # min semi on v's forest path
    x = v
    while anc[x] != x:
      path.append(x)
      x = anc[x]
    while path:
      y = path.pop()
      a = anc[y]
      if semi[best[a]] < semi[best[y]]:
        best[y] = best[a]
      anc[y] = x
    return best[v]

  for w in range(k - 1, 0, -1):
    sw = semi[w]
    for j in range(start[w], start[w + 1]):
      v = pred[j]
      if anc[v] != v:
        v = eval_(v)
      if semi[v] < sw:
        sw = semi[v]
    semi[w] = sw
    nxt[w], head[sw] = head[sw], w  # bucket sw, linked
    p = anc[w] = par[w]
    v = head[p]
    while v >= 0:
      u = eval_(v)
      idom[v] = u if semi[u] < semi[v] else p
      v = nxt[v]
    head[p] = -1
  res = [-1] * n
  res[s] = s
  for i in range(1, k):
    if idom[i] != semi[i]:
      idom[i] = idom[idom[i]]
    res[order[i]] = order[idom[i]]
  return res
