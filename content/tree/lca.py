"""
Title: LCA by binary lifting
Description: Lowest common ancestor, distance in edges, k-th
  ancestor, and `jump(u, v, k)`: the vertex `k` edges from `u` on
  the path to `v`. Also gives `depth`, the parent `up[0][v]` (the
  root's is itself) and a DFS preorder: `order[tin[v]] = v`, and
  every subtree is a contiguous block of `tin` values.
Precondition: `adj` holds both directions of the $n - 1$ edges of
  a tree on vertices `0..n-1`, $n >= 1$. `k >= 0`.
Time: $O(n log n)$ time and memory to build (iterative), then
  $O(log n)$ per query. PyPy, $n = 2 dot 10^5$: about 0.2 s to
  build and 0.3 s for $2 dot 10^5$ `lca` queries.
Usage:
  lc = LCA(adj)  # rooted at 0; LCA(adj, root) otherwise
  lc.lca(u, v); lc.dist(u, v); lc.depth[v]
  lc.kth_ancestor(u, k)  # -1 if k > depth[u]
  lc.jump(u, v, k)  # -1 if k > dist(u, v); k = 0 gives u
Tweak: Weighted distance: also keep `wd[v]`, the weighted depth
  (set `wd[v] = wd[u] + w` next to `depth`), and return `wd[u] +
  wd[v] - 2 * wd[lca(u, v)]`. Path aggregates (max edge on a
  path): lift a second table `mx[k][v]` beside `up` and fold it
  in during the same jumps.
Example:
  >>> adj = [[1, 2], [0, 3, 4], [0], [1], [1, 5], [4]]
  >>> lc = LCA(adj)  # edges 0-1 0-2 1-3 1-4 4-5
  >>> lc.lca(3, 5), lc.lca(5, 2), lc.dist(3, 5), lc.dist(2, 5)
  (1, 0, 3, 4)
  >>> lc.kth_ancestor(5, 2), lc.kth_ancestor(5, 4)
  (1, -1)
  >>> [lc.jump(2, 5, k) for k in range(6)]  # path 2 0 1 4 5
  [2, 0, 1, 4, 5, -1]
  >>> # max edge weight on a path; w[v]: edge v - parent
  >>> w = [0, 7, 1, -2, -5, -3]
  >>> mx = [w]  # mx[k][v]: max of the 2^k edges above v
  >>> for k in range(1, len(lc.up)):
  ...   m, up = mx[-1], lc.up[k - 1]
  ...   mx.append([max(m[v], m[up[v]]) for v in range(6)])
  >>> def climb(v, d):  # max of the d edges above v
  ...   r, k = -(1 << 62), 0  # -INF: weights may be negative
  ...   while d:
  ...     if d & 1: r, v = max(r, mx[k][v]), lc.up[k][v]
  ...     d, k = d >> 1, k + 1
  ...   return r
  >>> def path_max(u, v):
  ...   a = lc.depth[lc.lca(u, v)]
  ...   du, dv = lc.depth[u] - a, lc.depth[v] - a
  ...   return max(climb(u, du), climb(v, dv))
  >>> path_max(3, 5), path_max(2, 5)  # 3-1-4-5, 2-0-1-4-5
  (-2, 7)
Notes: Cues: "distance between two vertices", "k-th vertex on
  the path", "is `x` on the path `a - b`" (iff `dist(a, x) +
  dist(x, b) == dist(a, b)`), "middle of a path" (`jump(u, v,
  dist(u, v) // 2)`; on a diameter `u - v` it is a center). `x`
  is an ancestor of `v` iff `depth[x] <= depth[v]` and
  `kth_ancestor(v, depth[v] - depth[x]) == x`, or in $O(1)$ with
  subtree sizes: `tin[x] <= tin[v] < tin[x] + size[x]` (sizes:
  add each `order[i]` into its parent for `i` from `n - 1` down
  to 1). Paths `a - b` and `c - d` meet iff `lca(a, b)` lies on
  `c - d` or `lca(c, d)` lies on `a - b`. The same lifting finds
  the highest ancestor with a property that is monotone along
  the root path (value below `x`, Kruskal reconstruction trees):
  for `k` from the top, jump while `up[k][v]` still satisfies
  it. Memory is $n log_2 n$ ints; for $10^6$ queries use
  tree/lca-rmq ($O(1)$ per query) or `HLD.lca` ($O(n)$ memory).
  For a forest, join one vertex of each tree to a new root `n`
  and test reachability separately.
Source: Folklore binary lifting; cp-algorithms.
Status: stress-tested against naive parent walking
"""


class LCA:
  def __init__(self, adj, root=0):
    n = len(adj)
    par, depth, tin, order = [root] * n, [0] * n, [0] * n, []
    st = [root]
    while st:
      u = st.pop()
      tin[u] = len(order)
      order.append(u)
      for v in adj[u]:
        if v != par[u]:
          par[v] = u
          depth[v] = depth[u] + 1
          st.append(v)
    up = [par]  # up[k][v]: 2^k-th ancestor, capped at the root
    for _ in range(1, max(n - 1, 1).bit_length()):
      prev = up[-1]
      up.append([prev[p] for p in prev])
    self.depth, self.tin = depth, tin
    self.order, self.up = order, up

  def kth_ancestor(self, u, k):
    if k > self.depth[u]:
      return -1
    up, j = self.up, 0
    while k:
      if k & 1:
        u = up[j][u]
      k >>= 1
      j += 1
    return u

  def lca(self, u, v):
    depth = self.depth
    if depth[u] < depth[v]:
      u, v = v, u
    u = self.kth_ancestor(u, depth[u] - depth[v])
    if u == v:
      return u
    for row in self.up[depth[v].bit_length() - 1 :: -1]:
      if row[u] != row[v]:
        u, v = row[u], row[v]
    return self.up[0][u]

  def dist(self, u, v):
    d = self.depth
    return d[u] + d[v] - 2 * d[self.lca(u, v)]

  def jump(self, u, v, k):
    depth = self.depth
    a = self.lca(u, v)
    du = depth[u] - depth[a]
    d = du + depth[v] - depth[a]
    if k > d:
      return -1
    if k <= du:
      return self.kth_ancestor(u, k)
    return self.kth_ancestor(v, d - k)
