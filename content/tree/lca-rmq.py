"""
Title: LCA by Euler tour and RMQ
Description: Lowest common ancestor in $O(1)$ after an
  $O(n log n)$ build. With a DFS preorder `tin`, entry `i` of the
  tour holds `tin` of the parent of the vertex entered at time
  `i + 1`; for `u != v`, the LCA is the vertex whose `tin` is the
  minimum of entries `[min(tin), max(tin))`. Like tree/lca it
  has `depth`, `tin` and `order` (`order[tin[v]] = v`); it also
  has `par`, -1 at the root (tree/lca has `up[0]` instead,
  which maps the root to itself).
Precondition: `adj` holds both directions of the $n - 1$ edges of
  a tree on vertices `0..n-1`, $n >= 1$.
Time: $O(n log n)$ time and memory to build (iterative), $O(1)$
  per query. PyPy, $n = 2 dot 10^5$: about 0.25 s to build and
  0.1 s for $2 dot 10^5$ `lca` queries.
Usage:
  lc = RmqLCA(adj)  # rooted at 0; RmqLCA(adj, root) otherwise
  lc.lca(u, v); lc.dist(u, v); lc.depth[v]
Tweak: Weighted distance: keep the weighted depth `wd` beside
  `depth` and return `wd[u] + wd[v] - 2 * wd[lca(u, v)]`.
Uses: ds/sparse-table
Example:
  >>> adj = [[1, 2], [0, 3, 4], [0], [1], [1, 5], [4]]
  >>> lc = RmqLCA(adj)  # edges 0-1 0-2 1-3 1-4 4-5
  >>> lc.lca(3, 5), lc.lca(5, 2), lc.lca(4, 4), lc.dist(2, 5)
  (1, 0, 4, 4)
  >>> # Steiner tree: edges needed to connect the set s
  >>> # = half the closed tour through s in preorder
  >>> s = sorted([3, 5, 2], key=lambda v: lc.tin[v])
  >>> sum(lc.dist(s[i - 1], s[i]) for i in range(len(s))) // 2
  5
Notes: Prefer it over binary lifting when queries dominate
  ($10^6$ LCA or distance queries) or inside another loop (virtual
  trees, sorting by distance). The tour has $n - 1$ entries, half
  the classic $2n - 1$ tour of depths. Memory is $n log_2 n$ ints:
  about 150 MB at $n = 10^6$ in PyPy. It answers only LCA and
  distance; for k-th ancestors and `jump` use tree/lca.
Source: kactl LCA (CC0), made iterative.
Status: stress-tested against naive parent walking
"""


class RmqLCA:
  def __init__(self, adj, root=0):
    n = len(adj)
    par, depth, tin, order = [-1] * n, [0] * n, [0] * n, []
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
    self.par, self.depth = par, depth
    self.tin, self.order = tin, order
    self.rmq = SparseTable([tin[par[v]] for v in order[1:]])

  def lca(self, u, v):
    if u == v:
      return u
    a, b = self.tin[u], self.tin[v]
    if a > b:
      a, b = b, a
    return self.order[self.rmq.query(a, b)]

  def dist(self, u, v):
    d = self.depth
    return d[u] + d[v] - 2 * d[self.lca(u, v)]
