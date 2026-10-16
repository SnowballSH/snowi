"""
Title: Heavy-light decomposition
Description: `pos[v]`: DFS order, heavy child (largest subtree)
  first; chains (top `head[v]`) and subtrees
  `[pos[v], pos[v] + size[v])` are contiguous. `path(u, v)`:
  $O(log n)$ half-open ranges listed from `u` to `v`; those from
  `u` up to the LCA run against position order (walk each right
  to left). `edges=True` leaves out the LCA (edge values at their
  lower endpoint). `climb(x, w, keep)`: the upward ranges from `x`
  to its ancestor `w` (`w` itself if `keep`). Also `par` (-1 at
  the root), `depth`, `heavy` (-1 at leaves).
Precondition: `adj` holds both directions of the $n - 1$ edges of
  a tree on vertices `0..n-1`, $n >= 1$.
Time: $O(n)$ build (iterative); `path` returns $O(log n)$ ranges
  and `lca` takes $O(log n)$. PyPy, $n = 2 dot 10^5$: build about
  0.3 s; $2 dot 10^5$ path sums with ds/segtree about 1 s.
Usage:
  h = HLD(adj)  # rooted at 0; HLD(adj, root) otherwise
  b = [0] * n  # ds/segtree (sum) over positions
  for v in range(n): b[h.pos[v]] = a[v]
  st = SegTree(b, add, 0)  # point update: st.set(h.pos[v], x)
  s = sum(st.query(l, r) for l, r in h.path(u, v))
  l, r = h.subtree(v)  # st.query(l, r)
  # ds/lazy-segtree lz: add x to every vertex on a path
  for l, r in h.path(u, v): lz.apply(l, r, x)
Tweak: Edge values: put the value of edge `(par[v], v)` at
  `pos[v]` and pass `edges=True`; the subtree's edges are then
  `[pos[v] + 1, pos[v] + size[v])`. Non-commutative `op`: with
  `w = lca(u, v)`, fold `path(u, w)` in a second ds/segtree with
  `lambda f, g: op(g, f)`, then `path(w, v, True)` in the normal
  one.
Example:
  >>> adj = [[1, 2], [0, 3, 4], [0], [1], [1, 5], [4]]
  >>> h = HLD(adj)  # edges 0-1 0-2 1-3 1-4 4-5
  >>> h.pos, h.head  # chain 0 1 4 5, then 3, then 2
  ([0, 1, 5, 4, 2, 3], [0, 0, 2, 3, 0, 0])
  >>> h.path(3, 2), h.path(3, 2, True), h.subtree(1)
  ([(4, 5), (0, 2), (5, 6)], [(4, 5), (1, 2), (5, 6)], (1, 5))
  >>> a = [10, 20, 30, 40, 50, 60]  # vertex values
  >>> b = [0] * 6
  >>> for v in range(6): b[h.pos[v]] = a[v]
  >>> sum(sum(b[l:r]) for l, r in h.path(5, 3))  # 60+50+20+40
  170
  >>> w = h.lca(3, 2)  # walk 3 -> 2: upward ranges reversed
  >>> seq = [b[l:r][::-1] for l, r in h.path(3, w)]
  >>> seq += [b[l:r] for l, r in h.path(w, 2, True)]
  >>> w, seq
  (0, [[40], [20, 10], [30]])
Notes: Order matters only for non-commutative operations.
  Positions grow downward along a chain, so `path(u, v)` is
  `climb(u, w)` (each range walked right to left) followed by
  `climb(v, w)` reversed (walked left to right), `w = lca(u, v)`.
  Hence the Tweak's split: `path(u, w)` is all upward,
  `path(w, v, True)` all downward without `w`. Cues: "change the
  values on the path from `u` to `v`", "sum or max on a path with
  updates", "add to a subtree, query a path". Only subtree
  queries: plain DFS order suffices. Static values with path
  sums: prefix sums from the root and an LCA are simpler. A root
  path crosses at most $log_2 n$ light edges, so any path crosses
  at most $2 log_2 n$, and with a segment tree a path query costs
  $O(log^2 n)$. Subtree of `v` when the root moves to `r`: if `r`
  is a proper descendant of `v`, it is everything outside the
  subtree of `v`'s child `c` towards `r` (climb from `r`: `c` is
  the last head below `v`, or `heavy[v]` once on `v`'s chain).
Source: kactl HLD (CC0); maspy's library (CC0).
Status: stress-tested against explicit path enumeration
"""


class HLD:
  def __init__(self, adj, root=0):
    n = len(adj)
    par, depth, order = [-1] * n, [0] * n, [root]
    for u in order:  # BFS: also visits vertices appended below
      for v in adj[u]:
        if v != par[u]:
          par[v] = u
          depth[v] = depth[u] + 1
          order.append(v)
    size, heavy = [1] * n, [-1] * n
    for i in range(n - 1, 0, -1):
      v = order[i]
      p = par[v]
      size[p] += size[v]
      if heavy[p] < 0 or size[v] > size[heavy[p]]:
        heavy[p] = v
    head, pos, t, st = [root] * n, [0] * n, 0, [root]
    while st:
      v = st.pop()
      pos[v] = t
      t += 1
      for c in adj[v]:
        if c != par[v] and c != heavy[v]:
          head[c] = c
          st.append(c)
      if heavy[v] >= 0:
        head[heavy[v]] = head[v]
        st.append(heavy[v])
    self.par, self.heavy, self.head = par, heavy, head
    self.pos, self.depth, self.size = pos, depth, size

  def path(self, u, v, edges=False):
    w = self.lca(u, v)
    a = self.climb(u, w, not edges)
    a.extend(reversed(self.climb(v, w, False)))
    return a

  def climb(self, x, w, keep):  # x up to its ancestor w
    head, pos, par, out = self.head, self.pos, self.par, []
    while head[x] != head[w]:
      out.append((pos[head[x]], pos[x] + 1))
      x = par[head[x]]
    if x != w or keep:
      out.append((pos[w] + (not keep), pos[x] + 1))
    return out

  def subtree(self, v):
    return self.pos[v], self.pos[v] + self.size[v]

  def lca(self, u, v):
    head, par, depth = self.head, self.par, self.depth
    while head[u] != head[v]:
      if depth[head[u]] > depth[head[v]]:
        u, v = v, u
      v = par[head[v]]
    return u if depth[u] < depth[v] else v
