"""
Title: Biconnected components
Description: Blocks of an undirected multigraph given as an edge
  list. `comp[e]` is the block of edge `e` (two edges share a block
  iff some simple cycle contains both); `bridge[e]` and `cut[v]`
  mark bridges and articulation points. `tree` is the block-cut
  forest: node `v < n` is vertex `v`, node `n + b` is block `b`, and
  each block is joined to all its vertices.
Precondition: Vertices 0-indexed. Self-loops are ignored
  (`comp[e] = -1`). Parallel edges form a cycle: they share a block
  and are never bridges. Isolated vertices are in no block.
Time: $O(n + m)$, iterative.
Usage:
  b = BCC(n, edges)  # edges: list of (u, v)
  b.comp[e]; b.bridge[e]; b.cut[v]  # b.nb blocks
  b.tree[n + b.comp[e]]  # the vertices of e's block
Tweak: 2-edge-connected components: join the endpoints of every
  non-bridge edge in a DSU. Vertex `v`'s blocks are its neighbours in
  `tree`; `v` is a cut vertex iff it has two or more.
Example:
  >>> # triangle 0-1-2, bridge 2-3, double edge 3-4, loop at 4
  >>> edges = [(0, 1), (1, 2), (2, 0), (2, 3)]
  >>> edges += [(3, 4), (4, 3), (4, 4)]
  >>> b = BCC(6, edges)
  >>> b.nb, b.comp
  (3, [2, 2, 2, 1, 0, 0, -1])
  >>> [e for e in range(7) if b.bridge[e]], sorted(b.tree[6 + 2])
  ([3], [0, 1, 2])
  >>> [v for v in range(6) if b.cut[v]], b.tree[5]  # 5 isolated
  ([2, 3], [])
Notes: Use the block-cut tree for questions about removing one
  vertex: the vertices other than `u` and `v` whose removal separates
  them are exactly the vertex nodes strictly inside the tree path
  between `u` and `v` (all cut vertices; path queries with LCA). For
  removing one edge use bridges and the
  bridge tree (2-edge-connected components joined by bridges). A
  graph has a strongly connected orientation iff it is connected and
  bridgeless (orient along the DFS). Cactus graphs (each edge on at
  most one cycle): every block is a single edge or a cycle, so DP on
  the block-cut tree handles cycles like vertices. Blocks with 3 or
  more vertices are 2-vertex-connected: any two of their vertices
  lie on a common cycle. Parallel edges must be told apart by edge
  id, which is why the DFS skips only the parent edge, not the
  parent vertex.
Source: Hopcroft-Tarjan, iterative; kactl BiconnectedComponents
  (CC0).
Status: stress-tested against vertex/edge removal and cycle
  enumeration
"""


class BCC:
  def __init__(self, n, edges):
    m = len(edges)
    adj = [[] for _ in range(n)]
    other = [0] * m  # u ^ v: the far end of e is other[e] ^ u
    for e, (u, v) in enumerate(edges):
      if u != v:
        other[e] = u ^ v
        adj[u].append(e)
        adj[v].append(e)
    tin, low, pe = [-1] * n, [0] * n, [-1] * n
    it, bid, head, vs = [0] * n, [-1] * n, [], []
    bridge, t = [False] * m, 0
    for r in range(n):
      if tin[r] >= 0:
        continue
      tin[r] = low[r] = t
      t += 1
      cs = [r]
      while cs:
        u = cs[-1]
        out, i = adj[u], it[u]
        if i < len(out):
          it[u] = i + 1
          e = out[i]
          v = other[e] ^ u
          if tin[v] < 0:
            tin[v] = low[v] = t
            t += 1
            pe[v] = e
            cs.append(v)
            vs.append(v)
          elif e != pe[u] and tin[v] < low[u]:
            low[u] = tin[v]
          continue
        cs.pop()
        if not cs:
          break
        p = cs[-1]
        if low[u] < low[p]:
          low[p] = low[u]
        if low[u] >= tin[p]:  # p cuts off u's subtree: a block
          bridge[pe[u]] = low[u] > tin[p]
          nb = len(head)
          while True:
            v = vs.pop()
            bid[v] = nb
            if v == u:
              break
          head.append(p)
    nb = len(head)
    tree = [[] for _ in range(n + nb)]
    for v in range(n):
      if bid[v] >= 0:
        tree[v].append(n + bid[v])
        tree[n + bid[v]].append(v)
    for b in range(nb):
      tree[head[b]].append(n + b)
      tree[n + b].append(head[b])
    comp = [-1] * m
    for e, (u, v) in enumerate(edges):
      if u != v:  # the deeper endpoint names the block
        comp[e] = bid[u if tin[u] > tin[v] else v]
    self.nb, self.comp, self.bridge = nb, comp, bridge
    self.tree = tree
    self.cut = [len(tree[v]) > 1 for v in range(n)]
