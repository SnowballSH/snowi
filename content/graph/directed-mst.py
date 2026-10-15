"""
Title: Directed MST
Description: Minimum arborescence rooted at `root`
  (Chu-Liu/Edmonds): every vertex except `root` takes one incoming
  edge `par[v] -> v` and all vertices are reached from `root`, at
  the smallest total weight. Returns the cost and `par`
  (`par[root] = -1`), or `(-1, [])` if some vertex cannot be reached
  from `root`.
Precondition: Vertices are 0-indexed, edges `a -> b` with integer
  weights. With negative weights a cost of -1 is possible, so test
  `not par`. Self-loops and parallel edges are fine.
Time: $O(m log m)$; PyPy needs about 1 s for $n = 2 dot 10^5$,
  $m = 5 dot 10^5$.
Usage:
  g = DirectedMST(n)
  g.add_edge(a, b, w)  # a -> b
  cost, par = g.solve(root)
Tweak: Maximum arborescence: negate the weights. Free root: add a
  vertex `n` with edges `n -> v` of weight $B = 1 + sum |w|$ and
  solve `DirectedMST(n + 1)` from `n`. One root suffices iff exactly
  one `v` has `par[v] == n`; that `v` is the best root and the cost
  is `cost - B`.
Uses: ds/dsu-rollback
Example:
  >>> g = DirectedMST(4)
  >>> for a, b, w in [(0, 1, 5), (1, 2, 1), (2, 3, 1), (3, 1, 1),
  ...                 (0, 3, 4)]:
  ...   g.add_edge(a, b, w)
  >>> g.solve(0)  # the cycle 1 -> 2 -> 3 -> 1 is entered at 3
  (6, [-1, 3, 1, 0])
  >>> g.solve(1)  # nothing enters 0
  (-1, [])
  >>> # a well at v for c[v], or a one-way pipe u -> v for w;
  >>> # vertex 3 is a super root
  >>> c, pipes = [3, 9, 9], [(0, 1, 2), (1, 2, 2), (2, 1, 1)]
  >>> g = DirectedMST(4)
  >>> for v in range(3): g.add_edge(3, v, c[v])
  >>> for a, b, w in pipes: g.add_edge(a, b, w)
  >>> g.solve(3)
  (7, [3, 0, 1, -1])
Notes: Think of it when every vertex needs exactly one supplier
  along directed edges and all must be reachable from a source:
  "broadcast from 1", "each node copies from one other node", "a
  one-way pipe from a neighbour". An undirected MST is wrong for
  one-way edges. "Build your own well for `c[v]` or take a pipe": a
  super root with edges of weight `c[v]` (several roots allowed, no
  `B` needed). Idea: each vertex picks its cheapest incoming edge; a
  cycle among the picks is contracted, and an edge entering it is
  discounted by the cycle edge it replaces (the lazy add on the
  heap). The number of arborescences comes from the directed
  matrix-tree theorem instead.
Source: kactl DirectedMST (CC0), with index-based skew heaps and
  an iterative merge.
Status: stress-tested against brute force over all parent choices
"""


class DirectedMST:
  def __init__(self, n):
    self.n = n
    self.src, self.dst, self.w = [], [], []

  def add_edge(self, a, b, w):
    self.src.append(a)
    self.dst.append(b)
    self.w.append(w)

  def solve(self, root):
    n, src, dst = self.n, self.src, self.dst
    m = len(src)
    key, lazy = self.w[:], [0] * m
    left, right = [-1] * m, [-1] * m

    def prop(x):
      d = lazy[x]
      if d:
        key[x] += d
        if left[x] >= 0:
          lazy[left[x]] += d
        if right[x] >= 0:
          lazy[right[x]] += d
        lazy[x] = 0

    def merge(a, b):  # skew heaps; nodes are edge ids
      st = []
      while a >= 0 and b >= 0:
        prop(a)
        prop(b)
        if key[a] > key[b]:
          a, b = b, a
        st.append(a)
        a = right[a]
      if a < 0:
        a = b
      while st:
        x = st.pop()
        right[x] = left[x]
        left[x] = a
        a = x
      return a

    uf = RollbackDSU(n)
    find = uf.find
    heap, seen, inn = [-1] * n, [-1] * n, [-1] * n
    path, q, cycles = [0] * n, [0] * n, []
    for e in range(m):
      heap[dst[e]] = merge(heap[dst[e]], e)
    cost, seen[root] = 0, root
    for s in range(n):
      u, k = s, 0
      while seen[u] < 0:
        e = heap[u]
        if e < 0:
          return -1, []
        prop(e)
        cost += key[e]
        lazy[e] = -key[e]
        prop(e)
        heap[u] = merge(left[e], right[e])
        q[k], path[k], seen[u] = e, u, s
        k += 1
        u = find(src[e])
        if seen[u] == s:  # contract the cycle into one vertex
          cyc, end, t = -1, k, uf.time()
          while True:
            k -= 1
            x = path[k]
            cyc = merge(cyc, heap[x])
            if not uf.join(u, x):
              break
          u = find(u)
          heap[u], seen[u] = cyc, -1
          cycles.append((u, t, q[k:end]))
      for i in range(k):
        inn[find(dst[q[i]])] = q[i]
    for u, t, cyc in reversed(cycles):
      uf.rollback(t)
      e = inn[u]
      for c in cyc:
        inn[find(dst[c])] = c
      inn[find(dst[e])] = e
    return cost, [src[e] if e >= 0 else -1 for e in inn]
