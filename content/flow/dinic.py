"""
Title: Dinic max flow
Description: Maximum flow from `s` to `t`. Afterwards `flow(id)`
  is the flow on edge `id` (negative if it runs backwards on an
  undirected edge), and `source_side(v)` tells whether `v` is on
  the source side of a minimum cut: reachable from `s` in the
  residual graph.
Precondition: Vertices 0-indexed, `s != t`, capacities `>= 0`.
Time: $O(n^2 m)$, far faster in practice. Unit capacities:
  $O(m min(sqrt(m), n^(2/3)))$, and $O(m sqrt(n))$ if every
  vertex but $s$, $t$ has a single edge in or a single edge out
  (matching, split vertices).
Usage:
  d = Dinic(n)
  i = d.add_edge(u, v, cap)  # directed u -> v
  d.add_edge(a, b, c, c)  # undirected, capacity c
  f = d.max_flow(s, t)
  d.flow(i); d.source_side(v)
Tweak: A second `max_flow` call keeps the flow already sent and
  adds to it. To start over on the same graph, copy `orig` into
  `cap`.
Example:
  >>> d = Dinic(4)
  >>> a = d.add_edge(0, 1, 3); b = d.add_edge(0, 2, 2)
  >>> c = d.add_edge(1, 2, 5); _ = d.add_edge(1, 3, 2)
  >>> _ = d.add_edge(2, 3, 3)
  >>> d.max_flow(0, 3), d.flow(a), d.flow(b), d.flow(c)
  (5, 3, 2, 1)
  >>> [v for v in range(4) if d.source_side(v)]
  [0]
  >>> # Undirected unit edges: two edge-disjoint paths 0 -> 2
  >>> u = Dinic(3)
  >>> for x, y in [(0, 1), (1, 2), (0, 2)]:
  ...   _ = u.add_edge(x, y, 1, 1)
  >>> u.max_flow(0, 2)
  2
Notes: Read a minimum cut as the edges from source-side to
  sink-side vertices, not as "saturated edges": flow paths can
  saturate edges anywhere. An undirected edge is one call with
  `rev_cap = c`, not two calls. Keep an "infinite" capacity such
  as `1 << 60` below $2^63$ in sums: bigger ints are correct but
  slow. The search is iterative and, like the C++ recursion,
  splits flow over many paths at once, so long paths and long
  shared prefixes (time chains) are fine. In PyPy a random
  matching network with $10^5$ edges takes about a second; on
  big grids with varied capacities Dinic can take hundreds of
  phases. Think of max flow for "can every demand be met", "pair
  things up, each usable $k$ times", "fewest edges or vertices
  to delete to separate", "pick a subset where taking $x$ forces
  $y$": see the modelling note. With $n$ up to $10^6$ the
  intended solution is usually a greedy or a formula proved by
  the cut condition, not a flow.
Source: kactl Dinic (CC0); edge ids as in the AtCoder Library.
Status: stress-tested against cut enumeration and Edmonds-Karp
"""


class Dinic:
  def __init__(self, n):
    self.n, self.head, self.lvl = n, [-1] * n, [-1] * n
    self.to, self.cap, self.orig, self.nxt = [], [], [], []

  def add_edge(self, u, v, cap, rev_cap=0):
    i, head, nxt = len(self.to), self.head, self.nxt
    self.to += [v, u]
    self.cap += [cap, rev_cap]
    self.orig += [cap, rev_cap]
    nxt.append(head[u])
    head[u] = i
    nxt.append(head[v])
    head[v] = i + 1
    return i

  def flow(self, i):
    return self.orig[i] - self.cap[i]

  def source_side(self, v):
    return self.lvl[v] >= 0

  def max_flow(self, s, t):
    n, head, nxt = self.n, self.head, self.nxt
    to, cap, total, lim = self.to, self.cap, 0, [0] * n
    while True:
      lvl = self.lvl = [-1] * n
      lvl[s], q = 0, [s]
      for u in q:
        e, lv = head[u], lvl[u] + 1
        while e >= 0:
          v = to[e]
          if cap[e] and lvl[v] < 0:
            lvl[v] = lv
            q.append(v)
          e = nxt[e]
        if lvl[t] >= 0:
          break
      if lvl[t] < 0:
        return total
      it, path, u, e = head[:], [], s, head[s]
      lim[s] = 0  # lim[v]: what v on the path may still send
      while e >= 0:
        lim[s] += cap[e]
        e = nxt[e]
      sent = lim[s]
      while True:  # path: edges s -> u; r: what u has left
        if u == t:
          r = 0
        else:
          e, lv, r = it[u], lvl[u] + 1, lim[u]
          while e >= 0 and not (cap[e] and lvl[to[e]] == lv):
            e = nxt[e]
          it[u] = e
          if e >= 0:
            c, u = cap[e], to[e]
            lim[u] = c if c < r else r
            path.append(e)
            continue
        while path:  # u is done: apply what it sent over e
          e = path.pop()
          u = to[e ^ 1]
          c, left = cap[e], lim[u]
          d = (c if c < left else left) - r  # u's limit - r
          cap[e] = c - d
          cap[e ^ 1] += d
          r = left - d
          if r:  # e is saturated or leads to a dead end
            lim[u] = r
            it[u] = nxt[e]
            break
        else:
          total += sent - r
          break
