"""
Title: Dijkstra
Description: Shortest paths from a set of sources in a graph with
  non-negative weights. `dist[v]` is `Dijkstra.INF` if `v` cannot
  be reached; `par[v]` is the vertex before `v` on a shortest path
  (-1 for sources and unreachable vertices). `path(t)` lists a
  shortest path from a source to `t`, empty if unreachable.
Precondition: `adj[u]` is a list of pairs `(v, w)` for edges
  `u -> v` (add both directions for an undirected graph). Vertices
  are 0-indexed, weights are integers `w >= 0`, and distances stay
  below `Dijkstra.INF` $= 2^62$.
Time: $O((n + m) log m)$.
Usage:
  adj = [[] for _ in range(n)]
  adj[u].append((v, w))
  dj = Dijkstra(adj, [s])  # or several sources
  d = dj.dist[t]; p = dj.path(t)
Tweak: The heap holds ints `d * n + v`, about 1.5 times faster than
  tuples while the largest distance times `n` is below $2^63$.
  Beyond that PyPy falls back to big ints (still correct, about 2
  times slower). Float weights need tuples `(d, v)` in the heap.
  Bottleneck paths (maximize the smallest edge): `dist` starts at -1
  and sources at `INF`; push `-(d * n + v)` and pop with
  `divmod(-key, n)` (a max-heap), skip if `d < dist[u]`, relax if
  `min(d, w) > dist[v]`, and test `dist[t] >= 0` in `path`.
Example:
  >>> adj = [[(1, 4), (2, 1)], [(3, 1)], [(1, 2), (3, 5)], [], []]
  >>> dj = Dijkstra(adj, [0])
  >>> dj.dist[:4], dj.path(3), dj.path(4)
  ([0, 3, 1, 4], [0, 2, 1, 3], [])
  >>> dj.dist[4] == Dijkstra.INF
  True
  >>> Dijkstra(adj, [1, 2]).dist[1:4]  # several sources
  [0, 0, 1]
Notes: Think of it for "minimum total cost/time" with non-negative
  costs. Several starting points: pass them all as sources
  (multi-source), which also gives "distance to the nearest special
  vertex" for all vertices at once. Distances to one target: run on
  the reversed graph. Extra state (fuel left, parity, coupons used,
  time modulo $k$): run on pairs (vertex, state) numbered
  `v * K + state`. Weights 0 and 1 only: 0-1 BFS with a deque (push
  front on weight 0, back on weight 1) is $O(n + m)$. Edge
  `(u, v, w)` lies on some shortest `s`-`t` path iff
  $d_s (u) + w + d_t (v) = d_s (t)$ (run from `s` and, reversed, from
  `t`). Count shortest paths by DP over the tight edges
  ($d(u) + w = d(v)$) in topological order: increasing `dist` is one
  if all weights are positive, but 0-weight edges also need ties
  ordered along them (a 0-weight cycle gives infinitely many shortest
  walks). A negative edge breaks it silently: use Bellman-Ford or
  potentials (Johnson).
Source: Classic; folklore contest implementation.
Status: stress-tested against Bellman-Ford style brute force
"""


class Dijkstra:
  INF = 1 << 62

  def __init__(self, adj, srcs):
    from heapq import heappop, heappush

    n, INF = len(adj), self.INF
    dist, par, pq = [INF] * n, [-1] * n, []
    for s in srcs:
      dist[s] = 0
      heappush(pq, s)
    while pq:
      d, u = divmod(heappop(pq), n)
      if d > dist[u]:
        continue
      for v, w in adj[u]:
        nd = d + w
        if nd < dist[v]:
          dist[v], par[v] = nd, u
          heappush(pq, nd * n + v)
    self.dist, self.par = dist, par

  def path(self, t):
    p = []
    if self.dist[t] < self.INF:
      while t >= 0:
        p.append(t)
        t = self.par[t]
    return p[::-1]
