"""
Title: Bellman-Ford
Description: Shortest paths from `s` when weights may be negative.
  `dist[v]` is `INF` if `v` cannot be reached and `-INF` if paths
  to `v` can be arbitrarily short (`v` is reachable from a negative
  cycle that `s` reaches), with `INF = BellmanFord.INF`. For finite
  `dist[v]`, `par[v]` is the vertex before `v` on a shortest path
  (-1 for `s`).
Precondition: Edges `(a, b, w)` are directed `a -> b`; vertices
  0-indexed. True distances are below `INF` $= 2^62$ in absolute
  value.
Time: $O(n m)$; stops as soon as a round changes nothing.
Usage:
  bf = BellmanFord(n, edges, s)  # edges: (a, b, w), a -> b
  if bf.dist[v] == -BellmanFord.INF: ...  # no shortest path
Tweak: Print a negative cycle (any one if `dist` starts all 0):
  relax plainly for $n$ rounds (no `-INF` marking); if round $n$
  still relaxes some `b`, follow `par` from `b` $n$ times, then
  once around. Real weights (arbitrage): relax only if
  `da + w < dist[b] - 1e-9`, or rounding fakes negative cycles.
Example:
  >>> INF = BellmanFord.INF
  >>> edges = [(0, 1, 4), (0, 2, 5), (2, 1, -3), (1, 3, 1)]
  >>> bf = BellmanFord(5, edges, 0)
  >>> bf.dist[:4], bf.par[:4], bf.dist[4] == INF
  ([0, 2, 5, 3], [-1, 2, 0, 1], True)
  >>> edges += [(3, 4, 0), (4, 3, -1)]  # negative cycle 3, 4
  >>> bf = BellmanFord(5, edges, 0)
  >>> bf.dist[2], bf.dist[3] == bf.dist[4] == -INF
  (5, True)
  >>> # difference constraints x[b] - x[a] <= w, source n = 3
  >>> cons = [(0, 1, 3), (1, 2, -2), (0, 2, 2), (2, 0, -1)]
  >>> bf = BellmanFord(4, cons + [(3, v, 0) for v in range(3)], 3)
  >>> bf.dist[:3]  # one solution: x = dist
  [-3, 0, -2]
Notes: Use it when weights are negative: arbitrage (weights
  $-log(r)$, see Tweak), difference constraints (Example: feasible iff no
  `-INF`), "is there a cycle that gains money", or potentials for
  Johnson's all-pairs algorithm (then Dijkstra on
  $w(u, v) + d(u) - d(v) >= 0$). An undirected negative edge is
  already a negative cycle. SPFA (queue-based) has the same worst
  case and is often targeted by anti-tests. PyPy checks about
  $2 dot 10^7$ edges per second, so $n m <= 10^7$ is safe.
Source: kactl BellmanFord (CC0), simplified (no edge sorting).
Status: stress-tested against walk-length dynamic programming
"""


class BellmanFord:
  INF = 1 << 62

  def __init__(self, n, edges, s):
    INF = self.INF
    dist, par = [INF] * n, [-1] * n
    dist[s] = 0
    # Rounds >= n - 1 only mark: any change there means -INF.
    for i in range(2 * n):
      changed = False
      for a, b, w in edges:
        da = dist[a]
        if da == INF or dist[b] == -INF:
          continue
        if da == -INF or da + w < dist[b]:
          dist[b] = da + w if i < n - 1 else -INF
          par[b] = a
          changed = True
      if not changed:
        break
    self.dist, self.par = dist, par
