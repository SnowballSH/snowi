"""
Title: Floyd-Warshall
Description: All-pairs shortest paths; weights may be negative.
  On input `d[i][j]` is the weight of edge `i -> j`, `INF` if there
  is none. On output it is the distance: `INF` if `j` is unreachable
  from `i`, `-INF` if a negative cycle on the way makes it
  arbitrarily small.
Precondition: `d` is $n times n$; keep the lightest of parallel
  edges. `INF` $> n dot max|w|$, at most $2^62$ to keep PyPy on
  machine ints.
Time: $O(n^3)$; PyPy needs about 0.3 s for $n = 400$.
Usage:
  d = [[INF] * n for _ in range(n)]
  d[u][v] = min(d[u][v], w)
  floyd_warshall(d, INF)
Tweak: Transitive closure: rows as int bitmasks,
  `if r[i] >> k & 1: r[i] |= r[k]` in the same loop order.
Example:
  >>> INF = 1 << 62
  >>> d = [[INF] * 4 for _ in range(4)]
  >>> edges = [(0, 1, 2), (1, 2, -1), (0, 2, 4), (2, 3, 1)]
  >>> for u, v, w in edges: d[u][v] = w
  >>> floyd_warshall(d, INF); d[0], d[3][0] == INF
  ([0, 2, 1, 2], True)
  >>> d[2][1] = -1  # cycle 1 -> 2 -> 1 of weight -2
  >>> floyd_warshall(d, INF); d[0][3] == -INF, d[3][3]
  (True, 0)
  >>> # Add edge a -> b of weight w in O(n^2), if there is no
  >>> # negative cycle before or after: d[b][a] + w >= 0 (with
  >>> # d[b][a] == INF fine), else rerun. Skip the sentinel.
  >>> d = [[INF] * 3 for _ in range(3)]; d[0][1] = 5
  >>> floyd_warshall(d, INF); a, b, w = 1, 2, -1
  >>> for i in range(3):
  ...   for j in range(3):
  ...     if d[i][a] != INF and d[b][j] != INF:
  ...       d[i][j] = min(d[i][j], d[i][a] + w + d[b][j])
  >>> d[0], d[1][2], d[1][0] == INF
  ([0, 5, 4], -1, True)
Notes: Use it for $n <= 500$ when many pairs are queried, for
  transitive closure, or when edges are added one at a time
  ($O(n^2)$ per edge while no negative cycle appears, Example). The
  loop over `k` must be the outermost. Processing vertices `k` in a chosen order
  (sorted by some value) answers "shortest path using only
  intermediate vertices with value $<= x$" offline. Shortest cycle
  in a directed graph: skip the diagonal reset, then take the
  minimum of `d[i][i]`. Negative cycles make unclamped values grow
  exponentially, so sums are clamped at `-INF`.
Source: kactl FloydWarshall (CC0).
Status: stress-tested against walk-length dynamic programming
"""


def floyd_warshall(d, INF):
  n = len(d)
  for i in range(n):
    d[i][i] = min(d[i][i], 0)
  for k in range(n):
    dk = d[k]
    for i in range(n):
      di = d[i]
      dik = di[k]
      if dik == INF:
        continue
      for j in range(n):
        nd = dik + dk[j]
        # dk[j] is tested only on improvement: faster in PyPy.
        if nd < di[j] and dk[j] != INF:
          di[j] = nd if nd > -INF else -INF
  for k in range(n):
    if d[k][k] < 0:
      dk = d[k]
      for di in d:
        if di[k] != INF:
          for j in range(n):
            if dk[j] != INF:
              di[j] = -INF
