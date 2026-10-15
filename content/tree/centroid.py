"""
Title: Centroid decomposition
Description: Repeatedly removes a centroid (a vertex whose
  removal leaves parts of at most half the piece) and recurses
  into the parts. `par[v]` is the parent of `v` in the centroid
  tree (-1 for the first centroid), `level[v]` its depth there,
  and `dist[k][v]` the distance in the tree from `v` to its
  centroid ancestor at level `k` (for `k <= level[v]`). The
  path between `u` and `v` passes through their lowest common
  centroid ancestor `c`, so `dist(u, v) = dist(u, c) + dist(c, v)`.
Precondition: `adj` holds both directions of the $n - 1$ edges of
  a tree on vertices `0..n-1`, $n >= 1$.
Time: $O(n log n)$ time and memory, iterative; at most
  $floor(log_2 n) + 1$ levels. PyPy: about 1 s for
  $n = 2 dot 10^5$.
Usage:
  cd = Centroid(adj)  # nearest marked vertex, O(log n) per call:
  INF = 1 << 30
  best = [INF] * n  # best[c]: nearest marked in c's piece
  def mark(v):
    c = v
    while c >= 0:
      best[c] = min(best[c], cd.dist[cd.level[c]][v])
      c = cd.par[c]
  def nearest(v):  # INF if none
    c, r = v, INF
    while c >= 0:
      r = min(r, best[c] + cd.dist[cd.level[c]][v])
      c = cd.par[c]
    return r
Tweak: Weighted edges: keep `adj` as plain vertex lists and pass
  a parallel `wt` (`wt[u][i]` weighs the edge to `adj[u][i]`);
  set `gw = [wt[u] for u in order]` next to `g`, and in the piece
  BFS loop over `zip(g[u], gw[u])` with `d[v] = d[u] + w`.
  Sum of distances from `v` to a changing marked set: per
  centroid `c` keep `cnt[c]` and `sum[c]`, the number and
  distance sum (to `c`) of marked vertices in its piece, and
  `sub[c]`, their distance sum to `par[c]`; walking up from `v`
  with `ch` the previous centroid, add `sum[c] - sub[ch] +
  (cnt[c] - cnt[ch]) * d(v, c)` (without the `ch` terms at `c =
  v`).
Example:
  >>> adj = [[1], [0, 2], [1, 3], [2, 4], [3, 5], [4, 6], [5]]
  >>> cd = Centroid(adj)  # the path 0-1-2-3-4-5-6
  >>> cd.par, cd.level
  ([1, 3, 1, -1, 5, 3, 5], [2, 1, 2, 0, 2, 1, 2])
  >>> [cd.dist[k][6] for k in range(cd.level[6] + 1)]  # 3, 5, 6
  [3, 1, 0]
  >>> # pairs {u, v} with dist(u, v) == K: cnt[c][d] counts the
  >>> # added vertices of c's piece at distance d from c, and
  >>> # sub[c][d] the same from par[c] (pairs inside one part of
  >>> # par[c] must not count there)
  >>> n, K = 7, 2
  >>> size = [0] * n  # piece sizes
  >>> for v in range(n):
  ...   c = v
  ...   while c >= 0: size[c], c = size[c] + 1, cd.par[c]
  >>> cnt = [[0] * size[c] for c in range(n)]
  >>> sub = [[0] * (size[c] + 1) for c in range(n)]
  >>> pairs = 0
  >>> for v in range(n):
  ...   c, ch = v, -1
  ...   while c >= 0:
  ...     d = cd.dist[cd.level[c]][v]
  ...     e = K - d
  ...     if 0 <= e < size[c]: pairs += cnt[c][e]
  ...     if ch >= 0 and 0 <= e <= size[ch]: pairs -= sub[ch][e]
  ...     cnt[c][d] += 1
  ...     if ch >= 0: sub[ch][d] += 1
  ...     c, ch = cd.par[c], c
  >>> pairs
  5
Notes: Every path is counted at exactly one centroid: the first
  one removed on it, the lowest common centroid ancestor of its
  ends. Each vertex lies in at most $floor(log_2 n) + 1$ pieces,
  so $O(s)$ work per piece of size $s$ costs $O(n log n)$ in
  total, and $O(s log s)$ work (sorting depths, an FFT of the
  depth counts) $O(n log^2 n)$. For nearest or minimum queries
  (Usage), meeting a pair again at a higher common centroid only
  overestimates their distance, so no correction is needed; for
  counts and sums it is (the Example's `sub`). Cues: "count
  paths of length `K` or at most `K`", "number of pairs whose
  path has a property", "distance to the nearest red vertex with
  recolouring", "sum of distances to a changing set", and DP
  transitions $min_v ("val"_v + d(u, v))$ over active `v` (keep
  $"val"_v + d(v, c)$ in a multiset per centroid `c`). Size
  arrays per centroid by the piece size, never by `n`. Paths of
  length at most `K`: sort each piece's depths and count pairs
  with two pointers, minus the same inside each part.
Source: Classic; Codeforces 342E pattern.
Status: stress-tested against brute force
"""


class Centroid:
  def __init__(self, adj):
    n = len(adj)
    # Relabel in BFS order from 0: memory locality, 2x faster.
    order, new = [0], [0] + [-1] * (n - 1)
    for u in order:
      for v in adj[u]:
        if new[v] < 0:
          new[v] = len(order)
          order.append(v)
    g = [[new[v] for v in adj[u]] for u in order]
    par, level, dist = [-1] * n, [-1] * n, []
    bpar, size = [-1] * n, [1] * n
    for v in range(n - 1, 0, -1):  # BFS tree from 0: sizes
      bpar[v] = min(g[v])
      size[bpar[v]] += size[v]

    def centroid(s):  # walk towards the big part
      half, p = size[s] // 2, -1
      while p != s:
        p = s
        for v in g[p]:
          if v != bpar[p] and level[v] < 0 and size[v] > half:
            s = v
      return s

    st = [centroid(0), -1]  # flat pairs (centroid, parent)
    while st:
      p, c = st.pop(), st.pop()
      par[c], level[c] = p, 0 if p < 0 else level[p] + 1
      if len(dist) == level[c]:
        dist.append([0] * n)
      d, q = dist[level[c]], [c]
      bpar[c], d[c] = -1, 0
      for u in q:  # BFS over c's piece: parents, sizes, depths
        du = d[u] + 1
        for v in g[u]:
          if v != bpar[u] and level[v] < 0:
            bpar[v], d[v] = u, du
            q.append(v)
      for v in q:
        size[v] = 1
      for i in range(len(q) - 1, 0, -1):
        size[bpar[q[i]]] += size[q[i]]
      for v in g[c]:
        if level[v] < 0:
          st.append(centroid(v))
          st.append(c)
    order.append(-1)  # order[-1] = -1 keeps the root's parent
    self.par = [order[par[i]] for i in new]
    self.level = [level[i] for i in new]
    self.dist = [[row[i] for i in new] for row in dist]
