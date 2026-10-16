"""
Title: Manhattan MST
Description: Minimum spanning tree of the complete graph on the
  points, with distance $|x_i - x_j| + |y_i - y_j|$: returns the
  total weight and the $n - 1$ tree edges `(i, j)`.
Precondition: Integer tuples `(x, y, ...)`, $|x|, |y| <= 2^30$, at
  most $2^20$ points; equal points are allowed (distance 0).
Time: $O(n log n)$; 1.4 s for $2 dot 10^5$ points.
Usage:
  total, tree = manhattan_mst(pts)
Tweak: Chebyshev distance $max(|Delta x|, |Delta y|)$: use the
  points $(x + y, x - y)$ and halve the total.
Example:
  >>> manhattan_mst([(0, 0), (5, 1), (1, 4), (6, 6)])
  (17, [(0, 2), (0, 1), (1, 3)])
Notes: Kruskal runs on at most $4n$ candidate edges: each point
  joined to its nearest point in four of the eight octants, one of
  each opposite pair: if an MST edge leaves $p$ into octant $O$,
  it enters $q$ from the opposite octant, and among the points of
  one octant of $q$ the nearest one is at least as close to every
  other point of that octant, so only nearest neighbours are
  needed. Per octant, a sweep by decreasing $y$ asks a Fenwick
  tree over $x - y$ for the least $x + y$. All sorts are on packed
  ints: key sorts and tuple sorts are about ten times slower in
  PyPy.
Source: The octant argument of Zhou, Shi and Feng; Fenwick sweep.
Status: stress-tested against brute force
"""


def manhattan_mst(ps):
  n, OFF, M = len(ps), 1 << 32, (1 << 20) - 1
  xs, ys = [p[0] for p in ps], [p[1] for p in ps]
  ew, eu, ev = [], [], []
  for k in range(4):
    # For each i the j with y_j >= y_i, x_j - y_j >= x_i - y_i
    # and the least x_j + y_j: decreasing y, Fenwick over x - y.
    by_d = sorted(
      (xs[i] - ys[i] + OFF) << 20 | i for i in range(n)
    )
    pos, at, m, last = [0] * n, [0] * n, 0, -1
    for r in range(n - 1, -1, -1):  # Fenwick index 1 = largest
      i, d = by_d[r] & M, by_d[r] >> 20
      m += d != last
      last, pos[i], at[i] = d, m, r
    order = sorted(
      (ys[i] + OFF) << 20 | at[i] for i in range(n)
    )
    best, who = [1 << 62] * (m + 1), [-1] * (m + 1)
    for v in reversed(order):  # ties in y: larger x - y first
      i = by_d[v & M] & M
      s, r = xs[i] + ys[i], pos[i]
      b, w, j = 1 << 62, -1, r
      while j:
        if best[j] < b:
          b, w = best[j], who[j]
        j &= j - 1
      if w >= 0:
        ew.append(b - s)
        eu.append(i)
        ev.append(w)
      while r <= m:
        if s < best[r]:
          best[r], who[r] = s, i
        r += r & -r
    if k & 1:
      xs = [-x for x in xs]
    else:
      xs, ys = ys, xs
  par, total, tree = list(range(n)), 0, []
  for key in sorted(w << 22 | e for e, w in enumerate(ew)):
    e = key & ((1 << 22) - 1)
    u, v = eu[e], ev[e]
    while par[u] != u:
      par[u] = u = par[par[u]]
    while par[v] != v:
      par[v] = v = par[par[v]]
    if u != v:
      par[u] = v
      total += key >> 22
      tree.append((eu[e], ev[e]))
  return total, tree
