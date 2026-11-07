from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("graph/dijkstra")
Dijkstra = lib.Dijkstra
INF = Dijkstra.INF


def brute(adj, srcs):
  """Relax every edge n times (Bellman-Ford without early exit)."""
  n = len(adj)
  d = [INF] * n
  for s in srcs:
    d[s] = 0
  for _ in range(n):
    for u in range(n):
      if d[u] < INF:
        for v, w in adj[u]:
          d[v] = min(d[v], d[u] + w)
  return d


def certify(adj, srcs, dj):
  """dist is feasible and every reached vertex has a tight path."""
  n = len(adj)
  for s in srcs:
    check_eq(dj.dist[s], 0)
  for u in range(n):
    if dj.dist[u] < INF:
      for v, w in adj[u]:
        check(dj.dist[v] <= dj.dist[u] + w)
  for t in range(n):
    p = dj.path(t)
    if dj.dist[t] == INF:
      check_eq(p, [])
      check_eq(dj.par[t], -1)
      continue
    check(p and p[-1] == t and p[0] in srcs)
    check_eq(dj.par[p[0]], -1)
    for a, b in zip(p, p[1:]):
      check(any(v == b and dj.dist[a] + w == dj.dist[b] for v, w in adj[a]))
      check_eq(dj.par[b], a)


for _ in range(iters(3000)):
  n, m = rnd(1, 8), rnd(0, 20)
  maxw = 3 if rnd(0, 1) else 10**12
  adj = [[] for _ in range(n)]
  for _ in range(m):
    adj[rnd(0, n - 1)].append((rnd(0, n - 1), rnd(0, maxw)))
  srcs = [rnd(0, n - 1) for _ in range(rnd(1, 3))]
  dj = Dijkstra(adj, srcs)
  check_eq(dj.dist, brute(adj, srcs))
  certify(adj, srcs, dj)
# Largest distance allowed by the Precondition: INF - 1.
big = [[(1, INF // 2)], [(2, INF // 2 - 1)], []]
check_eq(Dijkstra(big, [0]).dist, [0, INF // 2, INF - 1])
check_eq(Dijkstra(big, [0]).path(2), [0, 1, 2])

n, m = 200_000, 400_000
adj = [[] for _ in range(n)]
for _ in range(m):
  u, v, w = rnd(0, n - 1), rnd(0, n - 1), rnd(0, 10**9)
  adj[u].append((v, w))
  adj[v].append((u, w))
with Timer("n = 2e5, 4e5 undirected edges"):
  dj = Dijkstra(adj, [0])
certify(adj, [0], dj)
line = [[(i + 1, 10**12)] for i in range(n - 1)] + [[]]
with Timer("path of 2e5 vertices, distances up to 2e17"):
  dl = Dijkstra(line, [0])
check_eq(dl.dist[n - 1], (n - 1) * 10**12)
check_eq(len(dl.path(n - 1)), n)
print("OK graph/dijkstra")
