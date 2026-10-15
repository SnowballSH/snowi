from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("graph/bellman-ford")
BellmanFord = lib.BellmanFord
INF = BellmanFord.INF


def brute(n, edges, s):
  """walk[k][x][v]: lightest walk with exactly k edges from x to v."""
  walk = [[[None] * n for _ in range(n)] for _ in range(n + 1)]
  for x in range(n):
    walk[0][x][x] = 0
  for k in range(n):
    for x in range(n):
      for a, b, w in edges:
        if walk[k][x][a] is not None:
          cand = walk[k][x][a] + w
          old = walk[k + 1][x][b]
          walk[k + 1][x][b] = cand if old is None else min(old, cand)

  def reach(x, v):
    return any(walk[k][x][v] is not None for k in range(n))

  on_neg = [
    any(
      walk[k][x][x] is not None and walk[k][x][x] < 0 for k in range(1, n + 1)
    )
    for x in range(n)
  ]
  d = []
  for v in range(n):
    if any(on_neg[y] and reach(s, y) and reach(y, v) for y in range(n)):
      d.append(-INF)
    else:
      vals = [walk[k][s][v] for k in range(n) if walk[k][s][v] is not None]
      d.append(min(vals) if vals else INF)
  return d


def check_par(n, edges, s, bf):
  """Finite vertices: par edges are tight and lead back to s."""
  dist, par = bf.dist, bf.par
  tight = [False] * n
  for a, b, w in edges:
    if abs(dist[b]) < INF and a == par[b] and dist[a] + w == dist[b]:
      tight[b] = True
  for v in range(n):
    if abs(dist[v]) < INF:
      u, steps = v, 0
      while u != s:
        check(steps < n and tight[u])
        u, steps = par[u], steps + 1
  if abs(dist[s]) < INF:
    check_eq(par[s], -1)


kinds = [0, 0, 0]
for _ in range(iters(2000)):
  n, m = rnd(1, 7), rnd(0, 14)
  lo, hi = (-5 if rnd(0, 2) else 0), 10
  if rnd(0, 4) == 0:
    lo, hi = -(10**12), 10**12
  edges = [(rnd(0, n - 1), rnd(0, n - 1), rnd(lo, hi)) for _ in range(m)]
  s = rnd(0, n - 1)
  bf = BellmanFord(n, edges, s)
  want = brute(n, edges, s)
  check_eq(bf.dist, want)
  check_par(n, edges, s, bf)
  for d in want:
    kinds[(d == INF) + 2 * (d == -INF)] += 1
check(min(kinds) > sum(kinds) // 10, f"unbalanced cases {kinds}")

# A cycle listed backwards with total weight -1: reaching every vertex
# takes n - 1 rounds, marking all of them -INF n - 1 more.
n = 300
edges = [(v - 1, v, -1) for v in range(n - 1, 0, -1)] + [(n - 1, 0, n - 2)]
check_eq(BellmanFord(n, edges, 0).dist, [-INF] * n)

# Speed: a chain listed backwards needs n - 1 full rounds.
n, m = 2000, 10000
edges = [(v - 1, v, -1) for v in range(n - 1, 0, -1)]
while len(edges) < m:
  a = rnd(0, n - 1)
  edges.append((a, rnd(0, a), 10**9))
with Timer("n = 2000, m = 1e4, n rounds: 2e7 relaxations"):
  bf = BellmanFord(n, edges, 0)
check_eq(bf.dist, [-v for v in range(n)])
check_par(n, edges, 0, bf)
print("OK graph/bellman-ford")
