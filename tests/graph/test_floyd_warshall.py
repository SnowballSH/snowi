from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("graph/floyd-warshall")
INF = 1 << 62


def brute(w):
  """walk[k][x][v]: lightest walk with exactly k edges from x to v."""
  n = len(w)
  walk = [[[INF] * n for _ in range(n)] for _ in range(n + 1)]
  for x in range(n):
    walk[0][x][x] = 0
  for k in range(n):
    for x in range(n):
      for a in range(n):
        for b in range(n):
          if walk[k][x][a] != INF and w[a][b] != INF:
            cand = walk[k][x][a] + w[a][b]
            walk[k + 1][x][b] = min(walk[k + 1][x][b], cand)

  def reach(x, v):
    return any(walk[k][x][v] != INF for k in range(n))

  on_neg = [any(walk[k][x][x] < 0 for k in range(1, n + 1)) for x in range(n)]
  d = [
    [min(walk[k][i][j] for k in range(n)) for j in range(n)] for i in range(n)
  ]
  for i in range(n):
    for j in range(n):
      if any(on_neg[y] and reach(i, y) and reach(y, j) for y in range(n)):
        d[i][j] = -INF
  return d


def add_edge(d, a, b, w):
  """The Example's O(n^2) update for a new edge a -> b of weight w."""
  for di in d:
    for j in range(len(d)):
      if di[a] != INF and d[b][j] != INF:
        di[j] = min(di[j], di[a] + w + d[b][j])


finite = neg = added = 0
for _ in range(iters(1500)):
  n = rnd(1, 6)
  lo, hi = (-4 if rnd(0, 2) else 0), 10
  if rnd(0, 4) == 0:
    lo, hi = -(10**12), 10**12
  w = [[INF] * n for _ in range(n)]
  for _ in range(rnd(0, 12)):
    a, b = rnd(0, n - 1), rnd(0, n - 1)
    w[a][b] = min(w[a][b], rnd(lo, hi))
  d = [row[:] for row in w]
  lib.floyd_warshall(d, INF)
  want = brute(w)
  check_eq(d, want)
  finite += sum(abs(x) < INF for row in want for x in row)
  neg += sum(x == -INF for row in want for x in row)
  # A new edge that closes no negative cycle: O(n^2) update.
  a, b, x = rnd(0, n - 1), rnd(0, n - 1), rnd(lo, hi)
  if any(-INF in row for row in d) or d[b][a] + x < 0:
    continue
  add_edge(d, a, b, x)
  added += 1
  w[a][b] = min(w[a][b], x)
  lib.floyd_warshall(w, INF)
  check_eq(d, w)
check(min(finite, neg) > (finite + neg) // 10, f"{finite} {neg}")
check(added > iters(1500) // 4, f"{added} edges added")
# Unreachable pairs stay unreachable (no INF - 1).
for x in (4, -1):
  d = [[INF] * 3 for _ in range(3)]
  d[0][1] = 5
  lib.floyd_warshall(d, INF)
  add_edge(d, 1, 2, x)
  check_eq(d, [[0, 5, 5 + x], [INF, 0, x], [INF, INF, 0]])

# Heavy negative cycles everywhere: clamped, all -INF.
n = 60
d = [[-(10**15)] * n for _ in range(n)]
lib.floyd_warshall(d, INF)
check_eq(d, [[-INF] * n for _ in range(n)])
# Another sentinel works too.
d = [[10**18, -1], [5, 10**18]]
lib.floyd_warshall(d, 10**18)
check_eq(d, [[0, -1], [5, 0]])

# Large, positive weights: certify against the edges directly.
n = 400
w = [[INF] * n for _ in range(n)]
for _ in range(8000):
  w[rnd(0, n - 1)][rnd(0, n - 1)] = rnd(1, 10**9)
d = [row[:] for row in w]
with Timer("n = 400"):
  lib.floyd_warshall(d, INF)
in_edges = [
  [(k, w[k][j]) for k in range(n) if w[k][j] != INF] for j in range(n)
]
for i in range(n):
  di = d[i]
  check_eq(di[i], 0)
  for j in range(n):
    if j != i:
      cands = [di[k] + x for k, x in in_edges[j] if di[k] != INF]
      check_eq(di[j], min(cands) if cands else INF)
print("OK graph/floyd-warshall")
