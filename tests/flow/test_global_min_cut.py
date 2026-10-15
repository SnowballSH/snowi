from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("flow/global-min-cut", "flow/dinic")
global_min_cut, Dinic = lib.global_min_cut, lib.Dinic


def cut_weight(w, inside):
  n = len(w)
  return sum(
    w[u][v] for u in range(n) for v in range(n) if inside[u] and not inside[v]
  )


def brute(w):
  """Minimum over every non-empty proper subset."""
  n = len(w)
  return min(
    cut_weight(w, [mask >> v & 1 for v in range(n)])
    for mask in range(1, (1 << n) - 1)
  )


def by_flows(w):
  """Minimum over t of the max flow from 0 to t."""
  n, best = len(w), None
  for t in range(1, n):
    d = Dinic(n)
    for u in range(n):
      for v in range(u):
        if w[u][v]:
          d.add_edge(u, v, w[u][v], w[u][v])
    f = d.max_flow(0, t)
    best = f if best is None else min(best, f)
  return best


def solve(w):
  """The cut, after checking that side is a proper part whose
  weight is the cut and that w is unchanged."""
  n, before = len(w), [row[:] for row in w]
  cut, side = global_min_cut(w)
  check_eq(w, before)
  check(0 < len(side) < n)
  inside = [False] * n
  for v in side:
    check(0 <= v < n and not inside[v])
    inside[v] = True
  check_eq(cut_weight(w, inside), cut)
  return cut


def random_graph(n, density, maxw):
  w = [[0] * n for _ in range(n)]
  for u in range(n):
    w[u][u] = rnd(0, maxw)  # ignored
    for v in range(u):
      if rnd(1, 100) <= density:
        w[u][v] = w[v][u] = rnd(0, maxw)
  return w


check_eq(solve([[0, 0], [0, 0]]), 0)
check_eq(solve([[5, 7], [7, 9]]), 7)
for _ in range(iters(1500)):
  n, density = rnd(2, 9), rnd(10, 100)
  maxw = rnd(1, 10) if rnd(0, 3) else 10**15
  w = random_graph(n, density, maxw)
  check_eq(solve(w), brute(w))

for _ in range(iters(40)):
  n = rnd(10, 30)
  w = random_graph(n, rnd(5, 100), rnd(1, 1000))
  check_eq(solve(w), by_flows(w))

# Large: two dense halves joined by a few light edges.
n = 500
w = random_graph(n, 50, 1000)
half = [rnd(0, 1) for _ in range(n)]
for u in range(n):
  for v in range(n):
    if half[u] != half[v]:
      w[u][v] = 0
planted = 0
for _ in range(10):
  u, v = rnd(0, n - 1), rnd(0, n - 1)
  if half[u] != half[v] and not w[u][v]:
    w[u][v] = w[v][u] = 1
    planted += 1
with Timer("n = 500, dense"):
  cut, side = global_min_cut(w)
check_eq(cut, planted)
check_eq(sorted(side), [v for v in range(n) if half[v] == half[side[0]]])
print("OK flow/global-min-cut")
