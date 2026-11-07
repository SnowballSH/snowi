import random
from itertools import product

from testlib import CONTENT, Timer, check, check_eq, iters, load, rnd

lib = load("flow/min-cost-flow")
MinCostFlow = lib.MinCostFlow
BIG = 1 << 62


def brute_exhaustive(n, arcs, s, t, limit):
  """Every integer flow: (min(limit, max value), min cost there)."""
  best = {}
  for f in product(*[range(a[2] + 1) for a in arcs]):
    net, c = [0] * n, 0
    for (u, v, _, w), x in zip(arcs, f):
      net[u] -= x
      net[v] += x
      c += x * w
    if net[t] >= 0 and all(net[v] == 0 for v in range(n) if v not in (s, t)):
      best[net[t]] = min(best.get(net[t], c), c)
  value = min(limit, max(best))
  return value, best[value]


def brute_ssp(n, arcs, s, t, limit):
  """Successive shortest paths with plain Bellman-Ford."""
  frm, to, cap, cost = [], [], [], []
  for u, v, c, w in arcs:
    frm += [u, v]
    to += [v, u]
    cap += [c, 0]
    cost += [w, -w]
  flow = total = 0
  while flow < limit:
    d, pe = [BIG] * n, [-1] * n
    d[s] = 0
    for _ in range(n):
      for e in range(len(to)):
        if cap[e] > 0 and d[frm[e]] < BIG and d[frm[e]] + cost[e] < d[to[e]]:
          d[to[e]] = d[frm[e]] + cost[e]
          pe[to[e]] = e
    if d[t] == BIG:
      break
    push, v = limit - flow, t
    while v != s:
      push = min(push, cap[pe[v]])
      v = frm[pe[v]]
    v = t
    while v != s:
      cap[pe[v]] -= push
      cap[pe[v] ^ 1] += push
      v = frm[pe[v]]
    flow += push
    total += push * d[t]
  return flow, total


def certify(mcf, n, arcs, ids, s, t, limit, f, c):
  """Feasible, right value and cost, no negative residual cycle,
  and no s-t residual path when below the limit."""
  net, cost, res = [0] * n, 0, []
  for (u, v, cap, w), i in zip(arcs, ids):
    x = mcf.flow(i)
    check(0 <= x <= cap)
    net[u] -= x
    net[v] += x
    cost += x * w
    if x < cap:
      res.append((u, v, w))
    if x > 0:
      res.append((v, u, -w))
  for v in range(n):
    if v != s and v != t:
      check_eq(net[v], 0)
  check_eq(net[t], f)
  check_eq(cost, c)
  d = [0] * n
  changed = True
  for _ in range(n + 1):
    changed = False
    for u, v, w in res:
      if d[u] + w < d[v]:
        d[v] = d[u] + w
        changed = True
    if not changed:
      break
  check(not changed, "negative residual cycle")
  if f < limit:
    seen = [False] * n
    seen[s] = True
    grew = True
    while grew:
      grew = False
      for u, v, _ in res:
        if seen[u] and not seen[v]:
          seen[v] = grew = True
    check(not seen[t], "augmenting path left")


def random_arcs(n, m, maxcap, maxslack):
  """Costs phi[v] - phi[u] + slack: negative edges, no negative cycle."""
  phi = [rnd(-5, 5) for _ in range(n)]
  arcs = []
  for _ in range(m):
    u, v = rnd(0, n - 1), rnd(0, n - 1)
    arcs.append((u, v, rnd(0, maxcap), phi[v] - phi[u] + rnd(0, maxslack)))
  return arcs


def build(n, arcs):
  mcf = MinCostFlow(n)
  return mcf, [mcf.add_edge(*a) for a in arcs]


def brute_any(n, arcs, s, t):
  """Least cost over every integer s-t flow of any value."""
  best = 0
  for f in product(*[range(a[2] + 1) for a in arcs]):
    net, c = [0] * n, 0
    for (u, v, _, w), x in zip(arcs, f):
      net[u] -= x
      net[v] += x
      c += x * w
    inner = all(net[v] == 0 for v in range(n) if v not in (s, t))
    if net[t] >= 0 and inner:
      best = min(best, c)
  return best


def cheapest_any_size():
  """MinCostFlow with the Tweak applied literally."""
  code = (CONTENT / "flow/min-cost-flow.py").read_text().split('"""', 2)[2]
  old = "        pot[v] += d[v] if d[v] < dt else dt\n"
  check_eq(code.count(old), 1)
  new = old + "      if pot[t] - pot[s] >= 0:\n        break\n"
  ns = {}
  exec(code.replace(old, new), ns)
  return ns["MinCostFlow"]


Cheapest = cheapest_any_size()
for _ in range(iters(700)):  # Tweak: cheapest flow of any size
  n, m = rnd(2, 4), rnd(0, 6)
  arcs = random_arcs(n, m, 2, 4)
  s = rnd(0, n - 1)
  t = rnd(0, n - 2)
  t += t >= s
  mcf = Cheapest(n)
  for a in arcs:
    mcf.add_edge(*a)
  check_eq(mcf.max_flow(s, t)[1], brute_any(n, arcs, s, t))

INF = float("inf")
for _ in range(iters(700)):
  n, m = rnd(2, 4), rnd(0, 6)
  arcs = random_arcs(n, m, 2, 0 if rnd(0, 1) else 4)
  s = rnd(0, n - 1)
  t = rnd(0, n - 2)
  t += t >= s
  limit = INF if rnd(0, 2) else rnd(0, 4)
  mcf, ids = build(n, arcs)
  f, c = mcf.max_flow(s, t, limit) if limit != INF else mcf.max_flow(s, t)
  check_eq((f, c), brute_exhaustive(n, arcs, s, t, limit))
  certify(mcf, n, arcs, ids, s, t, limit, f, c)

for _ in range(iters(700)):
  n, m = rnd(2, 12), rnd(0, 30)
  arcs = random_arcs(n, m, rnd(1, 6), rnd(0, 6))
  limit = INF if rnd(0, 2) else rnd(0, 10)
  mcf, ids = build(n, arcs)
  f, c = mcf.max_flow(0, n - 1, limit)
  check_eq((f, c), brute_ssp(n, arcs, 0, n - 1, limit))
  certify(mcf, n, arcs, ids, 0, n - 1, limit, f, c)

for _ in range(iters(700)):  # costs at the documented bound
  n, m = rnd(2, 8), rnd(0, 7)
  C = (1 << 62) // (8 * n)
  phi = [rnd(-C // 4, C // 4) for _ in range(n)]
  arcs = []
  for _ in range(m):
    u, v = rnd(0, n - 1), rnd(0, n - 1)
    arcs.append((u, v, rnd(0, 1), phi[v] - phi[u] + rnd(0, C // 2)))
  mcf, ids = build(n, arcs)
  f, c = mcf.max_flow(0, n - 1)
  check_eq((f, c), brute_ssp(n, arcs, 0, n - 1, INF))
  certify(mcf, n, arcs, ids, 0, n - 1, INF, f, c)

# Extremes: total cost 4e18, negative edge on the cheap path.
C, K = 10**9, 2 * 10**9
arcs = [(0, 1, K, C), (1, 3, K, C), (0, 2, K, -C), (2, 3, K, C), (1, 2, K, -C)]
mcf, ids = build(4, arcs)
f, c = mcf.max_flow(0, 3)
check_eq((f, c), (2 * K, 2 * K * C))
certify(mcf, 4, arcs, ids, 0, 3, INF, f, c)

# Planted assignment: rows -> columns with tight duals.
k = 200
n, s, t = 2 * k + 2, 2 * k, 2 * k + 1
perm = list(range(k))
random.Random(rnd(0, 10**9)).shuffle(perm)
u = [rnd(-(10**6), 10**6) for _ in range(k)]
v = [rnd(-(10**6), 10**6) for _ in range(k)]
arcs = [(s, i, 1, 0) for i in range(k)] + [(k + i, t, 1, 0) for i in range(k)]
for i in range(k):
  for j in range(k):
    slack = 0 if perm[i] == j else rnd(1, 10**6)
    arcs.append((i, k + j, 1, u[i] + v[j] + slack))
mcf, ids = build(n, arcs)
with Timer("assignment 200 x 200 (4e4 edges, 200 Dijkstras)"):
  f, c = mcf.max_flow(s, t)
check_eq((f, c), (k, sum(u) + sum(v)))
certify(mcf, n, arcs, ids, s, t, INF, f, c)

n = 2000
arcs = random_arcs(n, 20_000, 20, 100)
mcf, ids = build(n, arcs)
with Timer("random network n = 2000, m = 2e4, negative costs"):
  f, c = mcf.max_flow(0, n - 1)
certify(mcf, n, arcs, ids, 0, n - 1, INF, f, c)

# Comb s -> a_i (cost i) -> x (cost -2i) -> tail -> t. Without the
# Bellman-Ford potentials the first Dijkstra improves x k times and
# rescans the tail each time: k * L = 2.5e9 steps.
k = L = 50_000
x, t = k + 1, k + L + 2
mcf = MinCostFlow(t + 1)
for i in range(1, k + 1):
  mcf.add_edge(0, i, 1, i)
  mcf.add_edge(i, x, 1, -2 * i)
for j in range(L):
  mcf.add_edge(x + j, x + j + 1, k, 0)
mcf.add_edge(x + L, t, 1, 0)
with Timer("comb with negative costs, 1e5 vertices"):
  check_eq(mcf.max_flow(0, t), (1, -k))
# Sparse assignment: 10 random edges per row, F = 1000 paths.
k = 1000
n, s, t = 2 * k + 2, 2 * k, 2 * k + 1
arcs = [(s, i, 1, 0) for i in range(k)] + [(k + i, t, 1, 0) for i in range(k)]
for i in range(k):
  arcs += [(i, k + rnd(0, k - 1), 1, rnd(0, 10**6)) for _ in range(10)]
mcf, ids = build(n, arcs)
with Timer("sparse assignment, k = 1000, 1e4 edges"):
  f, c = mcf.max_flow(s, t)
certify(mcf, n, arcs, ids, s, t, INF, f, c)

# K cheap paths s -> x_i -> t next to a costly side chain with 4e5
# edges. Dijkstra must stop at t, not scan the chain K times.
K, L, y = 600, 20_000, 602
mcf = MinCostFlow(K + L + 2)
for i in range(K):
  mcf.add_edge(0, 2 + i, 1, i)
  mcf.add_edge(2 + i, 1, 1, 0)
mcf.add_edge(0, y, 1, 10**6)
for j in range(L - 1):
  for p in range(20):
    mcf.add_edge(y + j, y + j + 1, 1, 1 + p)
with Timer("600 paths next to a side chain of 4e5 edges") as timer:
  check_eq(mcf.max_flow(0, 1), (K, K * (K - 1) // 2))
check(timer.seconds < 3, "Dijkstra does not stop at t")  # else ~12 s
print("OK flow/min-cost-flow")
