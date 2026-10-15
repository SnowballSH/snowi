from collections import deque

from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("flow/dinic")
Dinic = lib.Dinic


def brute_cut(n, arcs, s, t):
  """Minimum s-t cut by enumerating every source side."""
  best = None
  for mask in range(1 << n):
    if not mask >> s & 1 or mask >> t & 1:
      continue
    c = 0
    for u, v, cap, rev in arcs:
      in_u, in_v = mask >> u & 1, mask >> v & 1
      if in_u and not in_v:
        c += cap
      if in_v and not in_u:
        c += rev
    best = c if best is None else min(best, c)
  return best


def edmonds_karp(n, arcs, s, t):
  c = [[0] * n for _ in range(n)]
  for u, v, cap, rev in arcs:
    if u != v:
      c[u][v] += cap
      c[v][u] += rev
  total = 0
  while True:
    par = [-1] * n
    par[s] = s
    q = deque([s])
    while q and par[t] < 0:
      u = q.popleft()
      for v in range(n):
        if par[v] < 0 and c[u][v] > 0:
          par[v] = u
          q.append(v)
    if par[t] < 0:
      return total
    f, v = None, t
    while v != s:
      f = c[par[v]][v] if f is None else min(f, c[par[v]][v])
      v = par[v]
    v = t
    while v != s:
      c[par[v]][v] -= f
      c[v][par[v]] += f
      v = par[v]
    total += f


def certify(d, n, arcs, ids, s, t, value):
  """Flow bounds, conservation, value, and a cut of equal capacity."""
  net, cut = [0] * n, 0
  for (u, v, cap, rev), i in zip(arcs, ids):
    f = d.flow(i)
    check(-rev <= f <= cap)
    net[u] -= f
    net[v] += f
    su, sv = d.source_side(u), d.source_side(v)
    if su and not sv:
      cut += cap
      check_eq(f, cap)
    if sv and not su:
      cut += rev
      check_eq(f, -rev)
  for v in range(n):
    if v != s and v != t:
      check_eq(net[v], 0)
  check_eq(net[t], value)
  check_eq(-net[s], value)
  check(d.source_side(s) and not d.source_side(t))
  check_eq(cut, value)


def random_arcs(n, m, maxc):
  arcs = []
  for _ in range(m):
    u, v, cap = rnd(0, n - 1), rnd(0, n - 1), rnd(0, maxc)
    kind = rnd(0, 3)
    rev = cap if kind == 1 else rnd(0, maxc) if kind == 2 else 0
    arcs.append((u, v, cap, rev))
  return arcs


def build(n, arcs):
  d = Dinic(n)
  return d, [d.add_edge(*a) for a in arcs]


for _ in range(iters(2000)):
  n, m = rnd(2, 7), rnd(0, 14)
  maxc = rnd(1, 10) if rnd(0, 3) else 10**15
  arcs = random_arcs(n, m, maxc)
  s = rnd(0, n - 1)
  t = rnd(0, n - 2)
  t += t >= s
  d, ids = build(n, arcs)
  f = d.max_flow(s, t)
  check_eq(f, brute_cut(n, arcs, s, t))
  certify(d, n, arcs, ids, s, t, f)
  check_eq(d.max_flow(s, t), 0)  # already maximum
  certify(d, n, arcs, ids, s, t, f)
  d.cap[:] = d.orig  # Tweak: start over on another pair
  f = d.max_flow(t, s)
  check_eq(f, brute_cut(n, arcs, t, s))
  certify(d, n, arcs, ids, t, s, f)

for _ in range(iters(200)):
  n, m = rnd(2, 25), rnd(0, 100)
  arcs = random_arcs(n, m, rnd(1, 1000))
  d, ids = build(n, arcs)
  f = d.max_flow(0, n - 1)
  check_eq(f, edmonds_karp(n, arcs, 0, n - 1))
  certify(d, n, arcs, ids, 0, n - 1, f)

# Extremes: capacities beyond 64 bits are fine in Python.
big = 1 << 70
arcs = [(0, 1, big, 0), (1, 3, big, big), (0, 2, big, 0), (2, 3, big, 0)]
d, ids = build(4, arcs)
check_eq(d.max_flow(0, 3), 2 * big)
certify(d, 4, arcs, ids, 0, 3, 2 * big)

side = 30_000
n, s, t = 2 * side + 2, 2 * side, 2 * side + 1
arcs = [(s, i, 1, 0) for i in range(side)]
arcs += [(side + i, t, 1, 0) for i in range(side)]
arcs += [
  (rnd(0, side - 1), side + rnd(0, side - 1), 1, 0) for _ in range(90_000)
]
d, ids = build(n, arcs)
with Timer("matching network, 3e4 + 3e4 vertices, 1.5e5 edges"):
  f = d.max_flow(s, t)
check(f > side * 9 // 10)
certify(d, n, arcs, ids, s, t, f)

w = 100  # grid, undirected random capacities, left column to right
n, s, t = w * w + 2, w * w, w * w + 1
arcs = [(s, r * w, 10**9, 0) for r in range(w)]
arcs += [(r * w + w - 1, t, 10**9, 0) for r in range(w)]
for r in range(w):
  for c in range(w):
    if c + 1 < w:
      x = rnd(1, 100)
      arcs.append((r * w + c, r * w + c + 1, x, x))
    if r + 1 < w:
      x = rnd(1, 100)
      arcs.append((r * w + c, r * w + c + w, x, x))
d, ids = build(n, arcs)
with Timer("grid 100 x 100, 2e4 undirected edges"):
  f = d.max_flow(s, t)
certify(d, n, arcs, ids, s, t, f)

n = 2000
arcs = random_arcs(n, 50_000, 10**12)
d, ids = build(n, arcs)
with Timer("random network, n = 2000, m = 5e4"):
  f = d.max_flow(0, n - 1)
certify(d, n, arcs, ids, 0, n - 1, f)

n = 200_000
arcs = [(i, i + 1, rnd(5, 9), 0) for i in range(n - 1)]
arcs += [(i, i + 3, 1, 0) for i in range(n - 3)]
d, ids = build(n, arcs)
with Timer("path of 2e5 vertices (deep augmenting paths)"):
  f = d.max_flow(0, n - 1)
certify(d, n, arcs, ids, 0, n - 1, f)
# Diamond chain of 2 x 60 vertices ending in a dead end, plus a
# separate path to t of the same length. Without current-edge
# pointers the dead end is explored once per path: 2^60 times.
L = 60
n, s, t = 3 * L + 2, 3 * L, 3 * L + 1
arcs = [(s, 0, 10, 0), (s, 1, 10, 0), (s, 2 * L, 1, 0)]
for i in range(L - 1):
  for a in (2 * i, 2 * i + 1):
    arcs += [(a, 2 * i + 2, 10, 0), (a, 2 * i + 3, 10, 0)]
  arcs.append((2 * L + i, 2 * L + i + 1, 1, 0))
arcs.append((3 * L - 1, t, 1, 0))
d, ids = build(n, arcs)
check_eq(d.max_flow(s, t), 1)
certify(d, n, arcs, ids, s, t, 1)

# Broom: a long chain shared by many unit paths, all in one phase.
# Augmenting one path at a time walks the chain each time.
L = K = 50_000
n, t = L + K + 2, L + K + 1
arcs = [(i, i + 1, 10**9, 0) for i in range(L)]
for j in range(L + 1, L + K + 1):
  arcs += [(L, j, 1, 0), (j, t, 1, 0)]
d, ids = build(n, arcs)
with Timer("broom: chain of 5e4 shared by 5e4 paths") as timer:
  check_eq(d.max_flow(0, t), K)
check(timer.seconds < 3, "augmenting path by path is quadratic")
certify(d, n, arcs, ids, 0, t, K)
print("OK flow/dinic")
