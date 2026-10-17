from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("graph/edge-colouring")


def check_colouring(L, R, edges):
  """Proper, with exactly the colours [0, Delta): Konig is tight."""
  col = lib.edge_colouring(L, R, edges)
  check_eq(len(col), len(edges))
  deg = [0] * (L + R)
  for u, v in edges:
    deg[u] += 1
    deg[L + v] += 1
  D = max(deg, default=0)
  seen = set()
  for (u, v), c in zip(edges, col):
    check(0 <= c < D)
    check((u, c) not in seen and (L + v, c) not in seen)
    seen.add((u, c))
    seen.add((L + v, c))


check_colouring(0, 0, [])
check_colouring(3, 2, [])
check_colouring(1, 1, [(0, 0)] * 3)
for _ in range(iters(2000)):
  L, R, m = rnd(1, 8), rnd(1, 8), rnd(0, 30)
  check_colouring(L, R, [(rnd(0, L - 1), rnd(0, R - 1)) for _ in range(m)])
for _ in range(iters(200)):  # k-regular: k perfect matchings
  n, k, edges = rnd(1, 30), rnd(1, 8), []
  for _ in range(k):
    p = list(range(n))
    rng.shuffle(p)
    edges += list(enumerate(p))
  check_colouring(n, n, edges)


def check_vizing(n, edges):
  """Proper, with colours in [0, Delta]."""
  col = lib.vizing(n, edges)
  check_eq(len(col), len(edges))
  deg = [0] * n
  for u, v in edges:
    deg[u] += 1
    deg[v] += 1
  D, seen = max(deg, default=0), set()
  for (u, v), c in zip(edges, col):
    check(0 <= c <= D)
    check((u, c) not in seen and (v, c) not in seen)
    seen.add((u, c))
    seen.add((v, c))


check_vizing(1, [])
check_vizing(3, [(0, 1), (1, 2), (2, 0)])
for n in range(1, 12):  # complete graphs: odd n needs Delta + 1
  check_vizing(n, [(a, b) for a in range(n) for b in range(a)])
for _ in range(iters(2000)):
  n, p = rnd(1, 12), rng.random()
  edges = [(a, b) for a in range(n) for b in range(a) if rng.random() < p]
  rng.shuffle(edges)
  edges = [(b, a) if rnd(0, 1) else (a, b) for a, b in edges]
  check_vizing(n, edges)

L = R = 2000
edges = [(rnd(0, L - 1), rnd(0, R - 1)) for _ in range(200_000)]
with Timer("2e5 edges on 2000 + 2000 vertices"):
  lib.edge_colouring(L, R, edges)
check_colouring(L, R, edges)
n = 2000
edges = list({(min(e), max(e)) for e in edges if e[0] != e[1]})
with Timer(f"vizing, {len(edges)} edges on 2000 vertices"):
  lib.vizing(n, edges)
check_vizing(n, edges)
print("OK graph/edge-colouring")
