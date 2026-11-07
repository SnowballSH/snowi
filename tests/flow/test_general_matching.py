from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("flow/general-matching")
general_matching = lib.general_matching


def adjacency(n, edges):
  g = [[] for _ in range(n)]
  for u, v in edges:
    g[u].append(v)
    if u != v:
      g[v].append(u)
  for a in g:
    rng.shuffle(a)
  return g


def brute(n, edges):
  """Maximum matching by DP over the set of vertices still free."""
  nb = [0] * n
  for u, v in edges:
    if u != v:
      nb[u] |= 1 << v
      nb[v] |= 1 << u
  f = [0] * (1 << n)
  for mask in range(1, 1 << n):
    i = (mask & -mask).bit_length() - 1
    rest = mask & ~(1 << i)
    best, m = f[rest], nb[i] & rest
    while m:
      best = max(best, 1 + f[rest & ~(m & -m)])
      m &= m - 1
    f[mask] = best
  return f[-1]


def tutte(n, edges):
  """Rank of a random Tutte matrix mod p: twice the maximum
  matching size, except with probability below n / p."""
  p = 1_000_000_007
  a = [[0] * n for _ in range(n)]
  for u, v in edges:
    if u != v:
      x = rnd(1, p - 1)
      a[u][v] = (a[u][v] + x) % p
      a[v][u] = (a[v][u] - x) % p
  rank = 0
  for c in range(n):
    r = next((r for r in range(rank, n) if a[r][c]), n)
    if r == n:
      continue
    a[r], a[rank] = a[rank], a[r]
    inv = pow(a[rank][c], p - 2, p)
    for i in range(rank + 1, n):
      f = a[i][c] * inv % p
      if f:
        a[i] = [(x - f * y) % p for x, y in zip(a[i], a[rank])]
    rank += 1
  return rank


def certify(n, edges, k, mate):
  """mate is a matching of k pairs that uses edges of the graph."""
  check_eq(len(mate), n)
  es = set(edges) | {(v, u) for u, v in edges}
  ends = 0
  for v in range(n):
    if mate[v] >= 0:
      check(mate[v] < n and mate[v] != v)
      check_eq(mate[mate[v]], v)
      check((v, mate[v]) in es)
      ends += 1
    else:
      check_eq(mate[v], -1)
  check_eq(ends, 2 * k)


def solve(n, edges):
  k, mate = general_matching(adjacency(n, edges))
  certify(n, edges, k, mate)
  return k


def random_edges(n, m):
  if not n:
    return []
  return [(rnd(0, n - 1), rnd(0, n - 1)) for _ in range(m)]


def flowers(n):
  """Odd cycles a, v, ..., end, a through an earlier vertex a:
  nested blossoms."""
  e, v = [], 1
  while v < n:
    a, end = rnd(0, v - 1), min(n, v + 2 * rnd(1, 3)) - 1
    e.append((a, v))
    e += [(x, x + 1) for x in range(v, end)]
    e.append((end, a))
    v = end + 1
  return e


check_eq(solve(0, []), 0)
check_eq(solve(1, [(0, 0)]), 0)
check_eq(solve(3, [(0, 1), (1, 2), (2, 0)]), 1)
# Two triangles joined by an edge: the blossom must be entered.
two = [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 5), (5, 3)]
check_eq(solve(6, two), 3)
# Both sides of a blossom must be shrunk: shrinking only the scanned
# vertex's side finds 2 pairs with these exact adjacency lists.
g = [[5, 3], [2, 5], [1, 4, 3], [2, 5, 0], [2], [1, 0, 3]]
check_eq(general_matching(g)[0], 3)

for _ in range(iters(2500)):
  n = rnd(0, 11)
  m = rnd(0, n * (n + 1) // 2 + 2)
  e = random_edges(n, m) if rnd(0, 4) else flowers(n)
  check_eq(solve(n, e), brute(n, e))

for _ in range(iters(150)):
  n = rnd(10, 50)
  if rnd(0, 1):
    e = random_edges(n, rnd(n // 2, 3 * n))
  else:
    e = flowers(n)
  check_eq(2 * solve(n, e), tutte(n, e))

for n in range(30):  # complete graphs
  e = [(i, j) for i in range(n) for j in range(i)]
  check_eq(solve(n, e), n // 2)

# Large: a planted perfect matching hidden among random edges and
# odd cycles; disjoint triangles; a sparse random graph; a dense
# one.
n = 100_000
p = list(range(n))
rng.shuffle(p)
e = random_edges(n, 100_000) + flowers(n)
e += [(p[i], p[i + 1]) for i in range(0, n, 2)]
g = adjacency(n, e)
with Timer("planted perfect matching, n = 1e5, m = 2.5e5"):
  k, mate = general_matching(g)
check_eq(k, n // 2)
certify(n, e, k, mate)

n = 60_000
e = []
for i in range(0, n, 3):
  e += [(i, i + 1), (i + 1, i + 2), (i + 2, i)]
g = adjacency(n, e)
with Timer("2e4 disjoint triangles"):
  k, mate = general_matching(g)
check_eq(k, n // 3)
certify(n, e, k, mate)

n = 100_000
e = random_edges(n, 200_000)
g = adjacency(n, e)
with Timer("random graph, n = 1e5, m = 2e5"):
  k, mate = general_matching(g)
certify(n, e, k, mate)

n = 1000
e = [(i, j) for i in range(n) for j in range(i) if rnd(0, 4) == 0]
g = adjacency(n, e)
with Timer("dense random graph, n = 1000, m = 1e5"):
  k, mate = general_matching(g)
check_eq(k, n // 2)
certify(n, e, k, mate)

# Broom: a path 0..L-1, then n - L leaves on L - 1. Each leaf's
# search fails after walking the path, unless the trees of failed
# searches are retired.
n = 100_000
L = n // 2
e = [(i, i + 1) for i in range(L - 1)] + [(L - 1, i) for i in range(L, n)]
g = adjacency(n, e)
with Timer("broom, n = 1e5") as timer:
  k, mate = general_matching(g)
check(timer.seconds < 3, "failed searches are repeated")  # else ~3 min
check_eq(k, (L + 1) // 2)
certify(n, e, k, mate)
print("OK flow/general-matching")
