from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("graph/max-clique")


def brute(adj):
  """Largest clique by DP over all subsets: a set is a clique iff its
  lowest vertex sees the rest and the rest is a clique."""
  n, best = len(adj), 0
  ok = [False] * (1 << n)
  ok[0] = True
  for mask in range(1, 1 << n):
    low = (mask & -mask).bit_length() - 1
    rest = mask & (mask - 1)
    ok[mask] = ok[rest] and adj[low] & rest == rest
    if ok[mask]:
      best = max(best, bin(mask).count("1"))
  return best


def bk(adj, p, x, size):
  """Bron-Kerbosch with a pivot: an independent check for larger n."""
  if not p:
    return size if not x else 0
  px = p | x
  piv = max(
    (u for u in range(len(adj)) if px >> u & 1),
    key=lambda u: bin(p & adj[u]).count("1"),
  )
  best, cand = 0, p & ~adj[piv]
  while cand:
    low = cand & -cand
    v = low.bit_length() - 1
    best = max(best, bk(adj, p & adj[v], x & adj[v], size + 1))
    p &= ~low
    x |= low
    cand ^= low
  return best


def random_graph(n, prob):
  adj = [0] * n
  for u in range(n):
    for v in range(u):
      if rng.random() < prob:
        adj[u] |= 1 << v
        adj[v] |= 1 << u
  return adj


def check_clique(adj, c, want):
  check_eq(len(c), want)
  check_eq(c, sorted(set(c)))
  for i, u in enumerate(c):
    check(0 <= u < len(adj))
    for v in c[:i]:
      check(adj[u] >> v & 1)


check_eq(lib.max_clique([]), [])
check_eq(lib.max_clique([0]), [0])
check_eq(len(lib.max_clique([0] * 5)), 1)
check_eq(lib.max_clique([127 & ~(1 << u) for u in range(7)]), list(range(7)))
for _ in range(iters(2000)):
  adj = random_graph(rnd(1, 11), rng.random())
  check_clique(adj, lib.max_clique(adj), brute(adj))
for _ in range(iters(30)):
  n = rnd(15, 32)
  adj = random_graph(n, rng.uniform(0.3, 0.85))
  check_clique(adj, lib.max_clique(adj), bk(adj, (1 << n) - 1, 0, 0))

# A planted clique of 24 in G(128, 1/2), whose own cliques have
# about 11 vertices.
n = 128
adj = random_graph(n, 0.5)
planted = rng.sample(range(n), 24)
for u in planted:
  for v in planted:
    if u != v:
      adj[u] |= 1 << v
with Timer("planted clique, n = 128"):
  check_eq(lib.max_clique(adj), sorted(planted))
adj = random_graph(100, 0.9)
with Timer("n = 100, p = 0.9"):
  c = lib.max_clique(adj)
check_clique(adj, c, len(c))
print("OK graph/max-clique")
