from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("graph/topo-sort")


def brute_has_cycle(adj):
  """A cycle exists iff some vertex reaches itself in >= 1 step."""
  n = len(adj)
  for s in range(n):
    seen, stack = set(), list(adj[s])
    while stack:
      v = stack.pop()
      if v == s:
        return True
      if v not in seen:
        seen.add(v)
        stack.extend(adj[v])
  return False


def check_order(adj, order):
  n = len(adj)
  check(order is not None)
  check_eq(sorted(order), list(range(n)))
  pos = [0] * n
  for i, v in enumerate(order):
    pos[v] = i
  for u in range(n):
    for v in adj[u]:
      check(pos[u] < pos[v])


def random_graph(n, m, dag):
  """Random graph; if dag, edges follow a random permutation."""
  rank = list(range(n))
  rng.shuffle(rank)
  adj = [[] for _ in range(n)]
  for _ in range(m):
    u, v = rnd(0, n - 1), rnd(0, n - 1)
    if dag and rank[u] >= rank[v]:
      if u == v:
        continue
      u, v = v, u
    adj[u].append(v)
  return adj


cyclic = 0
for _ in range(iters(4000)):
  n = rnd(0, 8)
  adj = random_graph(n, rnd(0, 16) if n else 0, rnd(0, 1))
  order = lib.topo_sort(adj)
  if brute_has_cycle(adj):
    check(order is None)
    cyclic += 1
  else:
    check_order(adj, order)
total = iters(4000)
check(total // 8 < cyclic < total * 7 // 8, f"{cyclic} cyclic")

n = 200_000
adj = random_graph(n, 400_000, True)
with Timer("2e5 vertices, 4e5 edges"):
  order = lib.topo_sort(adj)
check_order(adj, order)
chain = [[v + 1] for v in range(n - 1)] + [[]]
check_eq(lib.topo_sort(chain), list(range(n)))
chain[-1].append(0)
check(lib.topo_sort(chain) is None)
print("OK graph/topo-sort")
