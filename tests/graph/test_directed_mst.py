from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("graph/directed-mst")


def brute(n, root, es):
  """Every parent choice; keep those where all parents lead to root."""
  ins = [[] for _ in range(n)]
  for a, b, w in es:
    if b != root and a != b:
      ins[b].append((a, w))
  best = None
  par = [-1] * n

  def ok():
    for s in range(n):
      x, steps = s, 0
      while x != root and steps <= n:
        x, steps = par[x], steps + 1
      if x != root:
        return False
    return True

  def rec(v, cost):
    nonlocal best
    if v == n:
      if ok() and (best is None or cost < best):
        best = cost
      return
    if v == root:
      rec(v + 1, cost)
      return
    for a, w in ins[v]:
      par[v] = a
      rec(v + 1, cost + w)

  rec(0, 0)
  return best


def slow(n, root, es):
  """Textbook O(n m) Chu-Liu/Edmonds with explicit contraction."""
  es = [list(e) for e in es]
  big, res = 1 << 62, 0
  while True:
    inw, pre = [big] * n, [-1] * n
    for a, b, w in es:
      if a != b and w < inw[b]:
        inw[b], pre[b] = w, a
    inw[root] = 0
    if big in inw:
      return None
    cid, vis, cnt = [-1] * n, [-1] * n, 0
    for v in range(n):
      res += inw[v]
      x = v
      while vis[x] != v and cid[x] < 0 and x != root:
        vis[x] = v
        x = pre[x]
      if x != root and cid[x] < 0:
        y = pre[x]
        while y != x:
          cid[y] = cnt
          y = pre[y]
        cid[x] = cnt
        cnt += 1
    if not cnt:
      return res
    for v in range(n):
      if cid[v] < 0:
        cid[v] = cnt
        cnt += 1
    for e in es:
      a, b = e[0], e[1]
      if cid[a] != cid[b]:
        e[2] -= inw[b]
      e[0], e[1] = cid[a], cid[b]
    n, root = cnt, cid[root]


def check_tree(n, root, es, cost, par):
  check_eq(len(par), n)
  check_eq(par[root], -1)
  cheapest = {}
  for a, b, w in es:
    if (a, b) not in cheapest or w < cheapest[a, b]:
      cheapest[a, b] = w
  total = 0
  for v in range(n):
    if v != root:
      check((par[v], v) in cheapest)
      total += cheapest[par[v], v]
  check_eq(total, cost)
  state = [0] * n  # 0 new, 1 on the current walk, 2 reaches root
  state[root] = 2
  for s in range(n):
    walk, x = [], s
    while state[x] == 0:
      state[x] = 1
      walk.append(x)
      x = par[x]
    check_eq(state[x], 2)
    for y in walk:
      state[y] = 2


def run(n, root, es, want):
  g = lib.DirectedMST(n)
  for a, b, w in es:
    g.add_edge(a, b, w)
  for _ in range(2):  # solve() twice: it must not change the graph
    cost, par = g.solve(root)
    if want is None:
      check_eq((cost, par), (-1, []))
    else:
      check_eq(cost, want)
      check_tree(n, root, es, cost, par)


def random_edges(n, m, lo, hi):
  return [(rnd(0, n - 1), rnd(0, n - 1), rnd(lo, hi)) for _ in range(m)]


def free_root(n, es):
  """Tweak: super vertex n with edges n -> v of weight 1 + sum |w|;
  a single root works iff exactly one par[v] == n."""
  big = 1 + sum(abs(w) for _, _, w in es)
  g = lib.DirectedMST(n + 1)
  for a, b, w in es:
    g.add_edge(a, b, w)
  for v in range(n):
    g.add_edge(n, v, big)
  cost, par = g.solve(n)
  costs = [brute(n, r, es) for r in range(n)]
  costs = [c for c in costs if c is not None]
  if par.count(n) == 1:
    check_eq(cost - big, min(costs))
  else:
    check_eq(costs, [])


run(1, 0, [], 0)
run(1, 0, [(0, 0, -5)], 0)
run(2, 0, [], None)
run(2, 1, [(0, 1, 3)], None)
run(2, 0, [(0, 1, -1)], -1)  # cost -1 with a valid answer
run(3, 0, [(1, 2, 1), (2, 1, 1)], None)  # unreachable cycle
for _ in range(iters(2000)):
  n = rnd(1, 7)
  es = random_edges(n, rnd(0, 3 * n), -6 if rnd(0, 1) else 0, 10)
  root = rnd(0, n - 1)
  want = brute(n, root, es)
  check_eq(slow(n, root, es), want)
  run(n, root, es, want)
for _ in range(iters(500)):
  n = rnd(1, 6)
  free_root(n, random_edges(n, rnd(0, 2 * n), -6 if rnd(0, 1) else 0, 10))
for _ in range(iters(100)):
  n = rnd(2, 60)
  es = random_edges(n, rnd(n, 6 * n), -1000, 1000)
  es += [(rnd(0, v - 1), v, rnd(0, 2000)) for v in range(1, n)]
  run(n, 0, es, slow(n, 0, es))

# One big cycle of cheap edges: the whole ring contracts.
n = 100_000
es = [(v, v % (n - 1) + 1, 1) for v in range(1, n)]
root_edges = [(0, v, rnd(1000, 10**6)) for v in range(1, n)]
with Timer("ring of 1e5 vertices"):
  run(n, 0, es + root_edges, min(e[2] for e in root_edges) + n - 2)
# Weights k, k+1, k+2, k-1, k+3, k-2, ..., then one more high one,
# into one vertex: a skew heap whose right spine has m/2 nodes, so
# a recursive merge would hit the recursion limit.
k = 100_000
es = [(0, 1, k), (0, 1, k + 1)]
for i in range(1, k):
  es += [(0, 1, k + 1 + i), (0, 1, k - i)]
es.append((0, 1, 2 * k + 1))
with Timer("skew heap with a spine of 1e5"):
  run(2, 0, es, 1)
n = 1000
es = random_edges(n, 20 * n, 0, 10**9)
es += [(v - 1, v, rnd(0, 10**9)) for v in range(1, n)]
run(n, 0, es, slow(n, 0, es))
n = 100_000
es = random_edges(n, 2 * n, 0, 10**9)
es += [(rnd(0, v - 1), v, rnd(0, 10**9)) for v in range(1, n)]
g = lib.DirectedMST(n)
for a, b, w in es:
  g.add_edge(a, b, w)
with Timer("n = 1e5, m = 3e5"):
  cost, par = g.solve(0)
check_tree(n, 0, es, cost, par)
print("OK graph/directed-mst")
