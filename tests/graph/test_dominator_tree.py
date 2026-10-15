from testlib import Timer, check_eq, iters, load, rnd, rng

lib = load("graph/dominator-tree")


def reach(adj, s, banned):
  """Vertices reachable from s without passing through `banned`."""
  seen = [False] * len(adj)
  if s == banned:
    return seen
  seen[s] = True
  st = [s]
  while st:
    for v in adj[st.pop()]:
      if v != banned and not seen[v]:
        seen[v] = True
        st.append(v)
  return seen


def brute(adj, s):
  """d strictly dominates v iff removing d cuts v off from s. The
  strict dominators form a chain; idom is the one with the most
  dominators itself."""
  n = len(adj)
  ok = reach(adj, s, -1)
  dom = [[] for _ in range(n)]
  for d in range(n):
    if ok[d]:
      r = reach(adj, s, d)
      for v in range(n):
        if ok[v] and v != d and not r[v]:
          dom[v].append(d)
  idom = [-1] * n
  for v in range(n):
    if ok[v]:
      idom[v] = s if v == s else max(dom[v], key=lambda d: len(dom[d]))
  return idom


def dag_idom(adj, s):
  """The Notes' DAG rule: in topological order, idom[v] is the LCA of
  the predecessors of v that are reachable from s."""
  n = len(adj)
  ok = reach(adj, s, -1)
  preds, indeg = [[] for _ in range(n)], [0] * n
  for u in range(n):
    for v in adj[u]:
      preds[v].append(u)
      indeg[v] += 1
  order = [v for v in range(n) if indeg[v] == 0]
  for u in order:
    for v in adj[u]:
      indeg[v] -= 1
      if indeg[v] == 0:
        order.append(v)
  idom, depth = [-1] * n, [0] * n
  idom[s] = s
  for v in order:
    if v == s or not ok[v]:
      continue
    lca = -1
    for a in preds[v]:
      if not ok[a]:
        continue
      b = lca if lca >= 0 else a
      while a != b:
        if depth[a] < depth[b]:
          a, b = b, a
        a = idom[a]
      lca = a
    idom[v], depth[v] = lca, depth[lca] + 1
  return idom


def random_graph(n, m):
  adj = [[] for _ in range(n)]
  for _ in range(m):
    adj[rnd(0, n - 1)].append(rnd(0, n - 1))
  return adj


check_eq(lib.dominator_tree([[]], 0), [0])
check_eq(lib.dominator_tree([[0, 0]], 0), [0])
check_eq(lib.dominator_tree([[1], [], [1]], 0), [0, 0, -1])
check_eq(lib.dominator_tree([[1, 2], [3], [3], []], 0), [0, 0, 0, 0])
for _ in range(iters(3000)):
  n = rnd(1, 10)
  adj = random_graph(n, rnd(0, 3 * n))
  s = rnd(0, n - 1)
  check_eq(lib.dominator_tree(adj, s), brute(adj, s))
# DAGs (edges go up in number); unreachable vertices may point into
# reachable ones, e.g. 0 -> 1 <- 2 from 0.
check_eq(dag_idom([[1], [], [1]], 0), [0, 0, -1])
for _ in range(iters(1000)):
  n = rnd(1, 12)
  adj = [[] for _ in range(n)]
  for _ in range(rnd(0, 3 * n)):
    u, v = sorted(rng.sample(range(n), 2)) if n > 1 else (0, 0)
    if u < v:
      adj[u].append(v)
  s = rnd(0, n - 1)
  check_eq(lib.dominator_tree(adj, s), dag_idom(adj, s))
for _ in range(iters(50)):
  n = rnd(10, 200)
  adj = random_graph(n, rnd(n, 3 * n))
  check_eq(lib.dominator_tree(adj, 0), brute(adj, 0))

# A path plus back edges to ancestors: idom is the tree parent.
# Deep DFS and long compression paths.
n = 200_000
parent = [-1] * n
g = [[] for _ in range(n)]
for v in range(1, n):
  parent[v] = v - 1 if v < n // 2 else rnd(0, v - 1)
  g[parent[v]].append(v)
for v in range(1, n):
  for _ in range(2):
    a = v
    for _ in range(rnd(1, 4)):
      a = parent[a] if a else 0
    g[v].append(a)
g[n // 2 - 1].append(1)
for out in g:
  rng.shuffle(out)
with Timer("n = 2e5, m = 6e5"):
  idom = lib.dominator_tree(g, 0)
parent[0] = 0
check_eq(idom, parent)
g = random_graph(n, 3 * n)
with Timer("random n = 2e5, m = 6e5"):
  lib.dominator_tree(g, 0)
# Star: n DFS children of the root share one bucket.
star = [list(range(1, n))] + [[] for _ in range(n - 1)]
check_eq(lib.dominator_tree(star, 0), [0] * n)
# Path with edges from its end to every vertex: each step evaluates
# the whole path unless it is compressed.
path = [[v + 1] for v in range(n - 1)] + [list(range(1, n))]
with Timer("path with 2e5 back edges"):
  check_eq(lib.dominator_tree(path, 0), [0] + list(range(n - 1)))
print("OK graph/dominator-tree")
