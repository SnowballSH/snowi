from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("graph/euler-path")


def brute(u, left, edges, directed, used):
  if not left:
    return True
  for e, (a, b) in enumerate(edges):
    if used[e]:
      continue
    for x, y in [(a, b)] if directed else [(a, b), (b, a)]:
      if x == u:
        used[e] = True
        ok = brute(y, left - 1, edges, directed, used)
        used[e] = False
        if ok:
          return True
  return False


def check_walk(edges, directed, src, res):
  check(res is not None)
  vs, es = res
  m = len(edges)
  check_eq((len(vs), len(es), vs[0]), (m + 1, m, src))
  check_eq(sorted(es), list(range(m)))
  for i, e in enumerate(es):
    a, b = edges[e]
    fwd = (a, b) == (vs[i], vs[i + 1])
    check(fwd or (not directed and (b, a) == (vs[i], vs[i + 1])))


def walk_edges(n, m, src, directed):
  """Edges of a random walk from src, shuffled: a path exists."""
  edges, u = [], src
  for _ in range(m):
    v = rnd(0, n - 1)
    edges.append((v, u) if not directed and rnd(0, 1) else (u, v))
    u = v
  rng.shuffle(edges)
  return edges


def check_orientation(n, edges):
  """Notes: orient an undirected graph so that |in - out| <= 1, with
  one circuit per component and dummy edges pairing odd vertices."""
  adj = [[] for _ in range(n)]
  for u, v in edges:
    adj[u].append(v)
    adj[v].append(u)
  comp, bal = [-1] * n, [0] * n
  for s in range(n):
    if comp[s] >= 0 or not adj[s]:
      continue
    comp[s], q = s, [s]
    for x in q:
      for y in adj[x]:
        if comp[y] < 0:
          comp[y] = s
          q.append(y)
    ids = [e for e, (u, v) in enumerate(edges) if comp[u] == s]
    odd = [x for x in q if len(adj[x]) % 2]
    es = [edges[e] for e in ids] + list(zip(odd[::2], odd[1::2]))
    res = lib.euler_path(n, es, False, s)
    check(res is not None and res[0][-1] == s)
    vs, ep = res
    for i, e in enumerate(ep):
      if e < len(ids):
        bal[vs[i]] += 1
        bal[vs[i + 1]] -= 1
  check(all(abs(b) <= 1 for b in bal), f"{edges} {bal}")


found, total = 0, iters(3000)
for _ in range(total):
  n, m = rnd(1, 5), rnd(0, 7)
  src, directed = rnd(0, n - 1), bool(rnd(0, 1))
  if rnd(0, 1):
    edges = walk_edges(n, m, rnd(0, n - 1), directed)
  else:
    edges = [(rnd(0, n - 1), rnd(0, n - 1)) for _ in range(m)]
  want = brute(src, m, edges, directed, [False] * m)
  res = lib.euler_path(n, edges, directed, src)
  check_eq(res is not None, want)
  if want:
    found += 1
    check_walk(edges, directed, src, res)
check(total // 5 < found < total * 4 // 5, f"{found} of {total}")
for _ in range(iters(1000)):
  n, m = rnd(1, 8), rnd(0, 10)
  check_orientation(n, [(rnd(0, n - 1), rnd(0, n - 1)) for _ in range(m)])
# Two triangles: no single walk, but orientable per component.
tri2 = [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)]
check(lib.euler_path(6, tri2, False, 0) is None)
check_orientation(6, tri2)

n, m = 100_000, 200_000
for directed in (False, True):
  edges = walk_edges(n, m, 7, directed)
  with Timer(f"random walk, m = 2e5, directed = {directed}"):
    res = lib.euler_path(n, edges, directed, 7)
  check_walk(edges, directed, 7, res)
  edges += [(n - 1, n - 2), (n - 3, n - 4)]
  check(lib.euler_path(n, edges, directed, 7) is None)
cyc = [(v, (v + 1) % m) for v in range(m)]
check_walk(cyc, True, 0, lib.euler_path(m, cyc, True, 0))
print("OK graph/euler-path")
