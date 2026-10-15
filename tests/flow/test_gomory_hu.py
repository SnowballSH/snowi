from collections import deque

from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("flow/gomory-hu")
gomory_hu, Dinic = lib.gomory_hu, lib.Dinic
INF = float("inf")


def cut_weight(edges, inside):
  return sum(c for u, v, c in edges if inside[u] != inside[v])


def brute(n, edges):
  """All-pairs minimum cut by enumerating every vertex subset."""
  lam = [[INF] * n for _ in range(n)]
  for mask in range(1 << n):
    inside = [mask >> v & 1 for v in range(n)]
    c = cut_weight(edges, inside)
    for u in range(n):
      if inside[u]:
        for v in range(n):
          if not inside[v] and c < lam[u][v]:
            lam[u][v] = c
  return lam


def edmonds_karp(n, edges, s, t):
  c = [[0] * n for _ in range(n)]
  for u, v, w in edges:
    if u != v:
      c[u][v] += w
      c[v][u] += w
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
    f, v = INF, t
    while v != s:
      f = min(f, c[par[v]][v])
      v = par[v]
    v = t
    while v != s:
      c[par[v]][v] -= f
      c[v][par[v]] += f
      v = par[v]
    total += f


def certify(n, edges, tree):
  """tree is a spanning tree whose every edge splits the vertices
  into a cut of exactly its weight. Returns the lightest edge on
  each tree path (INF on the diagonal)."""
  check_eq(len(tree), max(n - 1, 0))
  adj = [[] for _ in range(n)]
  for i, (u, v, w) in enumerate(tree):
    check(0 <= u < n and 0 <= v < n and u != v and w >= 0)
    adj[u].append((v, i))
    adj[v].append((u, i))
  low = []
  for s in range(n):
    row = [-1] * n
    row[s] = INF
    st = [s]
    while st:
      u = st.pop()
      for v, i in adj[u]:
        if row[v] < 0:
          row[v] = min(row[u], tree[i][2])
          st.append(v)
    check(min(row) >= 0, "tree not connected")
    low.append(row)
  for i, (a, b, w) in enumerate(tree):
    inside = [False] * n
    inside[a] = True
    st = [a]
    while st:
      u = st.pop()
      for v, j in adj[u]:
        if j != i and not inside[v]:
          inside[v] = True
          st.append(v)
    check(not inside[b])
    check_eq(cut_weight(edges, inside), w)
  return low


def random_edges(n, m, maxc):
  return [(rnd(0, n - 1), rnd(0, n - 1), rnd(0, maxc)) for _ in range(m)]


def clustered(n, m):
  """Heavy random clusters joined by light edges: deep, uneven trees."""
  cl = [rnd(0, 4) for _ in range(n)]
  e = []
  for _ in range(m):
    u, v = rnd(0, n - 1), rnd(0, n - 1)
    e.append((u, v, rnd(5, 20) if cl[u] == cl[v] else rnd(0, 2)))
  return e


check_eq(gomory_hu(0, []), [])
check_eq(gomory_hu(1, [(0, 0, 5)]), [])
for _ in range(iters(1000)):
  n, m = rnd(1, 7), rnd(0, 12)
  maxc = rnd(1, 10) if rnd(0, 3) else 10**15
  e = random_edges(n, m, maxc)
  low, lam = certify(n, e, gomory_hu(n, e)), brute(n, e)
  for u in range(n):
    for v in range(n):
      if u != v:
        check_eq(low[u][v], lam[u][v])

for _ in range(iters(25)):
  n = rnd(10, 25)
  m = rnd(n, 4 * n)
  e = random_edges(n, m, rnd(1, 20)) if rnd(0, 1) else clustered(n, m)
  low = certify(n, e, gomory_hu(n, e))
  for _ in range(10):
    s, t = rnd(0, n - 1), rnd(0, n - 2)
    t += t >= s
    check_eq(low[s][t], edmonds_karp(n, e, s, t))

# Large: certify the cut property, then compare pairs with fresh
# max flows.
n = 500
e = random_edges(n, 5000, 10**9)
with Timer("n = 500, m = 5000"):
  tree = gomory_hu(n, e)
low = certify(n, e, tree)
for _ in range(20):
  s, t = rnd(0, n - 1), rnd(0, n - 2)
  t += t >= s
  d = Dinic(n)
  for u, v, c in e:
    d.add_edge(u, v, c, c)
  check_eq(low[s][t], d.max_flow(s, t))
print("OK flow/gomory-hu")
