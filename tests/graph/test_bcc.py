from testlib import Timer, check_eq, iters, load, rnd, rng

lib = load("graph/bcc")


def labels(n, edges, skip_v=-1, skip_e=-1):
  """Connected components avoiding vertex skip_v and edge skip_e."""
  lab = [-1] * n
  for s in range(n):
    if s == skip_v or lab[s] >= 0:
      continue
    lab[s] = s
    stack = [s]
    while stack:
      u = stack.pop()
      for e, (a, b) in enumerate(edges):
        if e == skip_e or skip_v in (a, b) or u not in (a, b):
          continue
        v = a ^ b ^ u
        if lab[v] < 0:
          lab[v] = s
          stack.append(v)
  return lab


def brute_bridge(n, edges, e):
  u, v = edges[e]
  if u == v:
    return False
  lab = labels(n, edges, skip_e=e)
  return lab[u] != lab[v]


def brute_cut(n, edges, v):
  lab = labels(n, edges, skip_v=v)
  nb = [a ^ b ^ v for a, b in edges if a != b and v in (a, b)]
  return len({lab[a] for a in nb}) > 1


def cycle_classes(n, edges):
  """Union the edges of every simple cycle (length >= 2)."""
  m = len(edges)
  parent = list(range(m))

  def find(x):
    while parent[x] != x:
      x = parent[x]
    return x

  def dfs(s, u, path, on):
    for e, (a, b) in enumerate(edges):
      if a == b or u not in (a, b) or e in path:
        continue
      v = a ^ b ^ u
      if v == s:
        for f in path:
          parent[find(f)] = find(e)
      elif v not in on:
        dfs(s, v, path + [e], on | {v})

  for s in range(n):
    dfs(s, s, [], {s})
  return [find(e) for e in range(m)]


def check_tree(n, edges, b):
  check_eq(len(b.tree), n + b.nb)
  want = [set() for _ in range(n + b.nb)]
  for e, (u, v) in enumerate(edges):
    if b.comp[e] >= 0:
      for x in (u, v):
        want[x].add(n + b.comp[e])
        want[n + b.comp[e]].add(x)
  for x in range(n + b.nb):
    check_eq(sorted(b.tree[x]), sorted(want[x]))
  check_eq(b.cut, [len(b.tree[v]) >= 2 for v in range(n)])


def check_separators(n, edges, b):
  """Notes: the vertices whose removal separates connected u != v
  are exactly the vertex nodes strictly inside their tree path."""
  lab = labels(n, edges)
  without = [labels(n, edges, skip_v=x) for x in range(n)]
  for u in range(n):
    par, q = {u: -1}, [u]
    for y in q:
      for z in b.tree[y]:
        if z not in par:
          par[z] = y
          q.append(z)
    for v in range(n):
      if v == u or lab[u] != lab[v]:
        continue
      inner, y = set(), par[v]
      while y != u:
        if y < n:
          inner.add(y)
        y = par[y]
      sep = {
        x
        for x in range(n)
        if x not in (u, v) and without[x][u] != without[x][v]
      }
      check_eq(sep, inner, f"separators of {u}, {v}")


for _ in range(iters(2000)):
  n, m = rnd(1, 7), rnd(0, 9)
  edges = [(rnd(0, n - 1), rnd(0, n - 1)) for _ in range(m)]
  b = lib.BCC(n, edges)
  cls = cycle_classes(n, edges)
  check_eq(b.bridge, [brute_bridge(n, edges, e) for e in range(m)])
  check_eq(b.cut, [brute_cut(n, edges, v) for v in range(n)])
  real = [e for e in range(m) if edges[e][0] != edges[e][1]]
  for e in range(m):
    if e not in real:
      check_eq(b.comp[e], -1)
  check_eq(sorted({b.comp[e] for e in real}), list(range(b.nb)))
  for e in real:
    for f in real:
      check_eq(b.comp[e] == b.comp[f], cls[e] == cls[f])
  check_tree(n, edges, b)
  check_separators(n, edges, b)
# Path 0-1-2-3: only 2 separates 1 and 3, though 1 is a cut vertex.
p4 = [(0, 1), (1, 2), (2, 3)]
check_separators(4, p4, lib.BCC(4, p4))

# k cycles of length L chained by bridges, with an extra parallel
# edge in each cycle and a self-loop; labels and edges shuffled.
k, L = 20000, 10
n = k * L
lab = list(range(n))
rng.shuffle(lab)
edges, kind = [], []  # kind: -1 loop, -2 bridge, else cycle id
for c in range(k):
  base = c * L
  for i in range(L):
    edges.append((lab[base + i], lab[base + (i + 1) % L]))
    kind.append(c)
  edges.append((lab[base], lab[base + 1]))
  kind.append(c)
  edges.append((lab[base + 2], lab[base + 2]))
  kind.append(-1)
  if c + 1 < k:
    edges.append((lab[base + L - 1], lab[base + L]))
    kind.append(-2)
order = list(range(len(edges)))
rng.shuffle(order)
edges = [edges[e] for e in order]
kind = [kind[e] for e in order]
with Timer(f"n = {n}, m = {len(edges)}"):
  b = lib.BCC(n, edges)
check_eq(b.nb, 2 * k - 1)
check_eq(b.bridge, [x == -2 for x in kind])
block_of = {}
for e, x in enumerate(kind):
  if x == -1:
    check_eq(b.comp[e], -1)
  elif x >= 0:
    check_eq(block_of.setdefault(x, b.comp[e]), b.comp[e])
check_eq(len(set(block_of.values())), k)
check_eq(sum(b.cut), 2 * (k - 1))
check_tree(n, edges, b)
path = [(v, v + 1) for v in range(n - 1)]
with Timer(f"path of {n} vertices"):
  p = lib.BCC(n, path)
check_eq(p.nb, n - 1)
check_eq(p.cut, [0 < v < n - 1 for v in range(n)])
print("OK graph/bcc")
