from bisect import bisect_left

from testlib import CONTENT, Timer, check, check_eq, iters, load, rnd, rng

lib = load("tree/centroid")
INF = 1 << 30


def weighted_centroid():
  """The Tweak for weighted edges, applied to the snippet's source:
  a parallel `wt`, relabelled as `gw`, added per edge in the BFS."""
  src = (CONTENT / "tree/centroid.py").read_text()
  edits = [
    ("def __init__(self, adj):", "def __init__(self, adj, wt):"),
    (
      "    g = [[new[v] for v in adj[u]] for u in order]\n",
      "    g = [[new[v] for v in adj[u]] for u in order]\n"
      "    gw = [wt[u] for u in order]\n",
    ),
    ("for v in g[u]:", "for v, w in zip(g[u], gw[u]):"),
    ("bpar[v], d[v] = u, du\n", "bpar[v], d[v] = u, d[u] + w\n"),
  ]
  for old, new in edits:
    check_eq(src.count(old), 1, f"Tweak target {old!r}")
    src = src.replace(old, new)
  namespace = {}
  exec(compile(src, "centroid-weighted", "exec"), namespace)
  return namespace["Centroid"]


WeightedCentroid = weighted_centroid()


def check_weighted(adj):
  """Weighted distances to every centroid ancestor; same tree."""
  n = len(adj)
  wt = [[0] * len(a) for a in adj]
  for u in range(n):
    for i, v in enumerate(adj[u]):
      if u < v:
        w = rnd(0, 10**12)
        wt[u][i] = w
        wt[v][adj[v].index(u)] = w
  cd, plain = WeightedCentroid(adj, wt), lib.Centroid(adj)
  check_eq((cd.par, cd.level), (plain.par, plain.level))
  for c in range(n):
    d, q = {c: 0}, [c]  # weighted distances from c
    for u in q:
      for v, w in zip(adj[u], wt[u]):
        if v not in d:
          d[v] = d[u] + w
          q.append(v)
    for v in range(n):
      x = v
      while x >= 0 and x != c:
        x = cd.par[x]
      if x == c:
        check_eq(cd.dist[cd.level[c]][v], d[v])


def random_tree(n, shape):
  """Shapes: 0 random, 1 path, 2 star, 3 long thin, 4 binary-ish."""
  perm = list(range(n))
  rng.shuffle(perm)
  adj = [[] for _ in range(n)]
  for i in range(1, n):
    if shape == 1:
      p = i - 1
    elif shape == 2:
      p = 0
    elif shape == 3:
      p = rnd(max(0, i - 3), i - 1)
    elif shape == 4:
      p = (i - 1) // 2
    else:
      p = rnd(0, i - 1)
    adj[perm[i]].append(perm[p])
    adj[perm[p]].append(perm[i])
  for a in adj:
    rng.shuffle(a)
  return adj


def bfs_from(adj, s, ok=None):
  """Distances from s, within vertices with ok[v] (all if None)."""
  d = {s: 0}
  q = [s]
  for u in q:
    for v in adj[u]:
      if v not in d and (ok is None or ok[v]):
        d[v] = d[u] + 1
        q.append(v)
  return d


def check_tree(adj):
  n = len(adj)
  cd = lib.Centroid(adj)
  check_eq(len(cd.par), n)
  d = [bfs_from(adj, s) for s in range(n)]
  piece = [[] for _ in range(n)]
  check_eq(cd.par.count(-1), 1)
  for v in range(n):
    p = cd.par[v]
    check_eq(cd.level[v], 0 if p < 0 else cd.level[p] + 1)
    check(cd.level[v] <= n.bit_length() - 1)
    c = v
    while c >= 0:
      piece[c].append(v)
      check_eq(cd.dist[cd.level[c]][v], d[c][v])
      c = cd.par[c]
  for c in range(n):
    ok = [cd.level[v] >= cd.level[c] for v in range(n)]
    check_eq(sorted(bfs_from(adj, c, ok)), sorted(piece[c]))
    ok[c] = False
    for v in adj[c]:
      if ok[v]:
        check(2 * len(bfs_from(adj, v, ok)) <= len(piece[c]))
  # Usage: nearest marked vertex.
  best = [INF] * n

  def mark(v):
    c = v
    while c >= 0:
      best[c] = min(best[c], cd.dist[cd.level[c]][v])
      c = cd.par[c]

  def nearest(v):
    c, r = v, INF
    while c >= 0:
      r = min(r, best[c] + cd.dist[cd.level[c]][v])
      c = cd.par[c]
    return r

  marked = []
  for _ in range(2 * n):
    v = rnd(0, n - 1)
    if rnd(0, 2) == 0:
      mark(v)
      marked.append(v)
    check_eq(nearest(v), min([d[v][x] for x in marked], default=INF))
  # Tweak: sum of distances to a changing marked set.
  cnt, tot, sub, on = [0] * n, [0] * n, [0] * n, [False] * n
  for _ in range(2 * n):
    v = rnd(0, n - 1)
    if rnd(0, 1) and not on[v]:
      on[v] = True
      c = v
      while c >= 0:
        cnt[c] += 1
        tot[c] += cd.dist[cd.level[c]][v]
        if cd.par[c] >= 0:
          sub[c] += cd.dist[cd.level[c] - 1][v]
        c = cd.par[c]
    got, c, ch = 0, v, -1
    while c >= 0:
      dv = cd.dist[cd.level[c]][v]
      got += tot[c] + cnt[c] * dv
      if ch >= 0:
        got -= sub[ch] + cnt[ch] * dv
      c, ch = cd.par[c], c
    check_eq(got, sum(d[v][x] for x in range(n) if on[x]))
  # Example: pairs at distance exactly K, for every K.
  size = [0] * n
  for v in range(n):
    c = v
    while c >= 0:
      size[c], c = size[c] + 1, cd.par[c]
  for K in range(1, n + 1):
    cnt = [[0] * size[c] for c in range(n)]
    sub = [[0] * (size[c] + 1) for c in range(n)]
    pairs = 0
    for v in range(n):
      c, ch = v, -1
      while c >= 0:
        dv = cd.dist[cd.level[c]][v]
        e = K - dv
        if 0 <= e < size[c]:
          pairs += cnt[c][e]
        if ch >= 0 and 0 <= e <= size[ch]:
          pairs -= sub[ch][e]
        cnt[c][dv] += 1
        if ch >= 0:
          sub[ch][dv] += 1
        c, ch = cd.par[c], c
    want = sum(d[u][v] == K for u in range(n) for v in range(u))
    check_eq(pairs, want)


for it in range(iters(500)):
  n = rnd(1, 4 if it < 100 else 18)
  check_tree(random_tree(n, it % 5))
  check_weighted(random_tree(n, it % 5))

# Large: n = 2e5 path (deep), nearest marked vertex against sorted
# marked positions; random tree and star: build time, level bound.
n = 200_000
perm = list(range(n))
rng.shuffle(perm)
path = [[] for _ in range(n)]
for i in range(1, n):
  path[perm[i]].append(perm[i - 1])
  path[perm[i - 1]].append(perm[i])
with Timer("path n = 2e5: build"):
  cd = lib.Centroid(path)
ops = [rnd(0, n - 1) for _ in range(n)]
best, got = [INF] * n, []
par, level, dist = cd.par, cd.level, cd.dist
with Timer("path: 1e5 marks + 1e5 nearest"):
  for t, i in enumerate(ops):
    v = c = perm[i]
    if t % 2:
      while c >= 0:
        best[c] = min(best[c], dist[level[c]][v])
        c = par[c]
    else:
      r = INF
      while c >= 0:
        r = min(r, best[c] + dist[level[c]][v])
        c = par[c]
      got.append(r)
marked, want = [], []
for t, i in enumerate(ops):
  if t % 2:
    marked.insert(bisect_left(marked, i), i)
  else:
    j = bisect_left(marked, i)
    w = INF
    if j < len(marked):
      w = marked[j] - i
    if j > 0:
      w = min(w, i - marked[j - 1])
    want.append(w)
check_eq(got, want)
for shape in (0, 2, 4):
  adj = random_tree(n, shape)
  with Timer(f"shape {shape}, n = 2e5: build"):
    cd = lib.Centroid(adj)
  size = [0] * n
  for v in range(n):
    c = v
    while c >= 0:
      size[c] += 1
      c = cd.par[c]
  check(all(size[c] <= n >> cd.level[c] for c in range(n)))
print("OK tree/centroid")
