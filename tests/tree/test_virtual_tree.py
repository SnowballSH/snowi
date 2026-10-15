from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("tree/virtual-tree", "tree/lca-rmq")


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


def bfs_parents(adj, root):
  n = len(adj)
  par, dep, seen, q = [-1] * n, [0] * n, [False] * n, [root]
  seen[root] = True
  for u in q:
    for v in adj[u]:
      if not seen[v]:
        seen[v] = True
        par[v], dep[v] = u, dep[u] + 1
        q.append(v)
  return par, dep


def naive_lca(par, dep, u, v):
  while dep[u] > dep[v]:
    u = par[u]
  while dep[v] > dep[u]:
    v = par[v]
  while u != v:
    u, v = par[u], par[v]
  return u


def check_vt(par, dep, tin, keys, vs, vpar):
  """Closure under pairwise LCA, sorted by tin; parent = nearest
  proper ancestor inside the set."""
  s = set(keys)
  for x in keys:
    for y in keys:
      s.add(naive_lca(par, dep, x, y))
  check_eq(sorted(vs), sorted(s))
  check_eq(len(vs), len(s))
  check(len(s) <= max(0, 2 * len(set(keys)) - 1))
  check_eq(len(vpar), len(vs))
  for i, v in enumerate(vs):
    if i:
      check(tin[vs[i - 1]] < tin[v])
    x = par[v]
    while x >= 0 and x not in s:
      x = par[x]
    check_eq(vs[vpar[i]] if vpar[i] >= 0 else -1, x)


def check_tree(adj, root):
  n = len(adj)
  lc = lib.LCA(adj, root)
  rl = lib.RmqLCA(adj, root)
  par, dep = bfs_parents(adj, root)
  for q in range(8):
    k = 0 if q == 0 else 1 if q == 1 else rnd(0, 2 * n)
    keys = [rnd(0, n - 1) for _ in range(k)]
    vs, vpar = lib.virtual_tree(lc, keys)
    check_vt(par, dep, lc.tin, keys, vs, vpar)
    vs2, vpar2 = lib.virtual_tree(rl, keys)
    check_vt(par, dep, rl.tin, keys, vs2, vpar2)
    # The header's Example: sum of distances over pairs of keys.
    ks = set(keys)
    cnt = [int(v in ks) for v in vs]
    total = 0
    for i in range(len(vs) - 1, 0, -1):
      d = lc.depth[vs[i]] - lc.depth[vs[vpar[i]]]
      total += d * cnt[i] * (len(ks) - cnt[i])
      cnt[vpar[i]] += cnt[i]
    want = 0
    for x in ks:
      for y in ks:
        if x < y:
          want += dep[x] + dep[y] - 2 * dep[naive_lca(par, dep, x, y)]
    check_eq(total, want)


for it in range(iters(800)):
  n = rnd(1, 4 if it < 100 else 30)
  check_tree(random_tree(n, it % 5), rnd(0, n - 1))

# Large: queries with sum k = 2e5 on a 2e5 path rooted at one end
# (nearest ancestor in the set = previous index) and a random tree.
n = 200_000
perm = list(range(n))
rng.shuffle(perm)
idx = [0] * n
for i, v in enumerate(perm):
  idx[v] = i
path = [[] for _ in range(n)]
for i in range(1, n):
  path[perm[i]].append(perm[i - 1])
  path[perm[i - 1]].append(perm[i])
lp = lib.LCA(path, perm[0])
queries, done = [], 0
while done < n:
  k = rnd(1, 2000)
  queries.append([rnd(0, n - 1) for _ in range(k)])
  done += k
with Timer("path n = 2e5: virtual trees, sum k = 2e5"):
  got = [lib.virtual_tree(lp, keys) for keys in queries]
for keys, (vs, vpar) in zip(queries, got):
  pos = sorted(set(idx[v] for v in keys))
  check_eq([idx[v] for v in vs], pos)
  check_eq(vpar, list(range(-1, len(pos) - 1)))

adj = random_tree(n, 0)
lc, rl = lib.LCA(adj, 0), lib.RmqLCA(adj, 0)
with Timer("random tree n = 2e5: virtual trees (LCA), sum k = 2e5"):
  got = [lib.virtual_tree(lc, keys) for keys in queries]
with Timer("random tree n = 2e5: virtual trees (RmqLCA), sum k = 2e5"):
  got2 = [lib.virtual_tree(rl, keys) for keys in queries]
check(got == got2)
for keys, (vs, vpar) in zip(queries, got):
  check(len(vs) <= 2 * len(keys) - 1)
  for i in range(1, len(vs)):
    p = vs[vpar[i]]
    check(lc.tin[p] < lc.tin[vs[i]] and lc.lca(p, vs[i]) == p)
print("OK tree/virtual-tree")
