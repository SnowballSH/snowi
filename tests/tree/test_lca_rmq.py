from testlib import Timer, check_eq, iters, load, rnd, rng

lib = load("tree/lca-rmq")


def random_tree(n, shape):
  """Shapes: 0 random, 1 path, 2 star, 3 long thin; labels shuffled."""
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


def is_anc(par, a, v):
  while v >= 0 and v != a:
    v = par[v]
  return v == a


def check_tree(adj, root):
  n = len(adj)
  lc = lib.RmqLCA(adj, root)
  par, dep = bfs_parents(adj, root)
  check_eq(len(lc.order), n)
  check_eq(lc.par, par)
  check_eq(lc.depth, dep)
  t = lc.tin
  for v in range(n):
    check_eq(lc.order[t[v]], v)
    inside = [is_anc(par, v, x) for x in range(n)]
    size = sum(inside)
    for x in range(n):
      check_eq(t[v] <= t[x] < t[v] + size, inside[x])
  for u in range(n):
    for v in range(n):
      a = naive_lca(par, dep, u, v)
      check_eq(lc.lca(u, v), a)
      check_eq(lc.dist(u, v), dep[u] + dep[v] - 2 * dep[a])
  # The header's Example: Steiner tree size of a random set s.
  s = [v for v in range(n) if rnd(0, 2) == 0]
  s = sorted(s, key=lambda v: lc.tin[v])
  got = sum(lc.dist(s[i - 1], s[i]) for i in range(len(s))) // 2
  want = 0
  for v in range(n):
    if v != root:
      inside = sum(is_anc(par, v, x) for x in s)
      want += 0 < inside < len(s)
  check_eq(got, want)


for it in range(iters(1000)):
  n = rnd(1, 4 if it < 100 else 24)
  check_tree(random_tree(n, it % 4), rnd(0, n - 1))

# Large: random recursive tree (shallow: naive walk is cheap) and a
# path rooted at one end (depth 2e5: the build must be iterative).
n = 200_000
adj = random_tree(n, 0)
root = rnd(0, n - 1)
qs = [(rnd(0, n - 1), rnd(0, n - 1)) for _ in range(n)]
with Timer("random tree n = 2e5: build + 2e5 lca"):
  lc = lib.RmqLCA(adj, root)
  got = [lc.lca(u, v) for u, v in qs]
par, dep = bfs_parents(adj, root)
check_eq(got, [naive_lca(par, dep, u, v) for u, v in qs])

perm = list(range(n))
rng.shuffle(perm)
path = [[] for _ in range(n)]
for i in range(1, n):
  path[perm[i]].append(perm[i - 1])
  path[perm[i - 1]].append(perm[i])
qs = [(rnd(0, n - 1), rnd(0, n - 1)) for _ in range(n)]
with Timer("path n = 2e5: build + 2e5 lca and dist"):
  lp = lib.RmqLCA(path, perm[0])
  got = [lp.lca(perm[i], perm[j]) for i, j in qs]
  dists = [lp.dist(perm[i], perm[j]) for i, j in qs]
check_eq(got, [perm[min(i, j)] for i, j in qs])
check_eq(dists, [abs(i - j) for i, j in qs])
print("OK tree/lca-rmq")
