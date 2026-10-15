from testlib import Timer, check_eq, iters, load, rnd, rng

lib = load("tree/lca")


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


def naive_path(par, dep, u, v):
  a = naive_lca(par, dep, u, v)
  left, right = [], []
  while u != a:
    left.append(u)
    u = par[u]
  while v != a:
    right.append(v)
    v = par[v]
  return left + [a] + right[::-1]


def is_anc(par, a, v):
  while v >= 0 and v != a:
    v = par[v]
  return v == a


def check_tree(adj, root):
  n = len(adj)
  lc = lib.LCA(adj, root)
  par, dep = bfs_parents(adj, root)
  check_eq(len(lc.order), n)
  for v in range(n):
    check_eq(lc.depth[v], dep[v])
    check_eq(lc.up[0][v], root if v == root else par[v])
    check_eq(lc.order[lc.tin[v]], v)
  for v in range(n):  # each subtree is a contiguous block of tin
    inside = [is_anc(par, v, x) for x in range(n)]
    size, t = sum(inside), lc.tin
    for x in range(n):
      check_eq(t[v] <= t[x] < t[v] + size, inside[x])
  for u in range(n):
    for k in range(dep[u] + 3):
      want = u
      for _ in range(k):
        want = par[want] if want >= 0 else -1
      check_eq(lc.kth_ancestor(u, k), want if k <= dep[u] else -1)
    for v in range(n):
      p = naive_path(par, dep, u, v)
      check_eq(lc.lca(u, v), naive_lca(par, dep, u, v))
      check_eq(lc.dist(u, v), len(p) - 1)
      for k in range(len(p) + 2):
        check_eq(lc.jump(u, v, k), p[k] if k < len(p) else -1)


def check_example(adj, root):
  """The header's Example: max edge weight on a path, with weights
  of both signs (a path of only negative edges must not give 0)."""
  n = len(adj)
  lc = lib.LCA(adj, root)
  par, dep = bfs_parents(adj, root)
  w = [rnd(-(10**9), 10**9) for _ in range(n)]
  mx = [w]
  for k in range(1, len(lc.up)):
    m, up = mx[-1], lc.up[k - 1]
    mx.append([max(m[v], m[up[v]]) for v in range(n)])

  def climb(v, d):
    r, k = -(1 << 62), 0
    while d:
      if d & 1:
        r, v = max(r, mx[k][v]), lc.up[k][v]
      d, k = d >> 1, k + 1
    return r

  def path_max(u, v):
    a = lc.depth[lc.lca(u, v)]
    du, dv = lc.depth[u] - a, lc.depth[v] - a
    return max(climb(u, du), climb(v, dv))

  for _ in range(3 * n):
    u, v = rnd(0, n - 1), rnd(0, n - 1)
    p = naive_path(par, dep, u, v)
    a = naive_lca(par, dep, u, v)
    if u != v:
      check_eq(path_max(u, v), max(w[x] for x in p if x != a))


for it in range(iters(1000)):
  n = rnd(1, 4 if it < 100 else 24)
  adj, root = random_tree(n, it % 4), rnd(0, n - 1)
  check_tree(adj, root)
  check_example(adj, root)

# Large: random recursive tree (shallow: naive walk is cheap) and a
# path rooted at one end (depth 2e5: the build must be iterative).
n = 200_000
adj = random_tree(n, 0)
root = rnd(0, n - 1)
qs = [(rnd(0, n - 1), rnd(0, n - 1)) for _ in range(n)]
with Timer("random tree n = 2e5: build + 2e5 lca"):
  lc = lib.LCA(adj, root)
  got = [lc.lca(u, v) for u, v in qs]
par, dep = bfs_parents(adj, root)
check_eq(got, [naive_lca(par, dep, u, v) for u, v in qs])

perm = list(range(n))
rng.shuffle(perm)
path = [[] for _ in range(n)]
for i in range(1, n):
  path[perm[i]].append(perm[i - 1])
  path[perm[i - 1]].append(perm[i])
qs = [(rnd(0, n - 1), rnd(0, n - 1), rnd(0, n)) for _ in range(n)]
with Timer("path n = 2e5: build + 2e5 lca, dist, jump"):
  lp = lib.LCA(path, perm[0])
  got = [lp.lca(perm[i], perm[j]) for i, j, _ in qs]
  dists = [lp.dist(perm[i], perm[j]) for i, j, _ in qs]
  jumps = [lp.jump(perm[i], perm[j], k) for i, j, k in qs]
check_eq(got, [perm[min(i, j)] for i, j, _ in qs])
check_eq(dists, [abs(i - j) for i, j, _ in qs])
want = []
for i, j, k in qs:
  s = i + k if i <= j else i - k
  want.append(perm[s] if k <= abs(i - j) else -1)
check_eq(jumps, want)
print("OK tree/lca")
