from operator import add

from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("tree/hld", "ds/segtree", "ds/lazy-segtree")
MOD = 998244353


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


class Brute:
  def __init__(self, adj, root):
    n = len(adj)
    self.par, self.dep, self.size = [-1] * n, [0] * n, [1] * n
    seen, q = [False] * n, [root]
    seen[root] = True
    for u in q:
      for v in adj[u]:
        if not seen[v]:
          seen[v] = True
          self.par[v], self.dep[v] = u, self.dep[u] + 1
          q.append(v)
    for v in reversed(q[1:]):
      self.size[self.par[v]] += self.size[v]

  def lca(self, u, v):
    par, dep = self.par, self.dep
    while dep[u] > dep[v]:
      u = par[u]
    while dep[v] > dep[u]:
      v = par[v]
    while u != v:
      u, v = par[u], par[v]
    return u

  def is_anc(self, a, v):
    while v >= 0 and v != a:
      v = self.par[v]
    return v == a

  def path(self, u, v):
    """Vertices u .. v and sides: -1 u's side, 0 the lca, 1 v's."""
    a = self.lca(u, v)
    left, right = [], []
    while u != a:
      left.append(u)
      u = self.par[u]
    while v != a:
      right.append(v)
      v = self.par[v]
    side = [-1] * len(left) + [0] + [1] * len(right)
    return left + [a] + right[::-1], side


def then(f, g):  # map y -> b y + c: f, then g
  return g[0] * f[0] % MOD, (g[0] * f[1] + g[1]) % MOD


def check_tree(adj, root):
  n = len(adj)
  h = lib.HLD(adj, root)
  b = Brute(adj, root)
  check_eq(h.par, b.par)
  check_eq(h.depth, b.dep)
  check_eq(h.size, b.size)
  check_eq(sorted(h.pos), list(range(n)))
  at = [0] * n
  for v in range(n):
    at[h.pos[v]] = v
  for v in range(n):
    kids = [c for c in adj[v] if c != b.par[v]]
    if not kids:
      check_eq(h.heavy[v], -1)
    else:
      check_eq(b.par[h.heavy[v]], v)
      check_eq(b.size[h.heavy[v]], max(b.size[c] for c in kids))
      check_eq(h.pos[h.heavy[v]], h.pos[v] + 1)
    top = v == root or h.heavy[b.par[v]] != v
    check_eq(h.head[v], v if top else h.head[b.par[v]])
    l, r = h.subtree(v)
    for x in range(n):
      check_eq(l <= h.pos[x] < r, b.is_anc(v, x))
  maps = [(rnd(0, MOD - 1), rnd(0, MOD - 1)) for _ in range(n)]
  fwd = lib.SegTree([maps[at[i]] for i in range(n)], then, (1, 0))
  bwd = lib.SegTree(
    fwd.t[fwd.size : fwd.size + n], lambda f, g: then(g, f), (1, 0)
  )
  for u in range(n):
    for v in range(n):
      check_eq(h.lca(u, v), b.lca(u, v))
      p, side = b.path(u, v)
      where = {x: i for i, x in enumerate(p)}
      for edges in (False, True):
        want = [x for x, s in zip(p, side) if not edges or s]
        got = []
        for l, r in h.path(u, v, edges):
          check(l < r)
          s = min(side[where[at[x]]] for x in range(l, r))
          if s < 0:  # u's side: walked right to left
            got += [at[x] for x in range(r - 1, l - 1, -1)]
          else:
            got += [at[x] for x in range(l, r)]
        check_eq(got, want)
      # Notes: at most 2 log2 n light edges (log2 n each side), so
      # path() returns at most that many ranges plus 2.
      a = b.lca(u, v)
      light = sum(x != a and h.head[x] == x for x in p)
      check(light <= 2 * (n.bit_length() - 1))
      check(len(h.path(u, v)) <= light + 2)
      # Tweak: non-commutative fold with a swapped-op segment tree.
      w = h.lca(u, v)
      f = (1, 0)
      for l, r in h.path(u, w):
        f = then(f, bwd.query(l, r))
      for l, r in h.path(w, v, True):
        f = then(f, fwd.query(l, r))
      g = (1, 0)
      for x in p:
        g = then(g, maps[x])
      check_eq(f, g)


def check_usage(adj, root, ops):
  """The header's Usage: path and subtree queries on positions."""
  n = len(adj)
  h = lib.HLD(adj, root)
  br = Brute(adj, root)
  a = [rnd(-(10**9), 10**9) for _ in range(n)]
  b = [0] * n
  for v in range(n):
    b[h.pos[v]] = a[v]
  st = lib.SegTree(b, add, 0)
  lz = lib.LazySegTree(b, min, 10**30, add, add, 0)
  for _ in range(ops):
    u, v, x, kind = rnd(0, n - 1), rnd(0, n - 1), rnd(-9, 9), rnd(0, 3)
    if kind == 0:
      a[v] = x
      st.set(h.pos[v], x)
      lz.set(h.pos[v], x)
    elif kind == 1:
      for l, r in h.path(u, v):
        lz.apply(l, r, x)
      for y in br.path(u, v)[0]:
        a[y] += x
        st.set(h.pos[y], a[y])
    elif kind == 2:
      p = br.path(u, v)[0]
      s = sum(st.query(l, r) for l, r in h.path(u, v))
      check_eq(s, sum(a[y] for y in p))
      m = min(lz.query(l, r) for l, r in h.path(u, v))
      check_eq(m, min(a[y] for y in p))
      se = sum(st.query(l, r) for l, r in h.path(u, v, True))
      check_eq(se, s - a[br.lca(u, v)])
    else:
      l, r = h.subtree(v)
      inside = [a[y] for y in range(n) if br.is_anc(v, y)]
      check_eq(st.query(l, r), sum(inside))
      check_eq(lz.query(l, r), min(inside))


# Reviewer's case: walked left to right, path(3, 2) is 3 0 1 2.
adj = [[1, 2], [0, 3, 4], [0], [1], [1, 5], [4]]
check_eq(lib.HLD(adj).path(3, 2), [(4, 5), (0, 2), (5, 6)])
check_tree(adj, 0)

for it in range(iters(600)):
  n = rnd(1, 4 if it < 100 else 18)
  adj = random_tree(n, it % 5)
  root = rnd(0, n - 1)
  check_tree(adj, root)
  check_usage(adj, root, 30)

# Large: 2e5 path sums with ds/segtree on a random tree and a path
# (deep), checked against prefix sums from the root.
n = 200_000
for shape in (0, 1):
  adj = random_tree(n, shape)
  root = rnd(0, n - 1)
  a = [rnd(-(10**9), 10**9) for _ in range(n)]
  qs = [(rnd(0, n - 1), rnd(0, n - 1)) for _ in range(n)]
  with Timer(f"shape {shape}, n = 2e5: build"):
    h = lib.HLD(adj, root)
  b = [0] * n
  for v in range(n):
    b[h.pos[v]] = a[v]
  st = lib.SegTree(b, add, 0)
  with Timer(f"shape {shape}: 2e5 path sums with ds/segtree"):
    got = []
    for u, v in qs:
      s = 0
      for l, r in h.path(u, v):
        s += st.query(l, r)
      got.append(s)
  pre, q = [0] * n, [root]
  pre[root] = a[root]
  for u in q:
    for v in adj[u]:
      if v != h.par[u]:
        pre[v] = pre[u] + a[v]
        q.append(v)
  want = []
  for u, v in qs:
    w = h.lca(u, v)
    check(h.depth[w] <= min(h.depth[u], h.depth[v]))
    want.append(pre[u] + pre[v] - 2 * pre[w] + a[w])
  check_eq(got, want)
print("OK tree/hld")
