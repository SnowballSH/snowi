from itertools import permutations

from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("tree/tree-hash")


def from_edges(n, es):
  adj = [[] for _ in range(n)]
  for u, v in es:
    adj[u].append(v)
    adj[v].append(u)
  if n <= 100:  # large trees get random labels instead
    for a in adj:
      rng.shuffle(a)
  return adj


def edges_of(adj):
  return [(u, v) for u in range(len(adj)) for v in adj[u] if u < v]


def random_tree(n, shape):
  """Shapes: 0 random, 1 path, 2 star, 3 long thin, 4 binary-ish."""
  perm = list(range(n))
  rng.shuffle(perm)
  es = []
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
    es.append((perm[i], perm[p]))
  return from_edges(n, es)


def random_perm(n):
  p = list(range(n))
  rng.shuffle(p)
  return p


def relabel(adj, perm):
  return from_edges(len(adj), [(perm[u], perm[v]) for u, v in edges_of(adj)])


def canon(adj, root):
  """Brute force over all relabellings: the smallest sorted edge
  list (root sent to 0 if root >= 0), and the automorphisms."""
  es = edges_of(adj)
  mine = sorted(es)
  best, aut = None, 0
  for p in permutations(range(len(adj))):
    m = sorted((min(p[u], p[v]), max(p[u], p[v])) for u, v in es)
    if root < 0 or p[root] == 0:
      best = m if best is None or m < best else best
    if root < 0 or p[root] == root:
      aut += m == mine
  return tuple(best), aut


def subtree(adj, r, v):
  """The subtree of v (tree rooted at r), relabelled with v -> 0."""
  par, q = {r: r}, [r]
  for x in q:
    for y in adj[x]:
      if y not in par:
        par[y] = x
        q.append(y)
  idx, st, es = {v: 0}, [v], []
  while st:
    x = st.pop()
    for y in adj[x]:
      if y != par[x]:
        idx[y] = len(idx)
        es.append((idx[x], idx[y]))
        st.append(y)
  return from_edges(len(idx), es)


def brute_centers(adj):
  n = len(adj)
  ecc = []
  for s in range(n):
    d, q = [-1] * n, [s]
    d[s] = 0
    for u in q:
      for v in adj[u]:
        if d[v] < 0:
          d[v] = d[u] + 1
          q.append(v)
    ecc.append(max(d))
  return [v for v in range(n) if ecc[v] == min(ecc)]


# One shared object; every id must match one isomorphism class.
th = lib.TreeHash()
seen = {"rooted": ({}, {}), "unrooted": ({}, {})}


def consistent(kind, c, ident):
  c2i, i2c = seen[kind]
  check(1 <= ident <= len(th.ids))
  check_eq(c2i.setdefault(c, ident), ident)
  check_eq(i2c.setdefault(ident, c), c)


def automorphisms(adj, r):
  """The header's Example, for any root r."""
  ids = th.rooted(adj, r)
  aut, cnt = 1, {}
  for v in range(len(adj)):
    if v != r:
      key = (th.par[v], ids[v])
      cnt[key] = cnt.get(key, 0) + 1
      aut *= cnt[key]
  return aut


def check_small(adj):
  n = len(adj)
  cs = sorted(th.centers(adj))
  check_eq(cs, brute_centers(adj))
  uc, uaut = canon(adj, -1)
  consistent("unrooted", uc, th.unrooted(adj))
  aut = automorphisms(adj, cs[0])
  if (
    len(cs) == 2
    and th.rooted(adj, cs[0])[cs[0]] == th.rooted(adj, cs[1])[cs[1]]
  ):
    aut *= 2
  check_eq(aut, uaut)
  r = rnd(0, n - 1)
  ids = th.rooted(adj, r)
  par, q = [-2] * n, [r]
  par[r] = -1
  for u in q:
    for v in adj[u]:
      if par[v] == -2:
        par[v] = u
        q.append(v)
  check_eq(th.par, par)
  rc, raut = canon(adj, r)
  consistent("rooted", rc, ids[r])
  check_eq(automorphisms(adj, r), raut)
  for v in range(n):
    consistent("rooted", canon(subtree(adj, r, v), 0)[0], ids[v])


def from_prufer(code, n):
  deg = [1] * n
  for x in code:
    deg[x] += 1
  es = []
  for x in code:
    v = deg.index(1)
    es.append((v, x))
    deg[v] -= 1
    deg[x] -= 1
  ends = [v for v in range(n) if deg[v] == 1]
  if len(ends) == 2:
    es.append(tuple(ends))
  return from_edges(n, es)


for it in range(iters(500)):
  n = rnd(1, 5 if it < 200 else 6)
  adj = random_tree(n, it % 5)
  check_small(adj)
  check_small(relabel(adj, random_perm(n)))
for it in range(iters(20)):
  adj = random_tree(7, it % 5)
  check_small(adj)
  check_small(relabel(adj, random_perm(7)))

# All labelled trees: non-isomorphic trees on n vertices are
# 1 1 1 2 3 6 11 (OEIS A000055), rooted 1 1 2 4 9 20 48.
unrooted_want, rooted_want = [1, 1, 1, 2, 3, 6, 11], [1, 1, 2, 4, 9, 20, 48]
for n in range(1, 8):
  un, ro = set(), set()
  code = [0] * max(0, n - 2)
  while True:
    adj = from_prufer(code, n)
    un.add(th.unrooted(adj))
    ro.add(th.rooted(adj, 0)[0])  # every rooted tree, some labelling
    i = 0
    while i < len(code):
      code[i] += 1
      if code[i] < n:
        break
      code[i] = 0
      i += 1
    if i == len(code):
      break
  check_eq(len(un), unrooted_want[n - 1])
  check_eq(len(ro), rooted_want[n - 1])

# Large: n = 2e5 random, path (two centers), complete binary tree,
# each against a relabelled copy; then a tree with fewer leaves.
n = 200_000
for shape in (0, 1, 4):
  a = random_tree(n, shape)
  perm = random_perm(n)
  b = relabel(a, perm)
  h = lib.TreeHash()
  r = rnd(0, n - 1)
  with Timer(f"shape {shape}, n = 2e5: rooted + unrooted"):
    ia = h.rooted(a, r)
    ua = h.unrooted(a)
  ib = h.rooted(b, perm[r])
  check_eq([ib[perm[v]] for v in range(n)], ia)
  check_eq(h.unrooted(b), ua)
  check_eq(sorted(perm[c] for c in h.centers(a)), sorted(h.centers(b)))
  if shape == 1:
    check_eq(len(h.centers(a)), 2)
  if shape == 0:  # move a leaf of x (degree >= 4) to another leaf
    x, leaf = next(
      (v, c) for v in range(n) if len(a[v]) >= 4 for c in a[v] if len(a[c]) == 1
    )
    y = next(v for v in range(n) if len(a[v]) == 1 and v != leaf)
    es = [e for e in edges_of(a) if set(e) != {x, leaf}] + [(y, leaf)]
    check(h.unrooted(from_edges(n, es)) != ua)
print("OK tree/tree-hash")
