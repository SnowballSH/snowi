from math import inf
from operator import add

from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("tree/link-cut")
NEG = -(1 << 62)
MSF_E = (-inf, -1)  # the Example's identity


def sum_max(a, b):
  return a[0] + b[0], max(a[1], b[1])


class Brute:
  """An explicit forest: edge sets, one marked root per tree."""

  def __init__(self, a):
    self.n, self.val = len(a), list(a)
    self.adj = [set() for _ in a]
    self.is_root = [True] * len(a)

  def comp(self, s):
    seen, q = {s}, [s]
    for u in q:
      for v in self.adj[u]:
        if v not in seen:
          seen.add(v)
          q.append(v)
    return q

  def root(self, v):
    return next(x for x in self.comp(v) if self.is_root[x])

  def path(self, u, v):
    par, q = {u: u}, [u]
    for x in q:
      for y in self.adj[x]:
        if y not in par:
          par[y] = x
          q.append(y)
    p = [v]
    while p[-1] != u:
      p.append(par[p[-1]])
    return p[::-1]

  def make_root(self, v):
    self.is_root[self.root(v)] = False
    self.is_root[v] = True

  def link(self, u, v):
    self.is_root[self.root(u)] = False
    self.adj[u].add(v)
    self.adj[v].add(u)

  def cut(self, u):
    r = self.root(u)
    if r != u:
      p = self.path(u, r)[1]
      self.adj[u].discard(p)
      self.adj[p].discard(u)
      self.is_root[u] = True

  def lca(self, u, v):
    a = self.path(u, self.root(u))
    b = set(self.path(v, self.root(v)))
    return next(x for x in a if x in b)


def stress(n, ops):
  """Random operations; path sums and (sum, max) pairs at once."""
  a = [rnd(-(10**9), 10**9) for _ in range(n)]
  lc = lib.LinkCut(a, add, 0)
  ls = lib.LinkCut([(x, x) for x in a], sum_max, (0, NEG))
  br = Brute(a)

  def fold(p):
    vals = [br.val[x] for x in p]
    return sum(vals), max(vals)

  for _ in range(ops):
    u, v, kind = rnd(0, n - 1), rnd(0, n - 1), rnd(0, 9)
    same = br.root(u) == br.root(v)
    if kind <= 1:  # link, or cut an edge when already connected
      if not same:
        lc.link(u, v)
        ls.link(u, v)
        br.link(u, v)
      elif u != v:
        p = br.path(u, v)
        k = rnd(0, len(p) - 2)
        x, y = p[k], p[k + 1]
        lc.cut(x, y)
        ls.cut(x, y)
        br.make_root(y)
        br.cut(x)
    elif kind == 2:  # cut u from its parent
      lc.cut(u)
      ls.cut(u)
      br.cut(u)
    elif kind == 3:
      lc.make_root(u)
      ls.make_root(u)
      br.make_root(u)
    elif kind == 4:
      x = rnd(-(10**9), 10**9)
      lc.set(u, x)
      ls.set(u, (x, x))
      br.val[u] = x
    elif kind == 5:
      check_eq(lc.connected(u, v), same)
      check_eq(ls.connected(u, v), same)
      check_eq(lc.find_root(u), br.root(u))
      check_eq(ls.find_root(v), br.root(v))
    elif kind == 6 and same:
      check_eq(lc.lca(u, v), br.lca(u, v))
      check_eq(ls.lca(u, v), br.lca(u, v))
    elif kind == 7 and same:  # path fold; the root becomes u
      want = fold(br.path(u, v))
      check_eq(lc.query(u, v), want[0])
      check_eq(ls.query(u, v), want)
      br.make_root(u)
    elif kind == 8:  # root path; the root stays
      want = fold(br.path(br.root(u), u))
      lc.access(u)
      ls.access(u)
      check_eq(lc.agg[u], want[0])
      check_eq(ls.agg[u], want)
    elif kind == 9:
      check_eq(lc.val[u], br.val[u])


def check_msf(n, edges):
  """The header's Example against Kruskal on the edges so far."""
  m = len(edges)
  lc = lib.LinkCut([MSF_E] * (n + m), max, MSF_E)
  total = 0

  def insert(i):  # the Example's loop body; returns the change
    u, v, w = edges[i]
    change = 0
    if lc.connected(u, v):
      mw, j = lc.query(u, v)
      if mw <= w:
        return 0
      lc.cut(n + j, edges[j][0])
      lc.cut(n + j, edges[j][1])
      change -= mw
    lc.set(n + i, (w, i))
    lc.link(n + i, u)
    lc.link(n + i, v)
    return change + w

  for i in range(m):
    total += insert(i)
    dsu = list(range(n))

    def find(x, dsu=dsu):
      while dsu[x] != x:
        dsu[x] = dsu[dsu[x]]
        x = dsu[x]
      return x

    want = 0
    for a, b, c in sorted(edges[: i + 1], key=lambda e: e[2]):
      if find(a) != find(b):
        dsu[find(a)] = find(b)
        want += c
    check_eq(total, want)
    for _ in range(3):
      x, y = rnd(0, n - 1), rnd(0, n - 1)
      check_eq(lc.connected(x, y), find(x) == find(y))


# n = 1: a lone vertex
lc = lib.LinkCut([7], add, 0)
check_eq((lc.query(0, 0), lc.lca(0, 0)), (7, 0))
lc.cut(0)
lc.set(0, -3)
check_eq((lc.find_root(0), lc.query(0, 0)), (0, -3))

for it in range(iters(700)):
  stress(rnd(1, 5 if it < 250 else 15), 60)
# Weights below the identity's -1 once returned the vertex value as
# the heaviest edge (j = -1).
check_msf(3, [(0, 1, -5), (1, 2, -3), (2, 0, -4)])
for _ in range(iters(150)):
  n, w = rnd(1, 8), rnd(1, 20)
  check_msf(n, [(rnd(0, n - 1), rnd(0, n - 1), rnd(-w, w)) for _ in range(25)])

# Large 1: a path 0 - 1 - ... - n-1 (deep) with values in a Fenwick
# tree; random everts, cuts and relinks of path edges. Path sums are
# interval sums and the lca with root r is r clamped to [u, v].
n = 100_000
a = [rnd(-(10**9), 10**9) for _ in range(n)]
fw = [0] * (n + 1)


def fw_add(i, x):
  i += 1
  while i <= n:
    fw[i] += x
    i += i & -i


def fw_pre(i):
  s = 0
  while i > 0:
    s += fw[i]
    i -= i & -i
  return s


for i in range(n):
  fw_add(i, a[i])
ops = [(rnd(0, 4), rnd(0, n - 1), rnd(0, n - 1)) for _ in range(n)]
with Timer("path of 1e5: build + 1e5 operations"):
  lc = lib.LinkCut(a, add, 0)
  for i in range(1, n):
    lc.link(i - 1, i)
  root = n - 1
  for kind, u, v in ops:
    if kind == 0:
      x = rnd(-(10**9), 10**9)
      lc.set(u, x)
      fw_add(u, x - a[u])
      a[u] = x
    elif kind == 1:
      lc.make_root(u)
      root = u
    elif kind == 2 and u + 1 < n:
      lc.cut(u, u + 1)
      root = u + 1
      check(not lc.connected(u, u + 1))
      lc.link(u, u + 1)
    elif kind == 3:
      check_eq(lc.lca(u, v), min(max(root, min(u, v)), max(u, v)))
    else:
      want = fw_pre(max(u, v) + 1) - fw_pre(min(u, v))
      check_eq(lc.query(u, v), want)
      root = u

# Large 2: a random tree rooted at 0 whose subtrees move under new
# parents (fixed root); lca and root path sums by climbing parents.
a = [rnd(-(10**9), 10**9) for _ in range(n)]
par = [-1] + [rnd(0, i - 1) for i in range(1, n)]
mark = [-1] * n
ops = [(rnd(1, n - 1), rnd(0, n - 1)) for _ in range(n)]
with Timer("random tree of 1e5: build + 1e5 operations"):
  lc = lib.LinkCut(a, add, 0)
  for i in range(1, n):
    lc.link(i, par[i])
  for it, (u, v) in enumerate(ops):
    y = u
    while y >= 0:
      mark[y] = it
      y = par[y]
    if it % 2 == 0:  # move u's subtree under v
      x = v
      while x >= 0 and x != u:
        x = par[x]
      if x == u:  # v is in the subtree of u
        continue
      lc.cut(u)
      lc.link(u, v)
      par[u] = v
    else:
      x = v
      while mark[x] != it:
        x = par[x]
      check_eq(lc.lca(u, v), x)
      s, y = 0, v
      while y >= 0:
        s += a[y]
        y = par[y]
      lc.access(v)
      check_eq(lc.agg[v], s)
check_eq(lc.find_root(rnd(0, n - 1)), 0)
print("OK tree/link-cut")
