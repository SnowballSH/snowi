from collections import Counter
from itertools import permutations
from math import sqrt

from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("contest/gen")
Gen = lib.Gen
SHAPES = [Gen.UNIFORM, Gen.PARENT, Gen.PATH, Gen.STAR, Gen.CATERPILLAR]


def degrees(n, e):
  deg = [0] * n
  for u, v in e:
    deg[u] += 1
    deg[v] += 1
  return deg


def simple(n, e):
  """Endpoints in range, no loops, no unordered pair twice."""
  seen = set()
  for u, v in e:
    if not (0 <= u < n and 0 <= v < n) or u == v:
      return False
    key = (min(u, v), max(u, v))
    if key in seen:
      return False
    seen.add(key)
  return True


def components(n, e):
  parent = list(range(n))

  def find(x):
    while parent[x] != x:
      parent[x] = parent[parent[x]]
      x = parent[x]
    return x

  count = n
  for u, v in e:
    a, b = find(u), find(v)
    if a != b:
      parent[a] = b
      count -= 1
  return count


def is_tree(n, e):
  return len(e) == n - 1 and simple(n, e) and components(n, e) == 1


def acyclic(n, e):
  adj = [[] for _ in range(n)]
  indeg = [0] * n
  for u, v in e:
    adj[u].append(v)
    indeg[v] += 1
  order = [v for v in range(n) if indeg[v] == 0]
  for u in order:
    for w in adj[u]:
      indeg[w] -= 1
      if indeg[w] == 0:
        order.append(w)
  return len(order) == n


def chi2_ok(counts, cells):
  """Pearson test against `cells` equally likely cells (unseen ones
  count too); fails with probability below about 1e-6 when uniform."""
  total = sum(counts.values())
  want = total / cells
  stat = want * (cells - len(counts))
  stat += sum((c - want) ** 2 / want for c in counts.values())
  dof = cells - 1
  return stat < dof + 7 * sqrt(2 * dof) + 10


def normalized(e):
  return tuple(sorted((min(u, v), max(u, v)) for u, v in e))


def check_shape(n, shape, e):
  check(is_tree(n, e), f"shape {shape}, n {n}: not a tree: {e}")
  deg = degrees(n, e)
  if shape == Gen.PATH:
    check(max(deg) <= 2)
  if shape == Gen.STAR:
    check(max(deg) == n - 1)
  if shape == Gen.CATERPILLAR:  # the non-leaves form a path
    inner = [0] * n
    for u, v in e:
      if deg[u] > 1 and deg[v] > 1:
        inner[u] += 1
        inner[v] += 1
    check(max(inner) <= 2)


def test_numbers():
  for _ in range(iters(300)):
    g = Gen(rnd(0, 10**18))
    lo = rnd(-5, 5)
    hi = lo + rnd(0, 6)
    xs = [g.rnd(lo, hi) for _ in range(200)]
    check(all(lo <= x <= hi for x in xs))
    check(lo in xs and hi in xs)
    big = g.rnd(-(10**18), 10**18)
    check(-(10**18) <= big <= 10**18)

    n = rnd(0, 12)
    hi = lo + (n - 1 if rnd(0, 1) else rnd(n, 3 * n + 3))
    for distinct in (False, True):
      a = g.array(n, lo, hi, distinct)
      check_eq(len(a), n)
      check(all(lo <= x <= hi for x in a))
      if distinct:
        check_eq(len(set(a)), n)
    check_eq(sorted(g.perm(n)), list(range(n)))

    alphabet = "ab" if rnd(0, 1) else "xyz"
    s = g.str(n, alphabet)
    check(isinstance(s, str))
    check_eq(len(s), n)
    check(set(s) <= set(alphabet))
  g = Gen(rnd(0, 10**18))
  a = g.array(1000, -(10**18), 10**18, True)
  check_eq(len(set(a)), 1000)
  check(all(-(10**18) <= x <= 10**18 for x in a))
  s = g.str(1000)
  check(s.count("a") > 400 and s.count("b") > 400)
  check_eq(s.count("a") + s.count("b"), 1000)
  letters = Counter(g.str(5000, "vwxyz"))  # every letter, equally often
  check_eq(len(letters), 5)
  check(chi2_ok(letters, 5), f"letters {letters}")
  tuples = Counter(tuple(g.array(2, 3, 6, True)) for _ in range(24000))
  pairs = Counter(tuple(g.array(2, 0, 1)) for _ in range(24000))
  perms = Counter(tuple(g.perm(4)) for _ in range(24000))
  check(chi2_ok(tuples, 12) and chi2_ok(pairs, 4) and chi2_ok(perms, 24))


def test_extremes():
  """Ranges of 2^63 values and more (random.sample needs len())."""
  lo, hi = -(2**63), 2**63 - 1
  g = Gen(rnd(0, 10**18))
  cases = [(3, 0, hi), (1, lo, hi), (4, lo, hi), (2, lo, lo + 1)]
  cases += [(3, -5 * 10**18, 5 * 10**18), (3, hi - 2, hi), (4, 0, 2**100)]
  for n, a, b in cases:
    for distinct in (False, True):
      v = g.array(n, a, b, distinct)
      check_eq(len(v), n)
      check(all(a <= x <= b for x in v), f"{v} not in [{a}, {b}]")
      if distinct:
        check_eq(len(set(v)), n)


def test_seeds():
  for _ in range(iters(50)):
    seed = rnd(0, 10**18)
    a, b = Gen(seed), Gen(seed)
    n = rnd(20, 40)
    m = rnd(n - 1, 2 * n)
    check_eq(a.tree(n), b.tree(n))
    check_eq(a.graph(n, m), b.graph(n, m))
    check_eq(a.array(n, 0, 10**9, True), b.array(n, 0, 10**9, True))
    check_eq(a.str(n), b.str(n))
    check(Gen(seed).perm(n) != Gen(seed + 1).perm(n))


def test_trees():
  for it in range(iters(2000)):
    g = Gen(rnd(0, 10**18))
    n = 1 + it % 3 if it < 50 else rnd(1, 30)
    for shape in SHAPES:  # with m = n - 1, graph is its tree
      check_shape(n, shape, g.tree(n, shape))
      check_shape(n, shape, g.graph(n, n - 1, True, shape))
  # Uniform labelled trees: n^(n-2) of them, all equally likely.
  g = Gen(rnd(0, 10**18))
  for n, samples in ((3, 3000), (4, 16000), (5, 37500)):
    cnt = Counter(normalized(g.tree(n)) for _ in range(samples))
    check_eq(len(cnt), n ** (n - 2))
    check(chi2_ok(cnt, n ** (n - 2)), f"uniform trees, n = {n}")
  # Relabelled, shuffled and flipped: star centres and path ends are
  # uniform, the centre is first in half of the edges, and 4 of the
  # 10 edge pairs of a 5-edge path are adjacent.
  centre, ends = Counter(), Counter()
  first = total = adjacent = 0
  for _ in range(6000):
    star = g.tree(6, Gen.STAR)
    deg = degrees(6, star)
    c = deg.index(5)
    centre[c] += 1
    first += sum(u == c for u, _ in star)
    total += len(star)
    path = g.tree(6, Gen.PATH)
    deg = degrees(6, path)
    ends.update(v for v in range(6) if deg[v] == 1)
    for (a, b), (x, y) in zip(path, path[1:]):
      adjacent += len({a, b} & {x, y}) > 0
  check(chi2_ok(centre, 6) and chi2_ok(ends, 6))
  check(abs(first - total // 2) < total // 50)
  check(abs(adjacent - 6000 * 4 * 4 // 10) < 6000 * 4 // 20)
  # Random-parent trees on 4 vertices: 2 of the 6 equally likely
  # parent arrays are stars, so 1/3 of them (uniform: 4/16).
  stars = sum(max(degrees(4, g.tree(4, Gen.PARENT))) == 3 for _ in range(12000))
  check(abs(stars - 4000) < 300, f"{stars} random-parent stars")
  # Caterpillars of 8 vertices: every spine (non-leaf) size from 1
  # (a star) to 6 (a path) occurs, and most are neither (~70%).
  spine = Counter()
  for _ in range(3000):
    deg = degrees(8, g.tree(8, Gen.CATERPILLAR))
    spine[sum(d > 1 for d in deg)] += 1
  check_eq(sorted(spine), [1, 2, 3, 4, 5, 6])
  check(sum(spine[k] for k in range(2, 6)) > 1500, f"spines {spine}")


def test_graphs():
  for it in range(iters(2000)):
    g = Gen(rnd(0, 10**18))
    n = 1 + it % 4 if it < 60 else rnd(1, 14)
    pairs = n * (n - 1) // 2
    connected = rnd(0, 1) == 1
    lo = n - 1 if connected else 0
    m = (pairs, lo, rnd(lo, pairs))[it % 3]
    e = g.graph(n, m, connected)
    check_eq(len(e), m)
    check(simple(n, e))
    if connected:
      check_eq(components(n, e), 1)
    d = g.dag(n, m)
    check_eq(len(d), m)
    check(simple(n, d) and acyclic(n, d))
  # m random pairs of n = 5 vertices: all C(10, 3) sets equally
  # likely, so the pair decoding is a bijection.
  g = Gen(rnd(0, 10**18))
  sets = Counter(normalized(g.graph(5, 3, False)) for _ in range(24000))
  check(chi2_ok(sets, 120))
  # With m = n - 1 a connected graph is its spanning tree, which is
  # uniform by default: all 16 labelled trees on 4 vertices.
  spanning = Counter(normalized(g.graph(4, 3)) for _ in range(16000))
  check_eq(len(spanning), 16)
  check(chi2_ok(spanning, 16), "spanning trees are not uniform")
  # Transitive tournaments on 4 vertices: the topological order is a
  # random one of the 24.
  orders = Counter()
  for _ in range(24000):
    out = [0] * 4
    for u, _ in g.dag(4, 6):
      out[u] += 1
    orders[tuple(out)] += 1
  check(set(orders) == set(permutations(range(4))))
  check(chi2_ok(orders, 24))
  forward = sum(u < v for _ in range(2000) for u, v in g.graph(6, 9))
  check(abs(forward - 9000) < 600)


def test_large():
  g = Gen(rnd(0, 10**18))
  n = 200_000
  # Idle: trees 0.8 s, graph and dag 1.9 s, the rest 0.4 s; the
  # bounds catch quadratic code and allow for a loaded machine.
  with Timer("5 trees, n = 2e5") as t:
    trees = [g.tree(n, shape) for shape in SHAPES]
  check(t.seconds < 5, "trees are slow")
  check(all(is_tree(n, e) for e in trees))
  with Timer("graph and dag, n = 2e5, m = 4e5") as t:
    e, d = g.graph(n, 2 * n), g.dag(n, 2 * n)
  check(t.seconds < 10, "graph or dag is slow")
  check(simple(n, e) and components(n, e) == 1 and acyclic(n, d))
  k = 700
  with Timer("complete graph, 244650 edges, and a tight array") as t:
    e = g.graph(k, k * (k - 1) // 2)
    a = g.array(n, 0, n - 1, True)
  check(t.seconds < 3, "dense graph or tight array is slow")
  check_eq(len(e), k * (k - 1) // 2)
  check_eq(sorted(a), list(range(n)))


test_numbers()
test_extremes()
test_seeds()
test_trees()
test_graphs()
test_large()
print("OK contest/gen")
