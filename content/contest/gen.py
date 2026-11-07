"""
Title: Random test generators
Description: Seeded random inputs. Edge lists are 0-indexed,
  randomly relabelled, shuffled and oriented (`mix(e, n)` does
  this to any list); `dag` orients along a random order; a
  connected `graph` contains `tree(n, shape)`. Shapes: `UNIFORM`
  (Prüfer), `PARENT`, `PATH`, `STAR`, `CATERPILLAR`. `UNIFORM`
  and `PARENT` trees are shallow (depth about 1100 and 40 at
  $n = 2 dot 10^5$): also test `PATH` and `CATERPILLAR`.
Precondition: `tree`, `graph`: `n >= 1`. `array(n, lo, hi,
  True)`: `hi - lo + 1 >= n`. `graph`, `dag`: `m <= n(n-1)/2`,
  and `m >= n - 1` if connected.
Time: $O((n + m) log n)$.
Usage:
  import sys
  g = Gen(int(sys.argv[1]))
  n = g.rnd(1, 6)
  m = g.rnd(n - 1, n * (n - 1) // 2)
  print(n, m)
  for u, v in g.graph(n, m):
    print(u + 1, v + 1, g.rnd(1, 9))  # weight: drop if unused
Example:
  t = g.tree(n, Gen.CATERPILLAR)  # or PATH, STAR, ...
  deep = g.graph(n, m, True, Gen.PATH)  # long paths
  d = g.dag(n, m)  # edges u -> v
  a = g.array(n, 1, 3)  # small range: many ties
  b = g.array(n, 1, 10**9, True)  # distinct
  p = g.perm(n)
  s = g.str(n, "abc")
  e = [(i, (i - 1) // 2) for i in range(1, n)]  # own shape,
  e = g.mix(e, n)  # e.g. a binary tree, then mix
  tests = g.rnd(1, 3)  # multi-test: print it, then the cases
Notes: `UNIFORM` and `PARENT` trees (rooted at vertex 0) miss
  deep-recursion crashes and the inputs where a wrong heavy child
  (HLD) or merge direction (small to large) is slow: use `PATH`
  and `CATERPILLAR`, also as spanning trees of `graph`; a binary
  tree (Example) has the most light edges. `PARENT` hangs each
  vertex below a random earlier one; `CATERPILLAR` is a path
  with leaves on it. Small value ranges and alphabets force
  ties, repeated keys and periodic strings.
  Keep `n` tiny (about 8) so the brute force is fast and a
  failing input is short to read. Multi-test input: print `T`,
  then each case from the same `g`, keeping the sum of sizes
  within the limits.
Status: stress-tested: tree shapes, simple, connected and acyclic
  graphs, 64-bit ranges, chi-square uniformity of trees,
  spanning trees, letters, edge sets and orders
"""


class Gen:
  UNIFORM, PARENT, PATH, STAR, CATERPILLAR = range(5)

  def __init__(self, seed):
    import random

    self.r = random.Random(seed)
    self.rnd = self.r.randint  # uniform in [lo, hi]

  def array(self, n, lo, hi, distinct=False):
    if not distinct:
      return [self.rnd(lo, hi) for _ in range(n)]
    s = set()  # Floyd's sampling without replacement
    for j in range(hi - n + 1, hi + 1):
      x = self.rnd(lo, j)
      s.add(j if x in s else x)
    return self.r.sample(list(s), n)  # in random order

  def perm(self, n):
    return self.r.sample(range(n), n)

  def mix(self, e, n, acyclic=False):
    lab, out = self.perm(n), []
    for u, v in e:
      if (u > v) if acyclic else self.r.getrandbits(1):
        u, v = v, u
      out.append((lab[u], lab[v]))
    self.r.shuffle(out)
    return out

  def tree(self, n, shape=UNIFORM):
    import heapq

    e, rnd = [], self.rnd
    if shape == Gen.UNIFORM:  # random Pruefer code, then n - 1
      code, deg = [rnd(0, n - 1) for _ in range(n - 2)], [1] * n
      for x in code:
        deg[x] += 1
      q = [i for i in range(n) if deg[i] == 1]  # sorted: a heap
      for v in code:
        e.append((heapq.heappop(q), v))
        deg[v] -= 1
        if deg[v] == 1:
          heapq.heappush(q, v)
      if n > 1:
        e.append((q[0], n - 1))
    else:  # path 0..k-1, the others below earlier vertices
      k = n if shape == Gen.PATH else 1
      if shape == Gen.CATERPILLAR:
        k = rnd(1, n)
      for i in range(1, n):
        below = i if shape == Gen.PARENT else k
        e.append((i, i - 1 if i < k else rnd(0, below - 1)))
    return self.mix(e, n)

  def graph(self, n, m, connected=True, shape=UNIFORM):
    e = self.tree(n, shape) if connected else []
    s = set(e)  # spanning tree edges
    # pair t = {u, u + d mod n}, u = t % n, d = t // n + 1
    for t in self.array(m, 0, n * (n - 1) // 2 - 1, True):
      u, v = t % n, (t % n + t // n + 1) % n
      if len(e) < m and (u, v) not in s and (v, u) not in s:
        e.append((u, v))
    return self.mix(e, n)

  def dag(self, n, m):
    return self.mix(self.graph(n, m, False), n, True)

  def str(self, n, alphabet="ab"):
    return "".join(self.r.choices(alphabet, k=n))
