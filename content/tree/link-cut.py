"""
Title: Link-cut tree
Description: A forest under edge insertions and deletions, with
  path queries, point updates. Each tree has a root, and
  `make_root(v)` moves it to `v` (evert). `link(u, v)` adds the
  edge `u - v` with `u` as a child of `v`; the root of `v`'s tree
  stays the root. `cut(u)` removes the edge from `u` to its parent
  (nothing at a root); `cut(u, v)` removes the edge `u - v` and
  leaves `v` as the root. `query(u, v)` folds the values on the
  path from `u` to `v` and leaves `u` as the root. `lca(u, v)` and
  `find_root(v)` use the current root. After `access(v)`, `agg[v]`
  is the fold of the root path of `v`.
Precondition: 0-indexed. `link`: `u` and `v` in different trees.
  `cut(u, v)`: the edge exists. `lca`, `query`: `u` and `v`
  connected. `op` is associative and commutative (paths get
  reversed) with identity `e`.
Time: Amortised $O(log n)$ per operation. PyPy: about 1.4 s for
  $10^5$ mixed operations (`query`, `set`, `lca`, `connected`) on
  $10^5$ vertices.
Usage:
  lc = LinkCut(a, add, 0)  # vertex values a[v], no edges
  lc.link(u, v); lc.cut(u, v); lc.connected(u, v)
  lc.query(u, v); lc.set(v, x); lc.lca(u, v)
  lc.access(v); s = lc.agg[v]  # root path; root unchanged
Tweak: Max: `LinkCut(a, max, -inf)`. Edge values: add a node per
  edge holding its value (vertices hold `e`) and link `u - x - v`.
  Non-commutative `op`: keep the reverse-order fold `ragg` next to
  `agg` everywhere (combining right child, `val`, left child) and
  swap `agg[c]` and `ragg[c]` the moment `flip[c]` is toggled (in
  the pushes and in `make_root`), not when `c` itself is pushed:
  the rotations read children with pending flips. For speed,
  inline `op`.
Example:
  >>> from operator import add
  >>> lc = LinkCut([5, 1, 4, 2, 3], add, 0)
  >>> lc.link(1, 0); lc.link(2, 1); lc.link(3, 1)  # root 0
  >>> lc.lca(2, 3), lc.find_root(3), lc.connected(2, 4)
  (1, 0, False)
  >>> _ = lc.access(3)  # root path 0 1 3: 5 + 1 + 2
  >>> lc.agg[3]
  8
  >>> lc.cut(1); lc.link(1, 4)  # move the subtree of 1 under 4
  >>> lc.lca(2, 3), lc.find_root(2), lc.query(2, 3)
  (1, 4, 7)
  >>> lc.lca(1, 4), lc.lca(3, 2)  # query made 2 the root
  (1, 2)
  >>> # Minimum spanning forest under edge insertions. Node n + i
  >>> # is edge i: value (weight, i), op = max, vertices hold e,
  >>> # which must be below every edge (weights may be negative).
  >>> from math import inf
  >>> n, e = 4, (-inf, -1)
  >>> edges = [(0, 1, 5), (1, 2, 3), (2, 0, 4), (2, 3, -2),
  ...          (3, 2, -3)]  # the last one replaces (2, 3, -2)
  >>> lc = LinkCut([e] * (n + len(edges)), max, e)
  >>> total = 0
  >>> for i, (u, v, w) in enumerate(edges):
  ...   if lc.connected(u, v):
  ...     mw, j = lc.query(u, v)  # heaviest edge on the cycle
  ...     if mw <= w: continue
  ...     lc.cut(n + j, edges[j][0]); lc.cut(n + j, edges[j][1])
  ...     total -= mw
  ...   lc.set(n + i, (w, i))
  ...   lc.link(n + i, u); lc.link(n + i, v)
  ...   total += w
  >>> total  # edges 3, 4, -3
  4
Notes: Cues: edges added and removed online with path questions
  ("are `u` and `v` connected now", "max weight on the path"), a
  spanning forest that must evict the heaviest or oldest edge of
  the cycle a new edge closes, parent pointers that change (a
  rooted forest where a subtree is moved under a new parent:
  `cut(u)`, then `link(u, p)`). For a fixed root, call
  `link(u, p)` only on a root `u` (after `cut(u)`) and never call
  `make_root`, `query` or `cut(u, v)`; use `access(v)` and
  `agg[v]` for root paths.
  Moving a subtree changes the depths of all its vertices at once:
  put that change on the new edge (a node between `u` and `p`, or
  the value of `u` itself) instead of on the vertices. To find the
  first vertex on a root path with a property of its fold (value
  at most `x`), `access(v)` and walk down from `v` using `agg` of
  the children, pushing as you go, then `access` the vertex found.
  Offline, with all edge changes known in advance and only
  connectivity or per-component questions,
  `offline/dynamic-connectivity` (DSU with rollback) is simpler;
  it knows no paths, so path queries, `lca` or a moving root still
  need this.
  Subtree folds need each node to also keep the fold of its
  virtual (path-parent) children.
Source: Sleator and Tarjan; array splay as in kactl LinkCutTree
  (CC0).
Status: stress-tested against an explicit forest
"""


class LinkCut:
  def __init__(self, a, op, e):
    n = len(a)  # node n is null and holds e
    self.n, self.op = n, op
    # the children of x are ch[2x] and ch[2x + 1]
    self.ch = [n] * (2 * n + 2)
    self.p, self.flip = [n] * (n + 1), [0] * (n + 1)
    self.val, self.agg = list(a) + [e], list(a) + [e]
    self.st = [0] * (n + 1)

  def access(self, x):
    """Makes root..x one path, x at the top of its splay tree.
    Returns where it met the root's old path: the lca with the
    vertex accessed before. Splays inline (faster in PyPy)."""
    ch, p, flip, st = self.ch, self.p, self.flip, self.st
    val, agg, op, n = self.val, self.agg, self.op, self.n
    last, nxt = n, x
    while True:
      y = x if nxt == n else nxt  # splay y
      st[0], h, z = y, 1, y
      while True:  # push flips from the top of the splay tree
        q = p[z]
        if ch[2 * q] != z and ch[2 * q + 1] != z:
          break
        st[h] = z = q
        h += 1
      while h:
        h -= 1
        z = st[h]
        if flip[z]:
          l, r = ch[2 * z], ch[2 * z + 1]
          ch[2 * z], ch[2 * z + 1] = r, l
          flip[l] ^= 1
          flip[r] ^= 1
          flip[z] = 0
      while True:
        q = p[y]
        if ch[2 * q] != y and ch[2 * q + 1] != y:
          break
        g, w, k = p[q], y, 1
        if ch[2 * g] == q or ch[2 * g + 1] == q:
          k = 2
          if (ch[2 * q] == y) == (ch[2 * g] == q):
            w = q
        while k:  # rotate w above its parent a
          a = p[w]
          b, d = p[a], 1 if ch[2 * a + 1] == w else 0
          if ch[2 * b] == a:
            ch[2 * b] = w
          elif ch[2 * b + 1] == a:
            ch[2 * b + 1] = w
          c = ch[2 * w + 1 - d]
          p[w], p[a], p[c] = b, w, a
          ch[2 * a + d], ch[2 * w + 1 - d] = c, a
          l, r = ch[2 * a], ch[2 * a + 1]
          agg[a] = op(op(agg[l], val[a]), agg[r])
          w, k = y, k - 1
      if nxt != n:
        ch[2 * y + 1] = last
      l, r = ch[2 * y], ch[2 * y + 1]
      agg[y] = op(op(agg[l], val[y]), agg[r])
      if nxt == n:
        return last
      last, nxt = y, p[y]

  def make_root(self, x):
    self.access(x)
    self.flip[x] ^= 1

  def find_root(self, x):
    ch, flip = self.ch, self.flip
    self.access(x)
    while True:
      if flip[x]:
        l, r = ch[2 * x], ch[2 * x + 1]
        ch[2 * x], ch[2 * x + 1] = r, l
        flip[l] ^= 1
        flip[r] ^= 1
        flip[x] = 0
      if ch[2 * x] == self.n:
        break
      x = ch[2 * x]
    self.access(x)
    return x

  def connected(self, u, v):
    return self.find_root(u) == self.find_root(v)

  def link(self, u, v):
    self.make_root(u)
    self.p[u] = v

  def cut(self, u, v=-1):
    if v >= 0:
      self.make_root(v)
    self.access(u)  # u is last on its path: no right child
    ch, n = self.ch, self.n
    self.p[ch[2 * u]] = n
    ch[2 * u] = n
    self.agg[u] = self.val[u]

  def query(self, u, v):
    self.make_root(u)
    self.access(v)
    return self.agg[v]

  def set(self, v, x):
    self.access(v)
    ch, agg, op = self.ch, self.agg, self.op
    self.val[v] = x
    agg[v] = op(op(agg[ch[2 * v]], x), agg[ch[2 * v + 1]])

  def lca(self, u, v):
    self.access(u)
    return self.access(v)
