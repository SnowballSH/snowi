"""
Title: Tree isomorphism
Description: Canonical ids for rooted trees (AHU): a vertex's id
  names the sorted tuple of its children's ids, so two rooted
  trees get the same id iff they are isomorphic. `rooted(adj, r)`
  returns the id of every subtree, rooted at `r` (equal ids mean
  isomorphic subtrees). `unrooted(adj)` is the smallest id over
  the one or two centers, equal for two trees iff they are
  isomorphic. `centers(adj)` gives the middle one or two vertices
  of every longest path. `par` holds the parents from the last
  `rooted` call. Ids are 1, 2, ... in order of first appearance.
Precondition: `adj` holds both directions of the $n - 1$ edges of
  a tree on vertices `0..n-1`, $n >= 1$. Compare only ids from the
  same `TreeHash` object.
Time: $O(n log n)$ per call; memory $O(n)$ per tree hashed. PyPy,
  $n = 2 dot 10^5$: about 0.3 s per call.
Usage:
  th = TreeHash()  # one object for all trees compared
  ids = th.rooted(adj, r)  # ids[v]: subtree of v
  iso = th.unrooted(adj1) == th.unrooted(adj2)
Tweak: Coloured vertices: `key = (col[q[i]],) + tuple(kids)`. Edge
  labels: sort `(label, child id)` pairs instead. Forests: hash
  each tree, then sort the ids.
Example:
  >>> th = TreeHash()
  >>> a = [[1, 2], [0], [0, 3], [2]]  # path 1 - 0 - 2 - 3
  >>> b = [[1], [0, 2], [1, 3], [2]]  # path 0 - 1 - 2 - 3
  >>> th.unrooted(a) == th.unrooted(b), th.centers(a)
  (True, [0, 2])
  >>> th.rooted(a, 0)[0] == th.rooted(b, 0)[0]  # middle vs end
  False
  >>> star = [[1, 2, 3], [0], [0], [0]]
  >>> th.unrooted(star) == th.unrooted(a)
  False
  >>> # complete binary tree: distinct subtrees, automorphisms
  >>> t = [[1, 2], [0, 3, 4], [0, 5, 6], [1], [1], [2], [2]]
  >>> ids = th.rooted(t, 0)
  >>> len(set(ids))  # leaf, cherry, whole tree
  3
  >>> # Fixing root r: siblings with equal ids can be permuted, k
  >>> # equal ones in k! ways. Unrooted: r = a center; with two
  >>> # centers c, d, double it if their rooted ids are equal.
  >>> aut, seen = 1, {}
  >>> for v in range(len(t)):
  ...   if v != 0:
  ...     key = (th.par[v], ids[v])
  ...     seen[key] = seen.get(key, 0) + 1
  ...     aut *= seen[key]
  >>> aut  # swap 1 and 2, swap 3 and 4, swap 5 and 6
  8
Notes: Cues: "are these trees the same up to relabelling", "number
  of different subtrees", "which roots give the same rooted tree"
  (group roots by `rooted(adj, r)[r]` in $O(n^2 log n)$, or by
  rerooted hashes below in $O(n)$), symmetric trees, counting
  automorphisms (Example) or distinct labellings ($n!$ over the
  automorphisms). An isomorphism maps centers to centers, hence
  the rooting at a center; two centers are the middle edge of
  every diameter. Random hashes instead of exact ids: $h(v) =
  product_c (x_k + h(c))$ mod a large prime $p$, with $k$ the
  height of `v` and a random $x_k$ per height; two different trees
  collide with probability about $n \/ p$ (Schwartz-Zippel). The
  sum form $h(v) = 1 + sum_c f(h(c))$, with $f$ a random 64-bit
  mixer (splitmix64 of $h$ plus a random constant), is a `merge`
  (+), `lift` ($f$) and `fin` (+1), so `tree/rerooting` gives the
  hash for every root in $O(n)$. Never use a plain sum or XOR of
  child hashes: different trees collide easily.
Source: Aho, Hopcroft, Ullman; centers by peeling leaves.
Status: stress-tested against isomorphism by permutations
"""


class TreeHash:
  def __init__(self):
    self.ids, self.par = {}, []

  def rooted(self, adj, r):
    n = len(adj)
    par, q, lo = [-1] * n, [r], [n] * (n + 1)
    for i in range(n):  # BFS; kids of q[i]: q[lo[i]:lo[i + 1]]
      u = q[i]
      lo[i] = len(q)
      for v in adj[u]:
        if v != par[u]:
          par[v] = u
          q.append(v)
    ids, res = self.ids, [0] * n
    for i in range(n - 1, -1, -1):
      kids = sorted([res[c] for c in q[lo[i] : lo[i + 1]]])
      key = tuple(kids)
      x = ids.get(key)
      if x is None:
        x = ids[key] = len(ids) + 1
      res[q[i]] = x
    self.par = par
    return res

  def centers(self, adj):
    left, deg = len(adj), [len(a) for a in adj]
    q = [v for v in range(left) if deg[v] <= 1]
    while left > 2:  # peel all leaves
      left -= len(q)
      nq = []
      for v in q:
        for u in adj[v]:
          deg[u] -= 1
          if deg[u] == 1:
            nq.append(u)
      q = nq
    return q

  def unrooted(self, adj):
    cs = self.centers(adj)
    return min(self.rooted(adj, c)[c] for c in cs)
