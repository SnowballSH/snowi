"""
Title: Gomory-Hu tree
Description: Builds a tree on the same vertices such that the
  minimum cut between any `u` and `v` equals the lightest edge on
  their tree path; removing that edge splits the vertices into the
  two sides of such a cut. Returns `n - 1` edges `(u, v, w)`.
Precondition: Undirected edges `(u, v, c)`, 0-indexed, `c >= 0`;
  parallel edges and self-loops are fine.
Time: $n - 1$ max flows; PyPy takes about 0.5 s for $n = 500$,
  $m = 5000$.
Usage:
  tree = gomory_hu(n, edges)  # edges: (u, v, capacity)
  for u, v, w in tree: ...  # w = min cut between u and v
Tweak: The tree is rooted at 0: `tree[v - 1]` holds `v`, its
  parent and the cut, ready for binary lifting of path minima.
Uses: flow/dinic
Example:
  >>> # A square 0-1-2-3 with a light diagonal 1-3
  >>> e = [(0, 1, 3), (1, 2, 2), (2, 3, 4), (3, 0, 1), (1, 3, 1)]
  >>> tree = gomory_hu(4, e)
  >>> tree  # the path 0 -4- 1 -4- 2 -6- 3
  [(1, 0, 4), (2, 1, 4), (3, 2, 6)]
  >>> # 1 and 3: lightest edge on 1-2-3 is 4; the edge (2, 1)
  >>> # splits {0, 1} from {2, 3}: 1 + 2 + 1 = 4.
  >>> # Sum of the minimum cut over all pairs: join tree edges
  >>> # heaviest first; sets of sizes a and b give a * b pairs.
  >>> comp, total = list(range(4)), 0
  >>> for u, v, w in sorted(tree, key=lambda x: -x[2]):
  ...   a, b = comp[u], comp[v]
  ...   total += w * comp.count(a) * comp.count(b)
  ...   comp = [a if c == b else c for c in comp]
  >>> total
  26
Notes: Gusfield's method: no graph contraction, each step one max
  flow on the original graph between `s` and its current tree
  parent. Answers all $binom(n, 2)$ pairwise min cuts, that is
  max flows (with unit capacities: numbers of edge-disjoint
  paths), with $n - 1$ flows; the lightest tree edge is the
  global minimum cut. Cues: "for every pair", "$q$ queries of
  the min cut between $u$ and $v$", "sum or count of pairs by
  connectivity". PyPy is about 6 times slower than C++:
  $n = 1000$, $m = 10^4$ takes 1.5 s and $n = 3000$, $m = 4500$
  takes 3 to 6 s, so for $n$ near 3000 use C++ (about 1 s).
  Undirected only: directed graphs have no such tree.
  Queries: min on a tree path (binary lifting), or offline in the
  heaviest-first DSU order above.
Source: Gomory and Hu (1961); Gusfield (1990), "Very simple
  methods for all pairs network flow analysis".
Status: stress-tested against cut enumeration and pairwise flows
"""


def gomory_hu(n, edges):
  d = Dinic(n)
  for u, v, c in edges:
    d.add_edge(u, v, c, c)
  par, cut = [0] * n, [0] * n
  for s in range(1, n):
    d.cap[:] = d.orig
    t = par[s]
    cut[s] = d.max_flow(s, t)
    side = d.lvl  # >= 0 on the source side
    for i in range(n):
      if i != s and par[i] == t and side[i] >= 0:
        par[i] = s
    if side[par[t]] >= 0:
      par[s], par[t] = par[t], s
      cut[s], cut[t] = cut[t], cut[s]
  return [(v, par[v], cut[v]) for v in range(1, n)]
