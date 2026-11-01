"""
Title: Global minimum cut
Description: Stoer-Wagner: the cheapest way to split an undirected
  weighted graph into two non-empty parts. Returns the cut weight
  and the vertices of one part.
Precondition: Adjacency matrix with `w[u][v] = w[v][u] >= 0` (sum
  parallel edges; the diagonal is ignored), $n >= 2$. `w` is not
  modified.
Time: $O(n^3)$; PyPy needs about 0.4 s for $n = 500$.
Usage:
  w = [[0] * n for _ in range(n)]
  w[u][v] += c; w[v][u] += c
  cut, side = global_min_cut(w)
Tweak: Sparse graphs: adjacency dicts and a max-heap of `key`
  (lazy deletion) give $O(n m log n)$.
Example:
  >>> w = [[0] * 6 for _ in range(6)]
  >>> for u, v, c in [(0, 1, 3), (1, 2, 3), (2, 0, 3), (2, 3, 1),
  ...                 (3, 4, 2), (4, 5, 2), (5, 3, 2)]:
  ...   w[u][v] += c; w[v][u] += c
  >>> cut, side = global_min_cut(w)
  >>> cut, sorted(side) in ([0, 1, 2], [3, 4, 5])
  (1, True)
  >>> # Edge connectivity of a 4-cycle: two roads must close.
  >>> w = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 0]]
  >>> global_min_cut(w)[0]
  2
Notes: Each phase adds vertices one by one, always the one most
  tightly connected to those already added; the last one, $t$, is
  separated from the one before, $s$, by exactly its `key`, which
  is a minimum $s$-$t$ cut. Then $s$ and $t$ merge. Cues:
  "disconnect the graph" with no $s$ and $t$ given, "fewest edges
  whose removal splits the network", "split into two non-empty
  groups, minimizing the weight between them". The same value is
  the minimum over $t$ of the $s$-$t$ max flow for any fixed $s$
  ($n - 1$ flows), and the lightest edge of a Gomory-Hu tree.
  Directed graphs have no Stoer-Wagner: take the minimum of the
  flows from a fixed $s$ to every $t$ and back, $2(n - 1)$ flows.
  To count or list every global minimum cut: see All minimum cuts
  in the modelling note.
Source: Stoer and Wagner (1997); kactl GlobalMinCut (CC0).
Status: stress-tested against enumeration of all cuts
"""


def global_min_cut(w):
  n = len(w)
  w = [row[:] for row in w]
  group = [[v] for v in range(n)]  # vertices merged into v
  alive, best, side = list(range(n)), 0, []
  while len(alive) > 1:
    key, rest, s, t = [0] * n, alive[:], -1, -1
    while rest:
      i, top = 0, key[rest[0]]
      for j in range(1, len(rest)):
        if key[rest[j]] > top:
          i, top = j, key[rest[j]]
      s, t = t, rest[i]
      rest[i] = rest[-1]
      rest.pop()
      row = w[t]
      for v in rest:
        key[v] += row[v]
    if not side or top < best:
      best, side = top, group[t]
    group[s] += group[t]
    ws, wt = w[s], w[t]
    for v in alive:
      ws[v] += wt[v]
      w[v][s] = ws[v]
    alive.remove(t)
  return best, side
