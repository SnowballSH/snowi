"""
Title: Maximum clique
Description: Returns the vertices of one maximum clique (sorted).
  Branch and bound: a greedy colouring of the candidates bounds the
  clique size (a clique uses each colour at most once), and branches
  start at the highest colours.
Precondition: `adj[u]` is an int bitmask of the neighbours of `u`,
  symmetric, without bit `u` itself.
Time: Exponential. Random graphs with edge probability 0.9: PyPy
  needs about 0.2 s for $n = 100$ and 1 s for $n = 128$; sparser or
  smaller ones are far faster.
Usage:
  adj = [0] * n
  adj[u] |= 1 << v; adj[v] |= 1 << u
  clique = max_clique(adj)
Tweak: Maximum independent set: run on the complement,
  `[~adj[u] & ~(1 << u) & (1 << n) - 1 for u in range(n)]`. Masks
  up to 62 bits stay machine words in PyPy; longer ones work but
  are slower.
Example:
  >>> def graph(n, edges):
  ...   adj = [0] * n
  ...   for u, v in edges: adj[u] |= 1 << v; adj[v] |= 1 << u
  ...   return adj
  >>> # K4 on 0..3 plus the edge 3-4
  >>> max_clique(graph(5, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3),
  ...                      (2, 3), (3, 4)]))
  [0, 1, 2, 3]
  >>> # independent set on the path 0-1-2-3-4: complement
  >>> adj = graph(5, [(0, 1), (1, 2), (2, 3), (3, 4)])
  >>> max_clique([~a & ~(1 << u) & 31 for u, a in enumerate(adj)])
  [0, 2, 4]
Notes: Clique, independent set and vertex cover are one problem: $S$
  is independent iff it is a clique of the complement iff $V without
  S$ is a vertex cover, so the minimum vertex cover is $n$ minus the
  maximum independent set. Check the structure before going
  exponential: in a bipartite graph the maximum independent set is
  $n$ minus the maximum matching (Konig, `flow/hopcroft-karp`);
  trees take a DP, interval graphs a greedy, partial orders
  Dilworth. Cues: "pairwise compatible" or "pairwise conflicting",
  "largest group where everyone knows everyone", a general graph
  with $n <= 100$. For $n <= 40$, meet in the middle ($2^(n\/2)$
  subsets per half) also counts the independent sets. Dense graphs
  are the slow case, and complements of sparse inputs are dense.
Source: Tomita and Seki (MCQ) with bitsets as in San Segundo's
  BBMC; compare kactl MaxClique (CC0).
Status: stress-tested against subset enumeration
"""


def max_clique(adj):
  n = len(adj)
  pc = lambda x: bin(x).count("1")
  ord_ = sorted(range(n), key=lambda v: -pc(adj[v]))
  g = [0] * n  # relabelled by decreasing degree
  for i in range(n):
    for j in range(n):
      if adj[ord_[i]] >> ord_[j] & 1:
        g[i] |= 1 << j
  best, cur = [], []

  def expand(cand):
    nonlocal best
    vs, col, k, left = [], [], 0, cand
    while left:  # colour classes
      k += 1
      q = left
      while q:
        low = q & -q
        v = low.bit_length() - 1
        q &= ~(g[v] | low)
        left ^= low
        vs.append(v)
        col.append(k)
    for i in range(len(vs) - 1, -1, -1):
      if len(cur) + col[i] <= len(best):
        return
      v = vs[i]
      cur.append(v)
      nxt = cand & g[v]
      if nxt:
        expand(nxt)
      elif len(cur) > len(best):
        best = cur[:]
      cur.pop()
      cand &= ~(1 << v)

  expand((1 << n) - 1)
  return sorted(ord_[v] for v in best)
