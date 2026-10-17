"""
Title: Edge colouring
Description: `edge_colouring`: a bipartite multigraph with $Delta$
  colours, the maximum degree (Kőnig). `vizing`: a simple graph with
  at most $Delta + 1$ colours (Misra-Gries).
Precondition: `edge_colouring`: edges `(u, v)` with `u` in `[0, L)`
  on the left and `v` in `[0, R)` on the right; parallel edges are
  fine. `vizing`: vertices `0..n-1`, no loops or parallel edges.
Time: $O(m (n + Delta))$ time and $8 n Delta$ bytes, $n$ the number
  of vertices.
Usage:
  col = edge_colouring(L, R, edges)  # in [0, Delta)
  col2 = vizing(n, edges)  # in [0, Delta]
Example:
  >>> edge_colouring(2, 2, [(0, 0), (0, 1), (1, 0), (1, 1)])
  [1, 0, 0, 1]
  >>> vizing(3, [(0, 1), (1, 2), (2, 0)])  # odd cycle: 3 colours
  [0, 1, 2]
  >>> # A matrix with all row and column sums k as k permutation
  >>> # matrices: one edge (i, j) per unit of a[i][j].
Notes: For each edge, take a colour `a` free at `u` and `b` free at
  `v`. If they differ, flip `a` and `b` along the alternating path
  from `v` that starts with colour `a`; in a bipartite graph it
  cannot reach `u`, and afterwards `a` is free at both ends. `vizing`
  recolours a fan of neighbours of `u` and one alternating path per
  edge. Deciding whether a general graph needs $Delta$ or
  $Delta + 1$ colours is NP-hard.
Source: Kőnig's proof by alternating paths; kactl EdgeColoring
  (CC0).
Status: stress-tested against brute force
"""


def edge_colouring(L, R, edges):
  deg, ends = [0] * (L + R), []  # ends[e]: XOR of the two ends
  for u, v in edges:
    deg[u] += 1
    deg[L + v] += 1
    ends.append(u ^ (L + v))
  D = max(deg, default=0)
  at = [-1] * ((L + R) * D)  # at[x * D + c]: edge coloured c
  for e, (u, v) in enumerate(edges):
    u, v, a, b = u * D, (L + v) * D, 0, 0
    while at[u + a] >= 0:
      a += 1
    while at[v + b] >= 0:
      b += 1
    x, c = v // D, a
    while a != b and x >= 0:  # flip a and b on the path from v
      f, i, j = at[x * D + c], x * D + a, x * D + b
      at[i], at[j] = at[j], at[i]
      x, c = -1 if f < 0 else ends[f] ^ x, c ^ a ^ b
    at[u + a] = at[v + a] = e
  col = [0] * len(edges)
  for i in range(L * D):
    if at[i] >= 0:
      col[at[i]] = i % D
  return col


def vizing(n, edges):
  cc = [0] * (n + 1)
  for a, b in edges:
    cc[a] += 1
    cc[b] += 1
  K = max(cc) + 1
  adj = [-1] * (n * K)  # adj[x * K + c]: neighbour via colour c
  fr, fan = [0] * n, [0] * n
  for u, v in edges:
    fan[0], loc, ind = v, [0] * K, 0
    # fan: next is u's neighbour via the free colour
    while True:
      d = fr[v]
      if loc[d] or adj[u * K + d] < 0:
        break
      v = adj[u * K + d]
      ind += 1
      loc[d], cc[ind], fan[ind] = ind, d, v
    c = fr[u]
    cc[loc[d]] = c
    at, end, cd = u, u, d
    while at != -1:  # swap c and d on the path from u
      i, j = at * K + cd, at * K + (cd ^ c ^ d)
      adj[i], adj[j] = adj[j], adj[i]
      end, cd = at, cd ^ c ^ d
      at = adj[at * K + cd]
    i = 0
    while adj[fan[i] * K + d] != -1:  # rotate the fan
      left, right = fan[i], fan[i + 1]
      i += 1
      f = cc[i]
      adj[u * K + f], adj[left * K + f] = left, u
      adj[right * K + f], fr[right] = -1, f
    adj[u * K + d], adj[fan[i] * K + d] = fan[i], u
    for y in (fan[0], u, end):
      z = 0
      while adj[y * K + z] != -1:
        z += 1
      fr[y] = z
  res = []
  for u, v in edges:
    c = 0
    while adj[u * K + c] != v:
      c += 1
    res.append(c)
  return res
