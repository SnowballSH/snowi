"""
Title: Hungarian algorithm
Description: Minimum-cost assignment: gives each of the `n` rows of
  the cost matrix `a` its own column, minimizing the total cost.
  Returns the cost and `match`, where `match[i]` is the column of
  row `i`.
Precondition: `n <= m` for an `n x m` matrix (transpose if not).
  Costs may be negative, with $|a_(i j)| < 2^60$.
Time: $O(n^2 m)$.
Usage:
  cost, match = hungarian(a)  # a: list of n rows of length m
Tweak: Maximum-weight assignment: negate every entry and the
  returned cost. Rows that may stay unassigned: add $n$ dummy
  columns of cost 0. Float costs work unchanged.
Example:
  >>> hungarian([[4, 1, 3], [2, 0, 5], [3, 2, 2]])
  (5, [1, 0, 2])
  >>> hungarian([[7, 3, 9, 1], [8, 2, 6, 5]])  # n < m
  (3, [3, 1])
  >>> c, match = hungarian([[-5, -2], [-1, -4]])  # maximize
  >>> -c, match  # scores 5 + 4
  (9, [0, 1])
Notes: This is successive shortest paths on the complete bipartite
  graph with dense $O(m)$ Dijkstra steps: rows are added one at a
  time, and potentials `u` (rows) and `v` (columns) keep
  $a_(i j) - u_i - v_j >= 0$, with equality on the matched pairs.
  Think of it for "assign workers to jobs", "match $n$ points to
  $n$ targets with minimum total distance", or a permutation
  minimizing $sum_i a_(i, p(i))$. In PyPy $n = m = 300$ takes
  about a second ($2.7 dot 10^7$ inner steps). Sparse or
  capacitated variants fit min-cost flow better.
Source: e-maxx (cp-algorithms) Hungarian with potentials; kactl
  (CC0).
Status: stress-tested against permutations and planted optima
"""


def hungarian(a):
  n = len(a)
  m = len(a[0]) if n else 0
  INF = 1 << 62  # an int: a list of floats slows PyPy
  u, v = [0] * n, [0] * (m + 1)  # column m is a virtual start
  p = [-1] * (m + 1)  # p[j]: row matched to column j
  way = [0] * (m + 1)
  for i in range(n):
    dist, done = [INF] * m, [False] * (m + 1)
    j0, p[m] = m, i
    while p[j0] >= 0:  # Dijkstra step: add the closest column
      done[j0] = True
      i0, delta, j1 = p[j0], INF, 0
      row, ui = a[i0], u[i0]
      for j in range(m):
        if not done[j]:
          cur = row[j] - ui - v[j]
          if cur < dist[j]:
            dist[j], way[j] = cur, j0
          if dist[j] < delta:
            delta, j1 = dist[j], j
      for j in range(m + 1):
        if done[j]:
          u[p[j]] += delta
          v[j] -= delta
        else:
          dist[j] -= delta
      j0 = j1
    while j0 != m:
      j1 = way[j0]
      p[j0], j0 = p[j1], j1
  match = [-1] * n
  for j in range(m):
    if p[j] >= 0:
      match[p[j]] = j
  return -v[m], match
