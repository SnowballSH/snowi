"""
Title: Strongly connected components
Description: Tarjan's algorithm, iterative. Returns `comp` and the
  number of components `k`. Components are numbered in reverse
  topological order: an edge `u -> v` between different components
  has `comp[u] > comp[v]`, so component 0 is a sink and
  `k - 1, ..., 0` is a topological order of the condensation.
  Edges needed to make the graph strongly connected: $max(s, t)$
  for $s$ sources and $t$ sinks of the condensation, 0 if $k = 1$.
Precondition: `adj[u]` lists the heads of the edges leaving `u`;
  vertices 0-indexed. Self-loops and parallel edges are fine.
Time: $O(n + m)$, no recursion; PyPy needs about 1 s for
  $n = m = 5 dot 10^5$.
Usage:
  comp, k = scc(adj)
Example:
  >>> adj = [[1], [2], [0, 3], [4], [5], [3], []]
  >>> comp, k = scc(adj); comp, k
  ([1, 1, 1, 0, 0, 0, 2], 3)
  >>> # condensation DAG without duplicate edges
  >>> dag = [set() for _ in range(k)]
  >>> for u in range(len(adj)):
  ...   for v in adj[u]:
  ...     if comp[u] != comp[v]: dag[comp[u]].add(comp[v])
  >>> dag
  [set(), {0}, set()]
Notes: Think SCC when the statement says "mutually reachable", "can
  return to", "cycles can be traversed many times" (collapse each
  component, then DP on the DAG: e.g. maximum total value collected
  on a walk), or for 2-SAT. A vertex reachable from all others
  exists iff the condensation has exactly one sink. In
  the DP, process components in increasing id for "what can be
  reached from here" and in decreasing id for "what reaches here".
Source: kactl SCC (CC0), made iterative.
Status: stress-tested against brute-force reachability
"""


def scc(adj):
  n = len(adj)
  comp, tin, low, it = [-1] * n, [-1] * n, [0] * n, [0] * n
  st, t, k = [], 0, 0
  for r in range(n):
    if tin[r] >= 0:
      continue
    cs = [r]
    st.append(r)
    tin[r] = low[r] = t
    t += 1
    while cs:
      u = cs[-1]
      out, i = adj[u], it[u]
      if i < len(out):
        it[u] = i + 1
        v = out[i]
        if tin[v] < 0:
          cs.append(v)
          st.append(v)
          tin[v] = low[v] = t
          t += 1
        elif comp[v] < 0 and tin[v] < low[u]:
          low[u] = tin[v]
        continue
      cs.pop()
      if cs and low[u] < low[cs[-1]]:
        low[cs[-1]] = low[u]
      if low[u] == tin[u]:
        while True:
          v = st.pop()
          comp[v] = k
          if v == u:
            break
        k += 1
  return comp, k
