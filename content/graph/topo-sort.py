"""
Title: Topological sort
Description: Kahn's algorithm: orders the vertices so that every
  edge `u -> v` has `u` before `v`. Returns `None` if the graph has
  a cycle.
Precondition: `adj[u]` lists the heads of the edges leaving `u`;
  vertices 0-indexed. Parallel edges are fine; a self-loop is a
  cycle.
Time: $O(n + m)$.
Usage:
  order = topo_sort(adj)  # None: cycle
Tweak: Lexicographically smallest order: keep the ready vertices in
  a heap (`heappush`, `heappop`) instead of scanning `order`,
  $O(n log n + m)$. To keep the acyclic part, return `order`
  anyway: it holds exactly the vertices not reachable from a cycle.
Example:
  >>> topo_sort([[1, 2], [3], [1, 3], []])
  [0, 2, 1, 3]
  >>> topo_sort([[1], [2], [0]]) is None
  True
  >>> # longest path in a DAG by DP in topological order
  >>> adj = [[1, 2], [3], [1, 3], []]
  >>> best = [0] * 4
  >>> for u in topo_sort(adj):
  ...   for v in adj[u]: best[v] = max(best[v], best[u] + 1)
  >>> best
  [0, 2, 1, 3]
Notes: Use it for "order tasks with prerequisites", for DP on a DAG
  (longest path, number of paths, reachability bitsets), and to
  detect cycles. The order is unique iff the queue never holds two
  vertices at once, i.e. consecutive vertices are joined by an edge.
  "Smallest label as early as possible" (vertex 0 first if it can
  be, then 1, ...) is not the lexicographically smallest order:
  reverse all edges, take the lexicographically largest order with a
  max-heap, and reverse the result. Constraints "`a` before `b`"
  from comparisons (sorting by partial information, dictionary order
  of an alien alphabet) are a topological sort.
Source: kactl TopoSort (CC0).
Status: stress-tested against brute force
"""


def topo_sort(adj):
  n = len(adj)
  indeg = [0] * n
  for out in adj:
    for v in out:
      indeg[v] += 1
  order = [u for u in range(n) if indeg[u] == 0]
  for u in order:  # also visits vertices appended below
    for v in adj[u]:
      indeg[v] -= 1
      if indeg[v] == 0:
        order.append(v)
  return order if len(order) == n else None
