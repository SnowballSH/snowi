"""
Title: Min-cost max flow
Description: Sends as much flow as possible from `s` to `t` (at
  most `limit`), and among maximum flows one of minimum total cost.
  Returns `(flow, cost)`. `flow(i)` is the flow on an edge.
Precondition: Vertices 0-indexed, `s != t`, `cap >= 0`. Costs may
  be negative but no cycle of positive capacity may have negative
  total cost, and $8 n max|"cost"| < 2^62$.
Time: $O(F m log m)$ for $F$ augmenting paths (at most the flow
  value), plus $O(n m)$ when costs are negative.
Usage:
  mcf = MinCostFlow(n)
  i = mcf.add_edge(u, v, cap, cost)
  f, c = mcf.max_flow(s, t)  # or max_flow(s, t, k)
  mcf.flow(i)
Tweak: Cheapest flow of any size: put
  `if pot[t] - pot[s] >= 0: break` right after the `pot` update.
  That difference is the cost of the path about to be used (in
  the loop condition it would be stale). Cheapest circulation with
  negative cycles: saturate every negative edge first (its
  residual reverse costs $> 0$), then route the excesses from a
  super source to a super sink. The heap holds ints `d * n + v`;
  past $2^63$ PyPy switches to slower big ints.
Example:
  >>> mcf = MinCostFlow(4)
  >>> a = mcf.add_edge(0, 1, 2, 1); b = mcf.add_edge(0, 2, 1, 5)
  >>> c = mcf.add_edge(1, 3, 1, 1); d = mcf.add_edge(1, 2, 1, 1)
  >>> e = mcf.add_edge(2, 3, 2, 1)
  >>> mcf.max_flow(0, 3)  # paths 0-1-3, 0-1-2-3, 0-2-3
  (3, 11)
  >>> [mcf.flow(i) for i in (a, b, c, d, e)]
  [2, 1, 1, 1, 2]
  >>> m2 = MinCostFlow(3)
  >>> _ = m2.add_edge(0, 1, 5, -2); _ = m2.add_edge(1, 2, 5, 3)
  >>> m2.max_flow(0, 2, 4)  # limit 4 units, negative cost edge
  (4, 4)
Notes: Model: capacities say how often something can be used,
  costs what each use costs. Cues: "minimum total cost to
  transport", "assign with capacities and prices", "choose $k$
  disjoint paths or cycles", or a weighted matching too sparse
  for Hungarian. A convex cost per unit on one edge (1st unit
  costs 1, 2nd costs 3, ...) becomes parallel unit edges with
  increasing costs. Each augmenting path costs a Dijkstra:
  $F = 1000$ paths over $2.4 dot 10^4$ arcs take about 0.4 s in
  C++ (1 to 2 s in PyPy), growing like $F m log m$. Dijkstra
  stops once `t` is settled; unsettled vertices gain only $d_t$
  in potential, which keeps reduced costs $>= 0$. On a path or tree the transport cost
  has a closed form: every edge carries $|$imbalance of one
  side$|$.
Source: Classic successive shortest paths with Johnson potentials;
  kactl MinCostMaxFlow (CC0); early stop at `t` as in the AtCoder
  Library (CC0).
Status: stress-tested against exhaustive flows and Bellman-Ford SSP
"""


class MinCostFlow:
  def __init__(self, n):
    self.n, self.adj = n, [[] for _ in range(n)]
    self.to, self.cap, self.cost = [], [], []

  def add_edge(self, u, v, cap, cost):
    i = len(self.to)
    self.to += [v, u]
    self.cap += [cap, 0]
    self.cost += [cost, -cost]
    self.adj[u].append(i)
    self.adj[v].append(i + 1)
    return i

  def flow(self, i):
    return self.cap[i ^ 1]

  def max_flow(self, s, t, limit=float("inf")):
    from heapq import heappop, heappush

    n, adj, to = self.n, self.adj, self.to
    cap, cost = self.cap, self.cost
    pot = self.pot = [0] * n
    for _ in range(n):  # Bellman-Ford: valid potentials
      changed = False
      for u in range(n):
        for e in adj[u]:
          if cap[e] and pot[u] + cost[e] < pot[to[e]]:
            pot[to[e]] = pot[u] + cost[e]
            changed = True
      if not changed:
        break
    INF, f, c, pe = 1 << 62, 0, 0, [0] * n
    while f < limit:
      d = [INF] * n
      d[s], pq = 0, [s]
      while pq:
        du, u = divmod(heappop(pq), n)
        if du > d[u]:
          continue
        if u == t:
          break
        base = du + pot[u]
        for e in adj[u]:
          if cap[e]:
            v = to[e]
            nd = base + cost[e] - pot[v]
            if nd < d[v]:
              d[v], pe[v] = nd, e
              heappush(pq, nd * n + v)
      dt = d[t]
      if dt == INF:
        break
      for v in range(n):
        pot[v] += d[v] if d[v] < dt else dt
      push, v = limit - f, t
      while v != s:
        push = min(push, cap[pe[v]])
        v = to[pe[v] ^ 1]
      v = t
      while v != s:
        cap[pe[v]] -= push
        cap[pe[v] ^ 1] += push
        v = to[pe[v] ^ 1]
      f += push
      c += push * (pot[t] - pot[s])
    return f, c
