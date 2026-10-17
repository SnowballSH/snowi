"""
Title: Offline dynamic connectivity
Description: Edges are added and removed in time order, and
  `query()` marks a moment to ask about the graph. `solve` then
  calls `answer(i)` for the queries `i = 0, 1, ...` in order, each
  time with the rollback DSU `d` holding exactly the edges present
  at query `i`.
Precondition: Vertices `0..n-1`. `remove(u, v)` only while that
  edge is present; parallel edges and loops are fine. Call `solve`
  once, after all operations.
Time: $O(m log q log n)$ for $m$ added edges and $q$ queries. In
  PyPy, $2 dot 10^5$ operations on $10^5$ vertices take about 1 s.
Usage:
  dc = DynConn(n)
  dc.add(u, v); dc.remove(u, v)  # in time order
  i = dc.query()  # ask about the graph at this moment
  def answer(i):
    same[i] = dc.d.find(a[i]) == dc.d.find(b[i])
    comps[i] = n - dc.d.time()
  dc.solve(answer)
Tweak: Bipartiteness, sums or counts per component: give the DSU
  that state, record it in
  `join` and restore it in `rollback`. Any structure with "insert"
  and "undo the last inserts" works instead of the DSU (XOR basis,
  knapsack array: save a copy per node of the segment tree over
  query times). Prune: once a node's edges decide all queries
  below it (an odd cycle), answer them there and do not descend.
Example:
  >>> dc, out = DynConn(4), []
  >>> dc.add(0, 1); dc.add(1, 2); dc.query()
  0
  >>> dc.remove(0, 1); dc.add(2, 3); dc.query()
  1
  >>> def answer(i):  # components, is 0 connected to 2?
  ...   d = dc.d
  ...   out.append((4 - d.time(), d.find(0) == d.find(2)))
  >>> dc.solve(answer); out
  [(2, True), (2, False)]
Notes: Think of it when edges (or items) appear and disappear,
  everything is known in advance, and the structure can insert and
  undo but not delete: connectivity, bipartiteness (DSU with
  parity), sums or counts per component, XOR basis, knapsack. An
  edge is alive on the half-open range from its addition to its
  removal: a query between them sees it, one after the removal
  does not. Only deletions: reverse time, so that they become
  insertions into a plain DSU. Forced online (queries decoded with
  the last answer): this does not apply. Memory is $O(m log q)$
  stored edges.
Uses: ds/dsu-rollback
Source: Segment tree over time; kactl RollbackUF (CC0); UCup
  S4-S5 Training Packet 14.
Status: stress-tested against brute force
"""


class DynConn:
  def __init__(self, n):
    self.n, self.q, self.d = n, 0, RollbackDSU(n)
    self.open = {}  # u * n + v -> start times of present copies
    self.ev = []  # flat l, r, u, v: alive for [l, r)

  def add(self, u, v):
    key = min(u, v) * self.n + max(u, v)
    self.open.setdefault(key, []).append(self.q)

  def remove(self, u, v):
    l = self.open[min(u, v) * self.n + max(u, v)].pop()
    self.ev += (l, self.q, u, v)

  def query(self):
    self.q += 1
    return self.q - 1

  def solve(self, answer):
    n, q, ev, d = self.n, self.q, self.ev, self.d
    for key, starts in self.open.items():
      for l in starts:
        ev += (l, q, key // n, key % n)
    S = 1
    while S < q:
      S *= 2
    at = [[] for _ in range(2 * S)]
    for k in range(0, len(ev), 4):
      l, r = ev[k] + S, ev[k + 1] + S
      u, v = ev[k + 2], ev[k + 3]
      while l < r:
        if l & 1:
          at[l] += (u, v)
          l += 1
        if r & 1:
          r -= 1
          at[r] += (u, v)
        l >>= 1
        r >>= 1
    join, stack = d.join, [1]
    while stack:
      x = stack.pop()
      if x < 0:  # leaving a node: undo its joins
        d.rollback(~x)
        continue
      t, e = d.time(), at[x]
      for k in range(0, len(e), 2):
        join(e[k], e[k + 1])
      if x < S:
        stack += (~t, 2 * x + 1, 2 * x)
      else:
        if x - S < q:
          answer(x - S)
        d.rollback(t)
