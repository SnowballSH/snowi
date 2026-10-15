#import "/typst/notes.typ": *

Offline means every operation is read before any answer is due, so
the answers may be computed in any order. Choose the order that
makes the structure static, insert-only or monotone. It does not
apply when queries are decoded with the previous answer (forced
online) or the problem is interactive. Keep each query's index and
print in input order.

=== Reordering
#fact[Sort the queries][Answer in an order in which the structure
  only grows. "Using only edges of weight $<= w$": sort edges and
  queries by weight, join into a DSU, answer when the weight is
  reached. "Count $a_j <= x$ in `[l, r)`": sort elements and
  queries by value, add positions to a Fenwick tree (`ds/fenwick`).
  A tree query "only vertices with $t_v <= x$" works the same way
  with Euler-tour ranges.]
#fact[Sweep line][Make one coordinate time: events at each $x$
  (insert, delete, ask) in order, with a Fenwick or segment tree
  over $y$. Static rectangle counts (four prefix corners with signs
  $plus.minus 1$), the area of a union of rectangles (min and its
  count), pairs $i < j$ with a condition on $(a_i, a_j)$: sweep
  $j$, query the inserted $i$. Fix whether events at equal $x$
  insert before they ask.]
#fact[Last occurrence][Distinct values in `a[l, r)`: sweep `r` and
  keep a mark only at the last position of each value seen; the
  answer is the number of marks in `[l, r)`. With `prev[i]`, the
  previous position of the same value, "first occurrence in the
  range" becomes `prev[i] < l`, a 2D dominance condition.]
#lang-code(
  cpp: "vector<vi> at(n + 1); // qs[i] = {l, r}, values in [0, V)
rep(i, 0, sz(qs)) at[qs[i].second].push_back(i);
Fenwick fw(n);
vi last(V, -1), res(sz(qs));
rep(r, 0, n + 1) {
  for (int i : at[r])
    res[i] = fw.query(r) - fw.query(qs[i].first);
  if (r == n) break;
  if (last[a[r]] >= 0) fw.add(last[a[r]], -1);
  fw.add(last[a[r]] = r, 1);
}",
  py: "at = [[] for _ in range(n + 1)]  # qs[i] = (l, r), a[i] < V
for i, (l, r) in enumerate(qs):
  at[r].append(i)
fw, last, res = Fenwick(n), [-1] * V, [0] * len(qs)
for r in range(n + 1):
  for i in at[r]:
    res[i] = fw.query(r) - fw.query(qs[i][0])
  if r < n:
    if last[a[r]] >= 0:
      fw.add(last[a[r]], -1)
    last[a[r]] = r
    fw.add(r, 1)",
)
#fact[Reverse time][Only deletions (cut edges, remove vertices,
  block cells): apply all of them, then walk backwards so that
  each deletion becomes an insertion into a DSU or a growing set,
  and reverse the answers. "Components after each removal" is the
  standard case; an element that is never deleted is inserted at
  the start.]

=== Divide and conquer
#fact[CDQ][Updates and queries interleaved in time, or a 3D partial
  order: split the operations in half; the left half's updates
  reach the right half's queries by one static sweep, so each
  element takes part in $O(log n)$ sweeps: $O(n log^2 n)$
  (`offline/cdq`). For a DP over the order, recurse left, then the
  cross step, then right.]
#fact[Segment tree over time][Objects alive on time intervals and a
  structure that can insert and undo but not delete (DSU, XOR
  basis, knapsack array): store each interval in $O(log q)$ nodes
  and run a DFS that inserts on entry and rolls back on exit
  (`offline/dynamic-connectivity`).]
#fact[Parallel binary search][Each query wants the first time (or
  the smallest value) at which its monotone condition holds, and a
  check needs the state after $t$ updates: run all binary searches
  together, one sweep over the updates per round, $O(log T)$
  rounds (`offline/parallel-binary-search`).]
#fact[Queries crossing a midpoint][Static range queries whose
  answer does not merge cheaply (best knapsack value on `[l, r)`,
  a non-invertible operation): recurse on `[lo, hi)` with middle
  `m`; queries inside one half go down, and a query crossing `m`
  combines a table for its suffix of `[lo, m)` with a table for its
  prefix of `[m, hi)`. Each level extends $O(n)$ tables by one
  element: $O(n log n)$ extensions plus one combination per
  query.]
#fact[Requests per segment-tree node][A range query splits into
  $O(log n)$ canonical nodes; when a node's part of the answer
  depends on values from the query (the maximum left of the node),
  attach the request to the node and answer all of a node's
  requests with one sorted sweep: $O((n + q) log^2 n)$.]

=== Other offline tools
- Mo's algorithm (`offline/mo`): range queries whose answer does
  not merge from halves but changes in $O(1)$ when an end moves.
- Blocks of $sqrt(q)$ operations: rebuild a static structure per
  block and handle the block's own updates by brute force
  (`ds2/sqrt`, periodic rebuild).
- Offline LCA, tree queries sorted by depth or by DFS order, and
  small-to-large merging are offline too: all queries of a vertex
  are answered when its subtree is complete.

#full[
=== Cues
- "Answer $q$ queries" with all of them in the input and no
  decoding: think offline first.
- "Removed one by one", "after each deletion": reverse time.
- "Edges appear and disappear", "toggle an edge": segment tree over
  time.
- "For each person, the first day when ...": parallel binary
  search. "K-th smallest in a range with updates": parallel binary
  search over values.
- "Count pairs $j < i$ with $x_j <= x_i$ and $y_j <= y_i$", "points
  inserted, rectangle asked": sweep, or CDQ when the inserts and
  asks interleave.

=== Pitfalls
- Equal keys: decide whether an element with the same key as a
  query counts, then sort events so that the rule holds.
- Time intervals are half-open: an edge removed at time $t$ is
  absent at time $t$.
- Undo exactly what was inserted; clearing a whole structure in
  every recursive call costs $O(n)$ each time.
#py[
- PyPy: recursion of depth $O(log n)$ is fine, a DFS over a tree is
  not: use an explicit stack.
]
]
