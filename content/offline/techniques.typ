#import "/typst/notes.typ": *

Offline means every operation is read before any answer is due, so
the answers may be computed in any order. Choose the order that
makes the structure static, insert-only or monotone. It does not
apply when queries are decoded with the previous answer (forced
online) or the problem is interactive, unless the answers are
predictable (Fake online below). Keep each query's index and
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
  insert before they ask. A contribution that is a polynomial in
  the query coordinate, $sum_k c_k x^k$, needs one Fenwick tree per
  power $k$ (the coefficients $w, w l, w l^2$ of each update).]
#fact[All left ends][Sweep `r` and keep $f(l, r)$ for every
  `l <= r` in a lazy segment tree over `l`. Appending `a[r]`
  changes $f(dot, r)$ on a few contiguous ranges of `l`, found with
  a monotone stack: `a[r]` is the new maximum exactly for `l` after
  the previous element `>= a[r]`, and each popped entry is a range
  whose old maximum is replaced. A query `(l, r)` reads `l` at time
  `r`; "how many good `l`" is a min and its count; sums over all
  subarrays inside `[L, R]` need historic sums (Linear tags in
  _Designing a segment tree node_).]
#fact[Rectangles in the $(l, r)$ plane][To count pairs or triples
  over all subarrays of `[L, R]`, describe the valid $(l, r)$ of
  each element as $O(1)$ rectangles $[a, b] times [c, d]$, $O(n)$ in
  total. Sweep $L$ downwards: a rectangle adds $(b + 1 - L) g(R)$
  while $a <= L <= b$ and $(b - a + 1) g(R)$ once $L < a$, with
  $g(R) = |[c, d] inter [1, R]|$. Range adds on $[c, d]$ into
  Fenwick trees for the count, $sum (b + 1)$ and $sum (b - a + 1)$
  give the answer as a combination linear in $L$.]
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
#fact[Fake online][Queries decoded with the previous answer, but
  the answers are monotone in time (connected until some deletion,
  then never again): decode as if every answer were the first
  kind, find the switch offline (reverse time), and print the
  other kind from there on. Queries decoded wrongly after the
  switch cannot change any answer.]
#fact[Final positions]["Insert after the `k`-th element", mixed
  with queries: walk the insertions backwards over $q$ free slots
  (a Fenwick tree of ones); each takes the $(k + 1)$-th free slot
  (`lower_bound`) and removes it. That slot is its final index, so
  the sequence becomes static: replay forwards, switching slots
  on, and find the `k`-th active slot with the same descent.]
#fact[Thresholds][$max - min$ of a set is the number of integers
  $t$ with some value $<= t$ and some value $> t$, and
  $w = sum_(t >= 0) [w > t]$. So costs built from $max - min$ or
  $|x - y|$ split into one 0/1 problem per threshold, on
  $[a_i > t]$; thresholds between consecutive relevant values
  share a problem, weighted by the gap. An adjacent pair is a 1
  followed by a 0 exactly for $t$ in $[a_(i + 1), a_i)$: one range
  add per pair on a segment tree over thresholds, two pairs per
  point update. MST weight by thresholds: _Graph theory facts_.]
#fact[Nested thresholds][To minimise such a cost over free cells,
  solve each 0/1 problem; the optima add up to the true optimum if
  they can be chosen nested (decreasing in $t$). That holds when
  every threshold has the same cost $f$, submodular:
  $f(c) + f(d) >= f(max(c, d)) + f(min(c, d))$ elementwise, and
  only the fixed cells, $[a_i > t]$, change with $t$. Divide and
  conquer over thresholds: solve the middle one; its free cells
  set to 1 stay 1 below it, those set to 0 stay 0 above, so each
  free cell recurses into one half only.]

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
#fact[All but one][For every $i$, a structure of all items except
  $i$ (or except group $i$), when items can be inserted but not
  removed (knapsack bitsets, products of polynomials): recurse on
  `[l, r)` carrying the structure of everything outside it; insert
  `[m, r)` into a copy before recursing left and `[l, m)` before
  recursing right. Each item is inserted $O(log n)$ times. For
  prefix times suffix products, carry only the coefficients the
  leaves below can still need.]

=== Other offline tools
- Mo's algorithm (`offline/mo`): range queries whose answer does
  not merge from halves but changes in $O(1)$ when an end moves.
- Mo with deletions only, when removing is cheap and adding is
  not (unlink from a sorted linked list; relinking in reverse
  order undoes it): per block of left ends, build the window
  `[block start, n)` once, take its queries by decreasing `r` and
  delete from the right for good; per query delete
  `[block start, l)`, answer, and undo. $O((n + q) sqrt(n))$.
- Blocks of $sqrt(q)$ operations: rebuild a static structure per
  block and handle the block's own updates by brute force
  (Periodic rebuild in _Sqrt decomposition and small tricks_).
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
- "Sum over all subarrays of `[l, r]`", "count good subarrays":
  all left ends, or rectangles in the $(l, r)$ plane.
- "For each item, the answer without it": all but one.

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
