#import "/typst/notes.typ": *

Statements rarely say "flow". Build a network, run `flow/dinic`,
and read the answer from the flow value or from the minimum cut:
the source side is what `s` still reaches in the residual graph.

=== Max flow and min cut
#fact[Max-flow min-cut][The maximum flow equals the minimum
  capacity of a cut: a vertex set $S$ with $s in S$, $t in.not S$,
  costing the capacities of edges from $S$ to the rest. Integer
  capacities give an integer maximum flow.]
#fact[Supply and demand][Supplies $b_v > 0$ and demands $b_v < 0$
  (summing to 0) can all be met iff every vertex set $X$ has net
  supply
  $sum_(v in X) b_v$ at most the capacity of the edges leaving
  $X$. When only a few $X$ can be tight (prefixes, subtrees, sides
  of a bridge, subsets of few types) check those directly.]
#fact[Disjoint paths][The most edge-disjoint $s$-$t$ paths, and
  the fewest edges whose removal separates $t$ from $s$, both
  equal the max flow with unit capacities (Menger). For
  vertex-disjoint paths also split every vertex except $s$ and $t$
  with capacity 1 (or run from $s_"out"$ to $t_"in"$).]

=== Building blocks
- Vertex capacity $c$: split $v$ into $v_"in" -> v_"out"$ with
  capacity $c$; edges enter $v_"in"$ and leave $v_"out"$.
- Undirected edge: one `add_edge(u, v, c, c)`.
- Several sources or sinks: a super source $S -> s_i$ (capacity
  the supply of $s_i$) and $t_j -> T$.
- Infinite edge $u -> v$: never cut, so "$u$ on the source side
  forces $v$ there too". Make it larger than every finite cut.
- "Each item usable $k$ times", "at most $k$ per group": a
  capacity-$k$ edge, through a group vertex if needed.
- Huge counts of few kinds: one vertex per kind with summed
  capacities; binary search an answer $K$ outside the flow.

=== Lower bounds
Edge $u -> v$ must carry between $l$ and $c$: add it with capacity
$c - l$ and record `ex[v] += l`, `ex[u] -= l`. Add $S -> v$ with
capacity `ex[v]` when it is positive and $v -> T$ with `-ex[v]`
when negative; for an $s$-$t$ flow also keep the id
`ts = add_edge(t, s, INF)` (`rev_cap` INF too if the value may be
negative). A feasible flow exists iff `max_flow(S, T)` saturates
every edge out of $S$. Then `flow(ts)` is a feasible value: read
it, then zero both capacities of `ts`. For the largest value add
`max_flow(s, t)` to it; for the smallest subtract
`max_flow(t, s)`. An edge carries $l$ plus its flow.

=== Choices as cuts
#fact[Project selection][Item $i$ has weight $w_i$ (profit if
  $w_i > 0$, cost if $w_i < 0$); taking $i$ requires taking $j$.
  Edges $s -> i$ with $w_i$ for profits, $i -> t$ with $-w_i$ for
  costs, $i -> j$ with $infinity$. Best total
  $= sum_(w_i > 0) w_i - "min cut"$; take the source side.]
#fact[Two labels][Each $x_v in {0, 1}$ (0 = source side). Cost $a$
  if $x_v = 1$: edge $s -> v$ with $a$. Cost $b$ if $x_v = 0$:
  $v -> t$ with $b$. Penalty $p >= 0$ if $x_u = 0$ and $x_v = 1$:
  $u -> v$ with $p$ (both directions for "differ"). Only
  penalties for disagreeing work; on a bipartite grid, flip the
  labels of one color class to turn "equal" into "differ".]

=== Bipartite matching (`flow/hopcroft-karp`)
#fact[Konig][Maximum matching $nu$ = minimum vertex cover. The
  rest of a cover is independent, so the maximum independent set
  is $|V| - nu$; so is the minimum edge cover (no isolated
  vertices).]
#fact[Hall][The left side $L$ can be fully matched iff
  $|N(X)| >= |X|$ for all $X subset.eq L$. In general
  $nu = |L| - max_X (|X| - |N(X)|)$. A $k$-regular bipartite
  multigraph has a perfect matching, so it splits into $k$ of
  them.]
#fact[Path cover][Vertex-disjoint paths covering a DAG: $n - nu$,
  matching out-copies to in-copies along edges. If paths may
  share vertices, match on the transitive closure.]
#fact[Dilworth][In a partial order, the fewest chains covering
  everything = the largest antichain = $n - nu$ on the strict
  order relation (transitive closure).]
- Costs on pairs: `flow/hungarian` (dense, $n <= 500$) or
  `flow/min-cost-flow`.

=== Min-cost flow (`flow/min-cost-flow`)
- $k$ disjoint paths of least total length: unit capacities,
  `max_flow(s, t, k)`; split vertices for vertex-disjoint.
- Maximum-weight matching of any size: costs $-w$, and stop
  when the next path costs $>= 0$ (see the Tweak: test it right
  after the potential update).
- Convex cost per unit: parallel unit edges with increasing cost.

=== When it is not a flow
- Plain implications between booleans with no costs, "not both",
  "at least one of two": `graph/two-sat`.
- $n$ up to $10^6$: use the cut condition to prove a criterion,
  then check it greedily. Transport on a path or tree: each edge
  carries the imbalance $b_e$ of the side it cuts off, so the cost
  is $sum_e c_e |b_e|$.

#full[
=== Worked examples
*Project selection.* P1 (+5) needs T1 (-3) and T2 (-4); P2 (+6)
needs T2; P3 (+2) needs T3 (-3). Profits total 13. The minimum
cut is 9: $s -> "P3"$ (2), $"T1" -> t$ (3), $"T2" -> t$ (4). Best
total $13 - 9 = 4$, from the source side {P1, P2, T1, T2}.

*Lower bounds.* Feasibility, then the largest and smallest value,
with `S = n`, `T = n + 1`:
#lang-code(
  cpp: "Dinic d(n + 2);
vector<ll> ex(n);
for (auto [u, v, l, c] : edges)
  d.add_edge(u, v, c - l), ex[v] += l, ex[u] -= l;
int ts = d.add_edge(t, s, 1e18); // only for an s-t flow
ll need = 0;
rep(v, 0, n) {
  if (ex[v] > 0) d.add_edge(n, v, ex[v]), need += ex[v];
  if (ex[v] < 0) d.add_edge(v, n + 1, -ex[v]);
}
bool ok = d.max_flow(n, n + 1) == need;
ll F = d.flow(ts); // a feasible value, if ok
auto [a, i] = d.pos[ts];
Dinic::E& x = d.g[a][i];
x.cap = d.g[x.to][x.rev].cap = 0; // remove t -> s
ll hi = F + d.max_flow(s, t), lo = hi - d.max_flow(t, s);",
  py: "d, ex = Dinic(n + 2), [0] * n
for u, v, l, c in edges:
  d.add_edge(u, v, c - l)
  ex[v] += l
  ex[u] -= l
ts = d.add_edge(t, s, 1 << 60)  # only for an s-t flow
for v in range(n):
  if ex[v] > 0:
    d.add_edge(n, v, ex[v])
  elif ex[v] < 0:
    d.add_edge(v, n + 1, -ex[v])
need = sum(x for x in ex if x > 0)
ok = d.max_flow(n, n + 1) == need
F = d.flow(ts)  # a feasible value, if ok
d.cap[ts] = d.cap[ts ^ 1] = 0  # remove t -> s
hi = F + d.max_flow(s, t)
lo = hi - d.max_flow(t, s)",
)

*Konig.* Left 1, 2, 3; right x, y; edges 1x, 1y, 2x, 3x. Hall
fails for $X = {2, 3}$ with $N(X) = {x}$, so $nu = 3 - 1 = 2$. A
minimum cover is {1, x}: every edge touches 1 or x.

*Vertex cut.* The fewest vertices (not $s$, $t$) whose removal
disconnects $t$ from $s$: split each vertex with capacity 1 (and
$infinity$ for $s$ and $t$), original edges with $infinity$; the
answer is the max flow, the cut vertices are those with $v_"in"$
on the source side and $v_"out"$ not.

=== Recognition
- "Assign each A to a B, B has capacity", "schedule jobs into
  slots": bipartite flow, or Hungarian with costs.
- "Fewest cells to block", "fewest roads to destroy": min cut,
  splitting vertices when cells (not borders) are blocked.
- "Choose a set, choosing $x$ forces $y$, maximize profit":
  closure. "Edges minus vertices" (densest subgraph): each edge
  is a profit item requiring its two endpoints.
- "Maximum set with no two conflicting" where conflicts are
  bipartite (cells of a chessboard, odd and even sums):
  $|V| - nu$.
- "Is it possible?" with a flow-like process and huge limits:
  write the cut capacity as a formula and find which sets can be
  tight.

=== Greedy instead of flow
- Points to intervals: sweep the points left to right; give each
  the free interval containing it that ends first (a heap).
- Neighborhoods that are prefixes $[0, r_j)$: sort the $r_j$; the
  $k$-th smallest must be $>= k$.
- Matching on a tree: match a leaf with its parent, bottom-up.
- Planar $s$-$t$ min cut with $s$ and $t$ on the outer face:
  split that face by an $s$-$t$ edge; the cut is a shortest path
  in the dual between the two halves (grid with $s$ left and $t$
  right: top border to bottom border). Otherwise it is a shortest
  dual cycle separating $s$ from $t$.

=== Pitfalls
- An "infinite" capacity must beat every finite cut, yet a few
  of them added together must #cpp[not overflow `ll`]#py[stay
  below $2^63$, or PyPy falls back to slow big ints].
- Read a cut from the residual reachability after the final
  flow, never from "saturated edges".
- Is $u$ in every maximum matching? Exactly when deleting $u$
  lowers $nu$.
- Capacities on vertices and edges at once: split the vertices
  and keep the edge capacities on the $v_"out" -> w_"in"$ edges.
]
