#import "/typst/notes.typ": *

A matroid is a "can these coexist?" rule for which greedy by weight
is always optimal. Recognise one, and Kruskal's proof covers your
greedy; with two such rules at once, use matroid intersection.

=== Definitions
#fact[Matroid][A finite ground set $E$ and a family of
  *independent* sets such that (1) $emptyset$ is independent;
  (2) subsets of independent sets are independent; (3) exchange:
  if $A$, $B$ are independent and $|A| > |B|$, some
  $x in A without B$ makes $B + x$ independent.]
- An *independence oracle* answers "is $S$ independent?"; the
  algorithms below only call it.
- A *basis* of $S$ is a maximal independent subset; all have the
  same size, the *rank* $r(S)$. Rank is submodular:
  $r(A union B) + r(A inter B) <= r(A) + r(B)$.
- A *circuit* is a minimal dependent set (a cycle of a graph). If
  $B$ is independent and $B + e$ is not, $B + e$ contains exactly
  one circuit $C(B, e)$, and $B + e - f$ is independent iff
  $f in C(B, e)$.
- Restrictions, truncations ("and at most $k$ elements"), direct
  sums and duals are matroids. The dual's bases are the complements
  of bases, $r^*(S) = |S| - r(E) + r(E without S)$.

=== Greedy
#fact[Rado-Edmonds][Scan the elements by decreasing weight and keep
  each one that leaves the kept set independent: the result is a
  maximum-weight basis for every weight function. Stop at the first
  weight $<= 0$ for a maximum-weight independent set. Conversely,
  if (1) and (2) hold and this greedy is optimal for all weights,
  then (3) holds.]
- The first $k$ picks form a best independent set of size $k$, so
  the best weight is concave in $k$.
- A basis $B$ is maximum iff no single swap improves it: for each
  $e in.not B$, every $f in C(B, e)$ has $w(f) >= w(e)$.
- Insert $e$: add it, and if that closes a circuit drop the
  circuit's lightest element. Delete $f in B$: add the heaviest
  element outside that fits.
- Latest basis: weight = position. Sweep $r$, keeping the
  maximum-weight basis $B_r$ of the prefix. Then $B_r$ restricted
  to $[l, r]$ is a basis of $[l, r]$: its rank is the number of
  elements of $B_r$ at positions $>= l$ (queries sorted by $r$, a
  Fenwick tree over positions). XOR vectors: at each pivot keep the
  newer vector and reduce the older one (`linalg/xor-basis`);
  graphs: an edge closing a cycle evicts the oldest edge on it
  (`tree/link-cut`).

=== Examples
#table(
  columns: 3,
  [*Matroid*], [*Independent sets*], [*Test*],
  [uniform], [at most $k$ elements], [count],
  [partition], [at most $c_j$ of colour $j$], [counts],
  [graphic], [forests], [`ds/dsu`],
  [linear], [independent vectors], [`linalg/xor-basis`, `linalg/gauss-mod`],
  [transversal], [vertices of one side matchable together], [`flow/hopcroft-karp`],
  [scheduling], [unit jobs that all meet their deadlines], [for all $D$: at most $D$ jobs due by $D$],
  [cographic], [edges whose removal keeps every component connected], [`ds/dsu` on the rest],
)
- Scheduling is transversal (job $i$ fits slots $1, dots, d_i$), and
  for prefix neighbourhoods Hall's condition is the prefix test.
  Most profit: by deadline, push profits into a min-heap and pop
  the cheapest whenever it holds more jobs than the deadline.
- Lower bounds "at least $c_j$ chosen among the first $p_j$" with
  $N$ chosen in all (caps on the unchosen prefixes: a matroid).
  Heaviest choice: sweep, push each weight into a max-heap, and at
  $p_j$ pop maxima until $c_j$ are chosen. The openings of a
  bracket sequence of length $2n$ are exactly the $n$-sets with at
  least $i$ among the first $2i - 1$ positions.
- Not matroids: matchings of a general graph as edge sets,
  knapsack, paths, independent sets of a graph ("no two adjacent").
  Check (3) on two tiny sets before trusting a greedy.

=== Matroid intersection
Find the largest (or heaviest) $S$ independent in two matroids
$M_1$ and $M_2$ on the same ground set. $S$ grows by one per round:
- Exchange graph: for $x in S$ and $y in.not S$, an arc $x -> y$ if
  $S - x + y$ is independent in $M_1$, and $y -> x$ if it is in
  $M_2$. Sources: $y in.not S$ with $S + y$ independent in $M_1$;
  sinks: $y$ with $S + y$ independent in $M_2$.
- Take a *shortest* path from a source to a sink (BFS) and flip
  membership along it. A longer path can break independence. No
  path means $S$ is maximum.
- At most $r + 1$ rounds of $O(r n)$ oracle calls: $O(r^2 n)$ calls
  for $n = |E|$ and answer size $r$.
- Edmonds: $max |S| = min_(A subset.eq E) (r_1 (A) + r_2 (E without A))$.
- Weighted: cost $w(x)$ on $x in S$ and $-w(y)$ on $y in.not S$;
  take the cheapest path, the fewest arcs among those
  (Bellman-Ford; there are no negative cycles). After $k$ rounds
  $S$ is a heaviest common independent set of size $k$, and that
  weight is concave in $k$.
- Fast exchange tests. Graphic: $S - x + y$ is a forest iff $S + y$
  is, or $x$ lies on the path in $S$ between the ends of $y$.
  Partition: $y$'s colour has room, or $x$ has $y$'s colour.
  Linear: $S + y$ is independent, or $x$ has a nonzero coefficient
  when $y$ is written in terms of $S$. Transversal (the right
  vertices matched by a matching $M$ of $S$): search from $y$ to
  its left neighbours, from each matched one to its partner in
  $M$, and on; $S + y$ is independent iff the search reaches a
  free left vertex, and otherwise $S - x + y$ is iff it reaches
  $x$. Scheduling: the slack tree in _Designing a segment tree
  node_.
- Pairs: bipartite matching (one partition per side); rainbow or
  colour-capped spanning forest (graphic, colours); spanning
  arborescence (graphic, in-degree $<= 1$ and 0 at the root); a
  common forest of two graphs on the same edges; edges whose
  removal keeps two graphs connected (two cographic). Three
  matroids are NP-hard (Hamiltonian path).
#lang-code(
  cpp: "// Largest set independent in two matroids. in1(S), in2(S) test
// a membership vector. O(r^2 n) oracle calls, r = answer size.
template<class F, class G>
vector<bool> intersect(int n, F in1, G in2) {
  vector<bool> S(n);
  auto test = [&](auto& f, int add, int rem) { // S + add - rem
    S[add] = 1;
    if (rem >= 0) S[rem] = 0;
    bool ok = f(S);
    S[add] = 0;
    if (rem >= 0) S[rem] = 1;
    return ok;
  };
  for (;;) {
    vi q, par(n, -2);
    rep(y, 0, n) {
      if (!S[y] && test(in1, y, -1))
        q.push_back(y), par[y] = -1;
    }
    int t = -1;
    for (int i = 0; i < sz(q) && t < 0; i++) { // BFS
      int u = q[i];
      if (!S[u] && test(in2, u, -1)) t = u;
      rep(v, 0, n) {
        if (t >= 0 || par[v] != -2 || S[v] == S[u]) continue;
        if (S[u] ? test(in1, v, u) : test(in2, u, v))
          q.push_back(v), par[v] = u;
      }
    }
    if (t < 0) return S;
    for (; t >= 0; t = par[t]) S[t] = !S[t];
  }
}",
  py: "def intersect(n, in1, in2):
  # Largest set independent in two matroids; in1(S), in2(S)
  # test a 0/1 list. O(r^2 n) oracle calls, r = answer size.
  S = [0] * n

  def test(f, add, rem):  # is S + add - rem independent?
    S[add] = 1
    if rem >= 0:
      S[rem] = 0
    ok = f(S)
    S[add] = 0
    if rem >= 0:
      S[rem] = 1
    return ok

  while True:
    q = [y for y in range(n) if not S[y] and test(in1, y, -1)]
    par = [-2] * n
    for y in q:
      par[y] = -1
    t = -1
    for u in q:  # BFS; q grows during the loop
      if not S[u] and test(in2, u, -1):
        t = u
        break
      for v in range(n):
        if par[v] == -2 and S[v] != S[u]:
          if test(in1, v, u) if S[u] else test(in2, u, v):
            par[v] = u
            q.append(v)
    if t < 0:
      return S
    while t >= 0:
      S[t] ^= 1
      t = par[t]",
)

=== Matroid union
For matroids $M_1, dots, M_k$ the unions $I_1 union dots union I_k$
of disjoint $I_j$ independent in $M_j$ form a matroid with rank
$r(S) = min_(T subset.eq S) (|S without T| + sum_j r_j (T))$.
Compute it as an intersection on the copies $(e, j)$: the direct sum
of the $M_j$ against "each $e$ used at most once".
- Nash-Williams, Tutte: a graph has $k$ edge-disjoint spanning trees
  iff every partition of $V$ into $p$ parts has at least $k(p - 1)$
  edges between parts. Its edges split into $k$ forests iff every
  vertex set $X$ spans at most $k(|X| - 1)$ edges.
- Two linear matroids with vectors $u_i$ and $v_i$ over a large
  prime field: with random nonzero $a_i, b_i$, the rank of the
  vectors $(a_i u_i, b_i v_i)$ is, with high probability, the most
  elements in disjoint $I_1 in M_1$, $I_2 in M_2$. Then plain
  elimination and the latest basis apply.

=== Recognition
- "Choose a subset of largest weight with no cycle / independent
  vectors / all matchable / all on time / at most $c_j$ per
  colour": sort and run the greedy.
- A heap greedy you cannot prove ("keep the best, throw out the
  cheapest when overfull"): check the exchange property (3).
- Many queries "rank, largest XOR or forest size of the elements
  in $[l, r]$": the latest basis.
- Two such rules at once: intersection. "$k$ disjoint spanning
  trees" or "two disjoint bases": union.
- "Spanning tree with exactly $k$ white edges": its least weight is
  convex in $k$, so use `dp/wqs`: add a price $lambda$ to every
  white edge and binary search $lambda$.
- Pitfalls: a maximum-weight basis (keep going past negative
  weights) differs from a maximum-weight independent set; loops
  (self-loop edges, zero vectors) are in no basis; break ties by
  index when a canonical basis is needed.
