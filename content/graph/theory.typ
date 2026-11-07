#import "/typst/notes.typ": *

$n$ vertices, $m$ edges, $Delta$ the maximum degree. $nu$, $alpha$,
$tau$, $rho$: maximum matching, maximum independent set, minimum
vertex cover, minimum edge cover.

=== Euler tours and parity
#full[#fact[Euler][Ignore isolated vertices; the rest must be
  connected. Undirected: a circuit iff every degree is even, a path
  iff exactly 0 or 2 degrees are odd (it joins the odd ones).
  Directed: a circuit iff $"in"(v) = "out"(v)$ everywhere, a path
  iff that holds except that one vertex may have $"out" = "in" + 1$
  (the start) and one $"in" = "out" + 1$ (the end).]]
#fact[Trails][Connected, $k > 0$ odd vertices: the edges split into
  $k\/2$ trails and no fewer (`graph/euler-path` after pairing the
  odd vertices by dummy edges). A connected digraph needs
  $max(1, sum_v max(0, "out"(v) - "in"(v)))$ trails. Chinese
  postman (closed walk using every edge): add a minimum-weight
  perfect matching of the odd vertices under shortest-path
  distances (bitmask DP for up to about 20 of them).]
#fact[Parity fixing][Connected graph, $T$ a vertex set of even
  size: in a spanning tree keep the edge above $v$ iff the subtree
  of $v$ holds an odd number of vertices of $T$. Exactly the
  vertices of $T$ get odd degree. Likewise, bottom-up
  multiplicities $0 <= m_e < k$ on the tree edges reach any
  degrees mod $k$ at all vertices but the root. The root's residue
  is then forced: both colour classes of the 2-coloured tree have
  degree sum $sum_e m_e$ (for $k = 2$: the total must be even).]

=== Matchings, covers and orders
#fact[Gallai, Kőnig][$alpha + tau = n$; $nu + rho = n$ without
  isolated vertices (extend a maximum matching greedily).
  Bipartite: $nu = tau$.]
#fact[Hall][Bipartite with sides $L, R$:
  $nu = |L| - max_(X subset.eq L) (|X| - |N(X)|)$. A $k$-regular
  bipartite multigraph (a nonnegative integer matrix with all row
  and column sums $k$) splits into $k$ perfect matchings
  (permutation matrices).]
#fact[Tutte][A perfect matching exists iff removing any set $U$
  leaves at most $|U|$ odd components. Tutte-Berge:
  $nu = min_U (n + |U| - "odd"(G - U))\/2$.]
#fact[Dilworth, Mirsky][In a partial order, the largest antichain
  equals the fewest chains covering everything ($n - nu$ on the
  strict relation), and the longest chain equals the fewest
  antichains covering everything (layer by height). On sequences:
  the fewest non-increasing subsequences covering a sequence is
  the length of its longest strictly increasing subsequence. A DAG
  is covered by $n - nu$ vertex-disjoint paths, matching out-copies
  to in-copies along edges; paths that may share vertices: the same
  on the transitive closure.]
#fact[Stable matching][Gale-Shapley: each free proposer proposes to
  the next on its list; a receiver keeps its best offer so far.
  $O(n^2)$ proposals; every proposer gets its best and every
  receiver its worst partner over all stable matchings, and all
  stable matchings match the same agents.]

=== Connectivity and orientations
- Kahn's algorithm leaves out exactly the vertices reachable
  from a cycle.
- A digraph in which every edge lies on a cycle, for example a
  balanced one ($"in" = "out"$ everywhere), has strongly connected
  weak components.
- Period: in a strongly connected digraph take BFS depths $d$
  from any vertex and $g = gcd_(u -> v) (d(u) + 1 - d(v))$. Every
  closed walk has length divisible by $g$, $g$ is the gcd of the
  cycle lengths, and labels $d mod g$ make each edge go from $t$ to
  $t + 1$. A closed walk of length not divisible by $k$ exists iff
  $k$ does not divide $g$.
- Give each edge to one of its ends, each vertex getting at most
  one: possible iff every component has $m_c <= n_c$; a tree
  component allows $n_c$ ways, one with a cycle 2 ways (1 if the
  cycle is a self-loop). Hakimi: an orientation with every
  out-degree $<= k$ exists iff every vertex set $X$ spans at most
  $k |X|$ edges (flow: source $->$ edge, edge $->$ both ends,
  vertex $->$ sink with capacity $k$). Mixed graph: an Euler
  circuit exists iff it is connected, every degree is even and
  that flow on the undirected edges, with vertex $->$ sink
  capacity $deg(v)\/2 - "out"(v)$ (directed out-edges only; all
  $>= 0$), saturates them all; orient each one out of the vertex
  it flows to.

=== Counting
#fact[Matrix-tree][Spanning trees of an undirected multigraph:
  $det$ of the Laplacian $L = D - A$ without the row and column of
  any one vertex (parallel edges count, self-loops never do); with
  weights, the sum over trees of the weight products. Directed,
  arborescences out of root $r$: in-degrees on the diagonal (into
  $r$: out-degrees), $L_(u v) = -c(u -> v)$, row and column $r$
  deleted. $det$ mod $p$: `linalg/gauss-mod`.]
#fact[BEST][Eulerian digraph: the Euler circuits (up to rotation)
  number $t_w product_v ("out"(v) - 1)!$, with $t_w$ the spanning
  trees directed into any fixed $w$. Starting at $s$: times
  $"out"(s)$.#full[ So there are $(k!)^(k^(n - 1)) \/ k^n$ cyclic
  de Bruijn sequences of order $n$ over $k$ letters.]]

=== Degree sequences
#fact[Erdős-Gallai][$d_1 >= dots >= d_n$ is the degree sequence of
  a simple graph iff $sum d_i$ is even and for every $k$:
  $sum_(i <= k) d_i <= k(k - 1) + sum_(i > k) min(d_i, k)$.]
#fact[Havel-Hakimi][Build it: join a vertex of largest remaining
  degree $d$ to the $d$ next largest, repeat; it fails iff the
  sequence is not graphic.]
#fact[Gale-Ryser][Bipartite with degrees $a_1 >= dots >= a_p$ and
  $b_1, dots, b_q$ (a 0/1 matrix with these row and column sums)
  iff $sum a = sum b$ and $sum_(i <= k) a_i <= sum_j min(b_j, k)$
  for every $k$.]
#fact[Landau][$s_1 <= dots <= s_n$ are the out-degrees of a
  tournament iff $sum_(i <= k) s_i >= binom(k, 2)$ for every $k$,
  with equality at $k = n$.]

=== Colouring and extremal facts
- Vizing: a simple graph needs $Delta$ or
  $Delta + 1$ edge colours (a multigraph up to
  $floor(3 Delta \/ 2)$); bipartite, even with parallel edges,
  exactly $Delta$.
- Chromatic number, $n <= 20$: $chi <= k$ iff
  $sum_(S subset.eq V) (-1)^(n - |S|) i(S)^k > 0$, with $i(S)$ the
  number of independent subsets of $S$ (subset sums; count modulo
  a random large prime).
- Turán: without $K_(r + 1)$ at most $(1 - 1\/r) n^2\/2$ edges
  (balanced complete $r$-partite); triangle-free:
  $floor(n^2\/4)$.
- Ramsey: $R(3, 3) = 6$, $R(3, 4) = 9$, $R(3, 5) = 14$,
  $R(4, 4) = 18$, $R(4, 5) = 25$;
  $R(s, t) <= binom(s + t - 2, s - 1)$.
- Planar: $n - m + f = 1 + c$ ($c$ components, $f$ faces with the
  outer one), so $c = n - m +$ (bounded faces); in a grid drawing
  the bounded faces are the regions of unit cells cut off from
  the outside. Simple, $n >= 3$: $m <= 3n - 6$, and $m <= 2n - 4$
  without triangles (outerplanar: $m <= 2n - 3$); some vertex has
  degree $<= 5$.
- Grids: free cells under 4-adjacency are cut apart exactly by
  chains of blocked cells under 8-adjacency, the outside counting
  as blocked.#full[ Planar minimum cuts: _Flow and cut
  modelling_.]

=== Shortest paths
- Multi-source Dijkstra with $"own"(v)$, the source that reaches
  $v$ first (undirected): the nearest other source of $s$ is at
  distance
  $min (d(u) + w + d(v))$ over edges $(u, v, w)$ with
  $"own"(u) = s != "own"(v)$; over all edges, the closest pair of
  sources.
- Event times: $T_v$ is the least of a start value and $f(T_u)$
  over rules $u -> v$, each $f$ non-decreasing with $f(x) >= x$
  (adding a cost $>= 0$, waiting for a departure). Dijkstra
  computes them, always finalising the smallest tentative value.
- Steiner tree, terminals $t_0, dots, t_(k - 1)$: $"dp"[S][v]$ is
  the lightest tree joining $v$ and the terminals in $S$, with
  $"dp"[{i}][t_i] = 0$. For each $S$ by increasing mask, first
  $"dp"[S][v] = min_(T subset S) "dp"[T][v] + "dp"[S without T][v]$,
  then Dijkstra seeded with all of $"dp"[S][dot]$. The answer is
  $min_v "dp"[{0, dots, k - 1}][v]$, $O(3^k n + 2^k m log m)$.

=== Minimum spanning trees
- All MSTs share the multiset of weights; the MST path minimizes
  the largest edge between its ends. An edge of weight $w$ is in
  some MST iff its ends are apart in the DSU of the edges lighter
  than $w$, and in every MST iff it is also a bridge of the graph
  on those DSU components formed by the edges of weight $w$.
- Borůvka: every component takes its cheapest outgoing edge (ties
  by index), so $O(log n)$ rounds. Implicit dense graphs (XOR MST,
  Manhattan MST) need per round one query per vertex for its
  cheapest edge out of its component (Another colour in
  _Designing a segment tree node_).
- Connected, integer weights in $[0, W]$: the MST weighs
  $sum_(v = 0)^(W - 1) (c_v - 1)$, with $c_v$ the number of
  components of the edges of weight $<= v$.
- Shifted edge groups: edges $(u + k, v + k)$ for all $k < ell$ at
  once. Keep one DSU per $j$ where $a tilde_j b$ means
  $a + k tilde b + k$ for all $k < 2^j$. Cover $ell$ by two blocks
  of length $2^j$; a union at level $j$ that merges pushes itself
  down as two unions at level $j - 1$ (offsets $0$, $2^(j - 1)$),
  one that finds the pair joined stops. Each level merges at most
  $n - 1$ times: $O(n log n)$ unions. In Kruskal (groups by cost)
  a level-0 merge is a tree edge.
- Two groups maximizing the lightest edge inside a group: add
  edges by increasing weight to a parity DSU; the first one that
  closes an odd cycle is the answer, and the parity colouring at
  that moment attains it.

=== Cycles
#full[
- Shortest cycle, unweighted: BFS from every vertex; a non-tree
  edge $(u, v)$ closes a walk of length $d(u) + d(v) + 1$, $O(n m)$
  (the minimum over all roots is the girth). Weighted: for each
  edge $(u, v, w)$, $w$ + shortest $u$-$v$ path without it.
  Unweighted, a shortest cycle through a given edge or vertex has
  no chord.
]
- Minimum mean cycle (Karp): $d_k (v)$ = lightest walk of exactly
  $k$ edges ending at $v$, with $d_0 = 0$ everywhere. The answer
  is $min_v max_(k < n) (d_n (v) - d_k (v))\/(n - k)$ over $v$
  with finite $d_n (v)$, $O(n m)$.
- Triangles in $O(m sqrt(m))$: orient each edge towards the
  endpoint of larger (degree, id), then intersect out-lists.
- Cycle space: $m - n + c$ independent cycles. XOR (or sum mod 2)
  of any $u$-$v$ path = tree path XOR a combination of
  fundamental cycles: XOR basis of the cycle values.

=== Building graphs
- Edges between a vertex $x$ and a whole index range: give a
  segment tree extra vertices. Down tree (node $->$ children,
  leaf $->$ vertex): $x ->$ the $O(log n)$ nodes covering the
  range. Up tree (vertex $->$ leaf, child $->$ parent): those
  nodes $-> x$. 2-SAT "$x$ forbids every $y$ in a range": $x ->$
  down nodes whose leaves lead to $not y$, and $y ->$ up leaf, up
  nodes $-> not x$, so every implication keeps its contrapositive.

#full[
=== More facts
- Whitney: vertex connectivity $<=$ edge connectivity $<=$
  minimum degree. $k$-connected iff every two vertices are joined
  by $k$ internally disjoint paths. A 2-connected graph with
  $n >= 4$ has a cycle of length $>= 4$, so if every cycle is a
  triangle, every block is an edge or a triangle.
- st-orientation: for $s != t$ the edges can be oriented without
  cycles so that $s$ is the only source and $t$ the only sink iff
  adding the edge $s t$ makes the graph 2-connected (so any
  $s != t$ in a 2-connected graph). An st-numbering (Even-Tarjan,
  $O(n + m)$) gives it: orient from smaller to larger number.
- Sparse certificate: $F_1$ a spanning forest of $G$, $F_2$ one of
  $G - F_1$, and so on. $F_1 union dots union F_k$ has at most
  $k(n - 1)$ edges and the same $min(k, lambda(u, v))$ for all
  pairs ($lambda$ = edge connectivity); $k = 2$ keeps the bridges
  and the 2-edge-connected components. Edges $1..r$, forests built
  newest edge first (_Matroids_, latest basis): their edges with
  index $>= l$ certify the edges in $[l, r]$.
- Series-parallel: a graph has no $K_4$ minor (treewidth $<= 2$)
  iff it vanishes under: delete a vertex of degree $<= 1$, delete a
  self-loop, merge parallel edges, replace a degree-2 vertex and
  its two edges by one edge. Heaviest simple cycle: series adds
  the weights; merging parallel $a$, $b$ records the 2-cycle
  $a + b$ and keeps $max(a, b)$; a self-loop is a cycle too.
- Vertex elimination (delete $v$, join its remaining neighbours
  pairwise): after eliminating a set $D$, in any order, two
  remaining vertices are adjacent iff some path joins them with
  all inner vertices in $D$ (fill path lemma). In a tree: iff
  every inner vertex of their path is in $D$.
- Inversion graph of a permutation ($i < j$ adjacent iff
  $p_i > p_j$): cliques are decreasing and independent sets
  increasing subsequences. It is bipartite iff no decreasing
  subsequence has length 3, and then the LIS is $n - nu$ (Kőnig);
  if it is a tree, it is a caterpillar (a path plus leaves).
- Grid cut vertices (4-adjacency): walk around the 8 neighbours of
  a free cell $v$; cells outside the grid are blocked, and so is a
  free diagonal cell between two blocked ones. $v$ is a cut vertex
  iff two different blocked stretches belong to the same
  8-connected blocked component, where all blocked cells on the
  border join the outside.
- Hamiltonian cycles ($n >= 3$): Dirac ($deg(v) >= n\/2$ for all
  $v$) and Ore ($deg(u) + deg(v) >= n$ for non-adjacent $u, v$)
  suffice; in general $O(2^n n^2)$ bitmask DP. Every tournament
  has a Hamiltonian path (insertion sort with the edges as
  comparator).
- Edge colouring a bipartite graph with $Delta$ colours:
  recolour along alternating paths, or split a $Delta$-regular
  supergraph into perfect matchings.
- Pairs (length, edges) compared lexicographically give the
  fewest edges among shortest paths.
- MST, more: the second-best MST swaps one edge: add a non-tree
  edge, drop the largest edge on its tree path (binary lifting).
  An edge is in every MST iff it is strictly lighter than every
  other edge crossing some cut, iff no cycle has it as a weakly
  heaviest edge. Kruskal reconstruction tree: _Tree techniques_.
- Directed shortest cycle through $s$: shortest $s -> u$ plus an
  edge $u -> s$; over all $s$, $O(n m)$ BFS or $n$ Dijkstras.
- Sources added one at a time (unweighted): BFS from each new
  source, pushing only strict improvements; a vertex is pushed at
  most $d_1 (v) + 1$ times in total, $d_1$ = distance to the first.
- Sequential overwrite on a functional graph $v -> A_v$, labels
  $l(v) in {1, dots, n}$: for $p = 1, dots, n$ in turn, all nodes
  labelled $p$ at that moment take the label of $A_v$ at once, read
  before any write of that turn. After this round $l'(v) = l(w)$,
  with $w$ the first node after $v$ on its forward path entered by
  a non-descent, $l(w) >= l("previous node")$. Every later round is
  the plain shift $l(v) <- l(A_v)$, so after round $t >= 1$ node $v$
  holds $l'(A^(t - 1)(v))$. Label $p$ survives round 1 iff some
  node labelled $p$ has a predecessor with label $<= p$. A cycle
  loses labels only then (same test along the cycle) and rotates
  after, so it can keep many forever: 1, 1, 2, 3 along a 4-cycle
  keeps all three. Turn games on pointers: one round, then jumps.

=== Cues
- Implicit dense graphs (XOR MST, Manhattan MST, "cheapest
  partner" queries): Borůvka, one data structure query per vertex
  and round that skips its own component (Another colour in
  _Designing a segment tree node_).
- A walk that must pass every edge, or fix degrees or visit
  counts mod 2: Euler tours plus parity fixing on a spanning tree.
- "Can return", closed walks of length mod $k$: strongly
  connected components and the period.
- Edges valid in an index range $[l, r]$, many queries: newest
  edge first spanning forests (latest basis).
]
