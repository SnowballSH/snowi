#import "/typst/notes.typ": *

$n$ vertices, $m$ edges, $Delta$ the maximum degree. $nu$, $alpha$,
$tau$, $rho$: maximum matching, maximum independent set, minimum
vertex cover, minimum edge cover.

=== Degrees and Euler tours
#fact[Handshake][$sum_v deg(v) = 2m$, so the number of odd-degree
  vertices is even. Directed: $sum "in" = sum "out" = m$.]
#fact[Euler][Ignore isolated vertices; the rest must be connected.
  Undirected: a circuit iff every degree is even, a path iff exactly
  0 or 2 degrees are odd (it joins the odd ones). Directed: a
  circuit iff $"in"(v) = "out"(v)$ everywhere, a path iff besides
  that one vertex has $"out" = "in" + 1$ (the start) and one
  $"in" = "out" + 1$ (the end). Connected, $k > 0$ odd vertices:
  the edges split into $k\/2$ trails and no fewer (`graph/euler-path`
  after pairing the odd vertices with dummy edges).]

=== Matchings, covers and orders
#fact[Gallai][$alpha + tau = n$ in every graph (the complement of
  a vertex cover is independent). $nu + rho = n$ if there is no
  isolated vertex (extend a maximum matching greedily).]
#fact[Kőnig][Bipartite: $nu = tau$, so $alpha = rho = n - nu$
  (no isolated vertices for $rho$).]
#fact[Hall][Bipartite with sides $L, R$: $L$ can be matched
  completely iff $|N(X)| >= |X|$ for all $X subset.eq L$;
  $nu = |L| - max_X (|X| - |N(X)|)$.]
#fact[Tutte][A perfect matching exists iff removing any set $U$
  leaves at most $|U|$ odd components. Tutte-Berge:
  $nu = min_U (n + |U| - "odd"(G - U))\/2$.]
#fact[Dilworth, Mirsky][In a partial order, the largest antichain
  equals the fewest chains covering everything ($n - nu$ on the
  strict relation), and the longest chain equals the fewest
  antichains covering everything (layer by height). On sequences:
  the fewest non-increasing subsequences covering a sequence is
  the length of its longest strictly increasing subsequence.]
#fact[Menger, max-flow min-cut][The most edge-disjoint $s$-$t$
  paths = the fewest edges separating $s$ from $t$; with
  vertex-disjoint paths and vertex cuts likewise ($s, t$ not
  adjacent). In general max flow = min cut (`flow/modelling`).]

=== Counting
#fact[Matrix-tree][Spanning trees of a multigraph: $det$ of the
  Laplacian $D - A$ without one row and its column. Directed, out
  of root $r$: diagonal = in-degrees, $L_(u v) = -c(u, v)$, delete
  row and column $r$; into $r$: out-degrees instead.]
#fact[BEST][Eulerian digraph: the Euler circuits (up to rotation)
  number $t_w product_v ("out"(v) - 1)!$, with $t_w$ the spanning
  trees directed into any fixed $w$. Starting at $s$: times
  $"out"(s)$.]
#full[
#fact[Cayley][$n^(n - 2)$ labelled trees; forests of $k$ trees
  rooted at $k$ given vertices: $k n^(n - k - 1)$; degrees $d_i$:
  $(n - 2)! \/ product (d_i - 1)!$.]
]

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
- Bipartite iff no odd cycle: 2-colour by BFS; an edge inside one
  colour plus the two tree paths is an odd cycle.
- Brooks: $chi <= Delta$ for a connected graph that is neither
  complete nor an odd cycle; greedy colouring never needs more
  than $Delta + 1$. Vizing: a simple graph needs $Delta$ or
  $Delta + 1$ edge colours (a multigraph up to
  $floor(3 Delta \/ 2)$); bipartite, even with parallel edges,
  exactly $Delta$.
- Turán: without $K_(r + 1)$ at most $(1 - 1\/r) n^2\/2$ edges
  (balanced complete $r$-partite); triangle-free:
  $floor(n^2\/4)$.
- Ramsey: $R(3, 3) = 6$, $R(3, 4) = 9$, $R(3, 5) = 14$,
  $R(4, 4) = 18$, $R(4, 5) = 25$;
  $R(s, t) <= binom(s + t - 2, s - 1)$. So any 6 vertices hold a
  triangle or an independent triple.
- Planar: $n - m + f = 1 + c$ ($c$ components, $f$ faces with the
  outer one). Simple, $n >= 3$: $m <= 3n - 6$, and
  $m <= 2n - 4$ without triangles; some vertex has degree
  $<= 5$; 4 colours suffice. $K_5$ and $K_(3,3)$ are not planar.

=== Trees and DAGs
- Any two of connected, acyclic, $m = n - 1$ imply the third; a
  tree has exactly one path between two vertices, and
  $sum deg = 2n - 2$.
- Diameter: the farthest vertex from any vertex is an end of a
  diameter; BFS (or DFS with weights $>= 0$) twice. Centre: the
  middle one or two vertices of a diameter.
- A topological order exists iff there is no directed cycle; it
  is unique iff consecutive vertices are joined by edges. Longest
  path, path counts: DP in topological order.

=== Shortest paths and MST
- Dijkstra needs weights $>= 0$. Johnson: potentials
  $h = $ Bellman-Ford distances from a new vertex joined by
  0-edges to all; $w(u, v) + h(u) - h(v) >= 0$ keeps shortest
  paths and $d(u, v) = d'(u, v) - h(u) + h(v)$.
- Weights 0/1: 0-1 BFS. Bellman-Ford: a shortest path has at most
  $n - 1$ edges; an improvement in round $n$ means a negative
  cycle. Floyd: negative cycle iff some $d(i, i) < 0$.
- Difference constraints $x_v - x_u <= w$: edge $u -> v$ of weight
  $w$; feasible iff no negative cycle, $x$ = distances from a
  new source joined by 0-edges.
- Cut property: a lightest edge across any cut is in some MST.
  Cycle property: a strictly heaviest edge of a cycle is in no
  MST. All MSTs share the multiset of weights; distinct weights
  give one MST. The MST path minimizes the largest edge between
  its ends (bottleneck).
- Borůvka: every component takes its cheapest outgoing edge;
  $O(log n)$ rounds. Use it when edges are implicit (XOR MST,
  Manhattan MST, "cheapest partner" queries per round).

=== Cycles
- Shortest cycle, unweighted: BFS from every vertex; a non-tree
  edge $(u, v)$ closes a walk of length $d(u) + d(v) + 1$, $O(n m)$
  (the minimum over all roots is the girth). Weighted: for each
  edge $(u, v, w)$, $w$ + shortest $u$-$v$ path without it.
- Triangles in $O(m sqrt(m))$: orient each edge towards the
  endpoint of larger (degree, id), then intersect out-lists.
- Cycle space: $m - n + c$ independent cycles. XOR (or sum mod 2)
  of any $u$-$v$ path = tree path XOR a combination of
  fundamental cycles: XOR basis of the cycle values.

#full[
=== More facts
- Whitney: vertex connectivity $<=$ edge connectivity $<=$
  minimum degree. $k$-connected iff every two vertices are joined
  by $k$ internally disjoint paths.
- Erdős-Szekeres: more than $r s$ distinct numbers contain an
  increasing subsequence of length $r + 1$ or a decreasing one of
  length $s + 1$.
- Minimum path cover of a DAG by vertex-disjoint paths:
  $n - nu$ in the split graph (out-copy to in-copy); paths may
  share vertices: the same on the transitive closure.
- Euler, more: the fewest trails covering a connected digraph's
  edges: $max(1, sum_v max(0, "out"(v) - "in"(v)))$. Chinese
  postman (closed walk using every edge): add a minimum-weight
  perfect matching of the odd vertices under shortest-path
  distances (bitmask DP for up to about 20 of them).
- De Bruijn: an Euler circuit of the graph on the $k^(n - 1)$
  words of length $n - 1$ gives a cyclic word of length $k^n$
  containing every word of length $n$ once; there are
  $(k!)^(k^(n - 1)) \/ k^n$ of them.
- Hamiltonian cycles ($n >= 3$): Dirac ($deg(v) >= n\/2$ for all
  $v$) and Ore ($deg(u) + deg(v) >= n$ for non-adjacent $u, v$)
  suffice; in
  general $O(2^n n^2)$ bitmask DP. Every tournament has a
  Hamiltonian path (insertion sort with the edges as comparator).
- Chromatic number for $n <= 20$: $chi <= k$ iff
  $sum_(S subset.eq V) (-1)^(n - |S|) i(S)^k > 0$, with $i(S)$ the
  number of independent subsets of $S$ (subset sums, modulo a
  prime). Edge
  colouring a bipartite graph with $Delta$ colours: recolour along
  alternating paths, or split a $Delta$-regular supergraph into
  perfect matchings.
- Planar duals: with $s$ and $t$ on the outer face, a minimum
  $s$-$t$ cut is a shortest dual path between the two parts of
  the outer face that $s$ and $t$ separate (grids: a blocking
  wall is a path of cells).
  Outerplanar: $m <= 2n - 3$.
- Trees with $n >= 2$ have $2 + sum_(deg(v) >= 3) (deg(v) - 2)$
  leaves. A centroid (every component after removing it has
  $<= n\/2$ vertices) exists; there are one or two, and two are
  adjacent. Pairwise distance sum $= sum_e s_e (n - s_e)$, $s_e$
  the size of one side. $K_(p, q)$ has $p^(q - 1) q^(p - 1)$
  spanning trees.
- MST, more: the second-best MST swaps one edge: add a non-tree
  edge, drop the largest edge on its tree path (binary lifting).
  An edge is in every MST iff it is strictly lighter than every
  other edge crossing some cut, iff no cycle has it as a weakly
  heaviest edge. The Kruskal reconstruction tree answers "lowest
  possible maximum edge" queries as an LCA.
- Minimum mean cycle (Karp): $d_k (v)$ = lightest walk of exactly
  $k$ edges ending at $v$, with $d_0 = 0$ everywhere. The answer
  is $min_v max_(k < n) (d_n (v) - d_k (v))\/(n - k)$ over $v$
  with finite $d_n (v)$, $O(n m)$.
- Directed shortest cycle through $s$: shortest $s -> u$ plus an
  edge $u -> s$; over all $s$, $O(n m)$ BFS or $n$ Dijkstras.
]
