#import "/typst/notes.typ": *

=== Orders and paths
#fact[Euler tour][Among the vertices of one depth sorted by `tin`,
  the descendants of $v$ form a range.]
#fact[Jump pointers][Over levels: store the ancestor $B$ levels
  up (downwards: the `tin` range $B$ levels below) and split a
  shift $t = t_1 B + t_0$ into $t_1$ big jumps and $t_0$ small
  ones.]
#fact[Upward sequences][To compare the values read upwards from
  $u$ and from $v$ (root paths, `nxt` chains), store $h_k (v)$, the
  polynomial hash of the $2^k$ values from $v$ up (while the chain
  has that many): $h_k (v) = h_(k - 1)(v) dot B^(2^(k - 1)) +
  h_(k - 1)("up"_(k - 1)(v))$. Jumping both while the hashes agree,
  largest $k$ first, gives the common prefix in $O(log n)$; the
  next values decide.]
#fact[Median][Exactly one vertex lies on all three paths between
  $u$, $v$ and $r$: the LCA of $u$ and $v$ with the tree rooted at
  $r$. With a fixed root it is the deepest of $"lca"(u, v)$,
  $"lca"(u, r)$, $"lca"(v, r)$ (two of them coincide), so queries
  under any root need no rerooting. For distinct $u, v, r$ it is the
  one vertex $w$ that puts them in three different parts: $w$
  itself and the components of $T - w$. With $A(x, y)$ the XOR of
  the labels on the path, both ends included,
  $A(u, v) plus.o A(v, r) plus.o A(r, u)$ is the median's label:
  every other vertex of the three paths lies on exactly two.]
#fact[Path intersection][Paths $a - b$ and $c - d$: let $p, q$ be
  the two deepest of $"lca"(a, c)$, $"lca"(a, d)$, $"lca"(b, c)$,
  $"lca"(b, d)$. The paths meet iff $"dep"(p) >=
  max("dep"("lca"(a, b)), "dep"("lca"(c, d)))$, and then they share
  exactly the path $p - q$.]

=== Subtrees and roots
#fact[Small to large][Keep the counts of the heavy child's subtree
  and add the light subtrees vertex by vertex: $O(n log n)$ calls
  of `add` and `remove` give any statistic of the multiset of each
  subtree (distinct colours, most frequent colour, depth counts).
  Iterative, with `tree/hld` (heavy child first in `pos`):]
#lang-code(
  cpp: "// HLD h(adj); at[h.pos[v]] = v for every vertex v.
for (int i = n; i--;) { // children before parents
  int v = at[i], hv = h.heavy[v]; // counts: subtree of hv
  rep(j, i + 1 + (hv < 0 ? 0 : h.size[hv]), i + h.size[v])
    add(at[j]); // the light subtrees
  add(v), ans[v] = get();
  if (h.head[v] == v) // not a heavy child: empty the counts
    rep(j, i, i + h.size[v]) remove(at[j]);
}",
  py: "at = [0] * n  # h = HLD(adj): heavy child first in pos
for v in range(n):
  at[h.pos[v]] = v
for i in range(n - 1, -1, -1):  # children before parents
  v = at[i]
  hv = h.heavy[v]  # the counts hold the subtree of hv
  lo = i + 1 + (h.size[hv] if hv >= 0 else 0)
  for j in range(lo, i + h.size[v]):  # light subtrees
    add(at[j])
  add(v)
  ans[v] = get()
  if h.head[v] == v:  # not a heavy child: empty the counts
    for j in range(i, i + h.size[v]):
      remove(at[j])",
)
#fact[Merging segment trees][Dynamic segment trees over one value
  range (children made on demand, 0 = empty): `merge(a, b)` returns
  the other if one is 0, combines two leaves, else merges the
  children pairwise and pulls. Each step destroys a node, so all
  merges of $n$ one-element trees cost $O(n log V)$ in total. For
  subtrees, merge the children into the parent; answer a vertex's
  queries before its tree is merged upwards (nodes are reused).]
#fact[Long paths][For DP arrays indexed by the depth below $v$
  (length $h_v + 1$, $h_v$ the height): $v$ inherits the array of
  its tallest child by an index shift (one buffer per long chain,
  $v$ one slot before that child; lazy tags for global adds), then
  adds every other child $c$ entry by entry in $O(h_c + 1)$. Each
  chain is added once, so these costs sum to at most $n$: $O(n)$ in
  total.]
#fact[Tree knapsack][Merging the tables of two parts with $s_1$
  and $s_2$ vertices in $s_1 s_2$ steps costs $O(n^2)$ in total:
  each pair of vertices meets once, at its LCA. Tables capped at
  $K$: $O(n K)$. Size each merged table $min(s_1 + s_2, K) + 1$
  and loop only over existing entries, never up to $n$ or $K$.]

=== Shape and covers
#fact[Diameter][#full[For a diameter $a - b$ (BFS twice; weights
  $>= 0$), the eccentricity of $v$ is $max(d(a, v), d(b, v))$, and
  every diameter holds a center. ]Joining two trees by an edge of
  length $w$ between their centers (vertices of least eccentricity
  $r$) gives the smallest diameter, $max(D_1, D_2, r_1 + w + r_2)$;
  unweighted, $r = ceil(D \/ 2)$.#full[ The centers are the middle
  one or two vertices of every diameter.]]
#fact[Centroid][With weights $w >= 0$ on the vertices, a vertex
  whose removal leaves parts of weight at most $W \/ 2$ minimises
  $sum_x w_x d(x, v)$ (unit weights: the sum of distances).]
#full[
#fact[Weighted medians][The medians, minimisers of
  $F(s) = sum_x w_x d(x, s)$ ($w >= 0$, lengths $> 0$), are exactly
  the vertices of Centroid. Crossing an edge from side $A$ to side
  $B$ changes $F$ by its length times $W(A) - W(B)$, which never
  falls along a path: $F$ is convex on paths and never rises towards
  a median. If $W > 0$ the medians form a path. Moving each vertex
  to its nearest one in a subtree never increases a distance, nor
  $F$ if the subtree holds a median. So
  $F_1(a) + F_2(b) + c dot d(a, b)$, $c >= 0$, has its closest
  optimal pair on the path between any medians $m_1$, $m_2$. Walk
  $a$ from $m_1$ and $b$ from $m_2$ along it, each while its next
  step costs it at most $c$ per unit of length. If they stop apart
  that is the pair; if they meet or cross, $a = b$ at the minimum of
  $F_1 + F_2$ on the stretch both walked. This works for any $G$,
  $H$ convex on a path, from minimisers $a_0 <= b_0$.]
]
#fact[Binary-tree split][At most two children per vertex and
  $W >= 2$ marked vertices, none an ancestor of another (say leaves,
  or one depth): some subtree holds between $W \/ 3$ and $2 W \/ 3$.
  Walk into the heavier child down to the first $v$ holding at most
  $W \/ 2$; take $v$ if it holds at least $W \/ 3$, else its parent.
  Questions "is it in the subtree of $x$?" then find a marked
  vertex in $ceil(log_(3\/2) W)$ steps.]
#fact[Hitting subtrees][Fewest vertices meeting every set of a
  family of connected vertex sets: go through the sets by
  decreasing depth of their top (shallowest) vertex and pick the
  top of each set not hit yet. The count equals the largest number
  of pairwise disjoint sets in the family.]
#fact[Distance-$k$ cover][Fewest centres with every vertex within
  distance $k$ of one: take the deepest uncovered vertex, put a
  centre at its $k$-th ancestor (the root if shallower), cover all
  within $k$ of it, repeat. In $O(n)$ bottom-up: keep the largest
  distance to an uncovered vertex and the smallest to a centre in
  the subtree; place a centre when the first reaches $k$ and the
  second cannot cover it, and one at the root if anything is left.]
#fact[Cube][List the vertices in DFS order, each on entry if its
  depth is even and on exit if it is odd: cyclically consecutive
  ones are at distance at most 3. So the cube of a connected graph
  has a Hamiltonian cycle (Sekanina; use a spanning tree).]
#fact[Joining components][Adding $k - 1$ edges to join components
  of sizes $s_1, dots, s_k$ ($n$ vertices) into one tree:
  $n^(k - 2) product_i s_i$ ways.#full[ Cayley, degrees, Prüfer
  codes: _Counting formulas_.]]

=== Many queries
#fact[Key sets][#icpc[For keys $v_1, dots, v_k$ sorted by `tin`,
  the smallest subtree containing them has
  $1/2 sum_i d(v_i, v_(i + 1))$ edges ($v_(k + 1) = v_1$), and the
  union of their root paths has
  $sum_i "dep"(v_i) - sum_(i < k) "dep"("lca"(v_i, v_(i + 1)))$
  edges (or weight, with weighted depths).]#full[The smallest
  subtree containing keys and the union of their root paths are
  sums over keys adjacent in `tin` order (`tree/virtual-tree`),
  also with weighted depths.]
  An ordered set by `tin` keeps either under insertions and
  deletions in $O(log n)$; with small to large over another tree,
  for each of its subtrees.]
#fact[Centroid queries][$min_x ("val"_x + d(v, x))$ over a changing
  set: each centroid $c$ keeps a multiset $M_c$ of
  $"val"_x + d(x, c)$ over its piece, and the answer is the minimum
  of $d(v, c) + min M_c$ over the centroid ancestors $c$ of
  $v$.#full[ Exact: every term is at least the truth, and the first
  centroid removed on the path to the best $x$ attains it.] Only
  $x$ with `tin` in a range (one branch of $v$, or the rest): a
  segment tree per centroid over its piece sorted by `tin`,
  $O(log^2 n)$ per query or update.]
#fact[Kruskal reconstruction tree][For "only edges of weight
  $<= x$": add edges by increasing weight; the union of roots $a$
  and $b$ creates a node $z$, parent of both, with $"val"(z) = w$.
  From $v$ with threshold $x$ one reaches exactly the leaves under
  its highest ancestor with $"val" <= x$ (binary lifting), and
  the bottleneck of $u$ and $v$ (the smallest possible largest
  edge on a path) is $"val"("lca"(u, v))$. Ids grow towards the
  root: decreasing ids are a top-down order.]
#lang-code(
  cpp: "// Leaves 0..n-1; es = {w, u, v} sorted by w. par = -1 at
// roots of the forest.
vi dsu(2 * n), par(2 * n, -1);
vector<ll> val(2 * n);
iota(all(dsu), 0);
auto find = [&](int x) {
  while (dsu[x] != x) x = dsu[x] = dsu[dsu[x]];
  return x;
};
int z = n;
for (auto [w, u, v] : es)
  if (int a = find(u), b = find(v); a != b)
    dsu[a] = dsu[b] = par[a] = par[b] = z, val[z++] = w;",
  py: "# Leaves 0..n-1; es = (w, u, v) triples. par = -1 at roots.
dsu, par = list(range(2 * n)), [-1] * (2 * n)
val = [0] * (2 * n)
def find(x):
  while dsu[x] != x:
    dsu[x] = dsu[dsu[x]]
    x = dsu[x]
  return x
z = n
for w, u, v in sorted(es):
  a, b = find(u), find(v)
  if a != b:
    dsu[a] = dsu[b] = par[a] = par[b] = z
    val[z] = w
    z += 1",
)
#fact[Records][Walking up and paying only at each new maximum:
  let `nxt[v]` be the first proper ancestor with a larger value.
  Values increase along `nxt` chains, so costs are prefix sums on
  the chain and the first record above a level is found by
  lifting on `nxt`. Build `nxt` top-down by lifting from the
  parent; walking the chain step by step is quadratic on trees.]

#full[
=== Changing trees
#fact[Edge, not vertices][When a subtree is cut and relinked
  elsewhere and all its vertices change alike (depth shift,
  subtree add), store the change once on the new edge or the
  subtree's root, and let root-path folds pick it up.]
#fact[Growing trees][Leaves attached over time to parents
  named in the input: build the final tree once (Euler tour, LCA,
  HLD) and give vertices not created yet neutral values.]
]

#full[
=== More facts
#fact[Leaf pairings][Pair up an even number of leaves and route
  each pair along its path. The possible edge loads $b_e$ are
  exactly those with $b_e = 1$ on leaf edges, $b_e$ of the parity
  of the number of leaves on one side of $e$, and at every vertex
  the largest incident load at most the sum of the others.]
#fact[Disjoint pieces][Two vertex-disjoint connected pieces are
  joined by a path with $j$ inner vertices and $j + 1$ edges. So
  the pairs that deleting $e$ separates, summed over the edges $e$,
  minus the pairs avoiding $v$ that deleting $v$ separates, summed
  over the vertices $v$, count every disjoint pair once (pairs of
  vertex sets whose spanning subtrees are disjoint).]

=== Cues
- A resource that only grows (rating, water level, time) and
  edges usable above a threshold: the reachable set is always a
  Kruskal reconstruction tree node; the process walks up it.
- "LCA when rooted at $r$", "the vertex where three paths meet":
  the median. "Kill every path or connected set with the fewest
  vertices": hitting subtrees. "Every vertex within distance $k$
  of a chosen one": distance-$k$ cover.
- Depth-indexed DP ("pairs at distance $k$", "vertices $d$ levels
  below"): long paths, not small to large.
- "Each vertex becomes the XOR of itself and its children, $t$
  times": $(I + C)^t$ over $"GF"(2)$; the descendant $d$ levels
  below contributes iff $binom(t, d)$ is odd, iff $d$ is a submask
  of $t$ (Lucas). Split $t$ into high and low bits and lift over
  levels.
- A walk goes to its largest unvisited neighbour, ties your way, and
  stops at the root $r$. It reaches $r$ from every start iff
  $w_x <= w_g$ for every $x$ of depth $>= 2$ with grandparent $g$
  ($<$ if ties go against you): the start at $y$ needs the parent of
  $y$ at least every child of $y$. Lowering values to fit, the best
  $w'_x$ is the minimum of $w$ on every other vertex of the path
  from $x$ to $r$, starting at $x$ (all roots: `tree/centroid`).
- Per centroid $c$: to count paths of length at most $K$ no second
  sort is needed, since one stable pass over the piece's sorted
  list, appending each vertex to its part's list, leaves every part
  sorted. With integer lengths the least weight on the path from
  $r$ to $c$ at even distance from $x$ is one of two values kept
  per $r$, the minima by the parity of $d(y, c)$: take the one
  matching $d(c, x)$.
- $x mod a_i$ along a path changes $x$ only $O(log x)$ times (each
  change at least halves it): jump to the next $a_i <= x$ with a
  path-minimum descent (binary lifting, or a link-cut tree when
  edges change).
- $n$ up to $10^6$: flat arrays, BFS orders, fast input, an
  $O(1)$ LCA (`tree/lca-rmq`).
- Shuffled labels make long paths jump around memory.#cpp[
  Relabelling in BFS or DFS order and rebuilding `adj` vertex by
  vertex in that order made `tree/centroid` 2.5 to 5 times faster
  on paths and 1.4 to 2 times on random trees, $n$ from
  $2 dot 10^5$ to $10^6$, with no gain on a star. Rebuilt from the
  edge list, paths gain about 2 times, random trees at most 1.3.]#py[
  `tree/centroid` relabels in BFS order inside, which gains 1.5 to
  2 times on paths and little on random trees at $n = 2 dot 10^5$.]

=== Pitfalls
- The depth reaches $n$: work in BFS or DFS order, never by
  recursion.
- Sums of distances reach $n^2 \/ 2$: use 64-bit.
- All children but one: prefix and suffix folds, never division
  (a factor can be $0$ mod $p$).
- Small to large: empty the whole subtree of a light child, not
  just the vertex; reset global maxima when the counts empty.
- Centroids: recompute subtree sizes inside the current piece;
  when counting pairs through a centroid, subtract the pairs
  inside one child part.
- Virtual tree: add the root as a key if the DP needs it; clear
  per-query arrays only on the virtual vertices.
- Kruskal reconstruction tree: same-component edges and
  self-loops create no node; a disconnected graph gives a forest.
]
