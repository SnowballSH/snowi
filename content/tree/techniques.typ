#import "/typst/notes.typ": *

Root the tree anywhere and work in BFS or DFS order, never by
recursion: the depth reaches $n$. Then ask what one vertex, edge or
path contributes.

=== Orders and paths
#fact[Euler tour][In DFS preorder `tin`, the subtree of $v$ is
  `[tin[v], tin[v] + size[v])`: subtree updates and sums become
  range operations, and $u$ is an ancestor of $v$ iff
  `tin[u] <= tin[v] < tin[u] + size[u]`. Among the vertices of one
  depth sorted by `tin`, the descendants of $v$ form a range.]
#fact[Root prefixes][With $S(v)$ the sum on the root path of $v$,
  the path $u - v$ sums to $S(u) + S(v) - 2 S(w) + a_w$ for
  $w = "lca"(u, v)$ (`tree/lca`). Changing $a_v$ changes $S$ on
  the subtree of $v$ (range add, point query); adding $x$ on the
  root path of $v$ is a point add at `tin[v]`, read by subtree
  sums. Updates and queries both on paths: `tree/hld`.]
#fact[Light edges][An edge to a child other than the largest one
  at least halves the subtree size, so a root path has at most
  $log_2 n$ light edges: `tree/hld` cuts any path into $O(log n)$
  ranges, and small-to-large moves each vertex $O(log n)$ times.]
#fact[Binary lifting][The $2^k$-th ancestors give the $k$-th
  ancestor, the LCA, and the highest ancestor where a property
  monotone along the root path still holds (jump while it holds,
  largest $k$ first). Jump pointers over levels: store the
  ancestor $B$ levels up (downwards: the `tin` range $B$ levels
  below) and split a shift $t = t_1 B + t_0$ into $t_1$ big jumps
  and $t_0$ small ones.]

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
#fact[Rerooting][The answer for every root is a DP on directed
  edges, each side of each edge once; prefix and suffix merges
  avoid inverses (`tree/rerooting`).]
#fact[Tree knapsack][Merging child tables of sizes $s_1$ and $s_2$
  in $s_1 s_2$ steps costs $O(n^2)$ in total: each pair of
  vertices meets once, at its LCA. Tables capped at $K$:
  $O(n K)$.]

=== Shape
#fact[Diameter][The farthest vertex $a$ from any vertex ends a
  diameter, and the farthest $b$ from $a$ ends it (BFS twice; DFS
  with weights $>= 0$). Eccentricity $= max(d(a, v), d(b, v))$.
  The center (middle one or two vertices of $a - b$) lies on every
  diameter; the radius is $ceil(D \/ 2)$. Joining two trees by an
  edge between their centers gives the smallest diameter,
  $max(D_1, D_2, ceil(D_1 \/ 2) + ceil(D_2 \/ 2) + 1)$.]
#fact[Centroid][Removing it leaves parts of at most $n \/ 2$
  vertices; there are one or two (adjacent). Walk from any vertex
  into a part larger than $n \/ 2$. It minimises the sum of
  distances. Removing centroids recursively takes at most
  $log_2 n + 1$ levels, and each path meets the first centroid
  removed on it (`tree/centroid`).]
#fact[Isomorphism][Rooted: ids of sorted child-id lists (AHU).
  Unrooted: root at the one or two centers (`tree/tree-hash`).]
#fact[Counting][Cayley: $n^(n - 2)$ labelled trees, $(n - 2)! \/
  product_i (d_i - 1)!$ with degrees $d_i$. Joining components of
  sizes $s_1, dots, s_k$ into a tree with $k - 1$ new edges:
  $n^(k - 2) product_i s_i$ ways. Prüfer codes: _Counting formulas_.]

=== Many queries
#fact[Virtual tree][Queries on key sets with $sum k$ bounded: the
  keys and the LCAs of keys adjacent in `tin` order ($<= 2k - 1$
  vertices) keep every branching; run the DP there with
  compressed edge lengths (`tree/virtual-tree`). Keys sorted by
  `tin` with $v_(k + 1) = v_1$: the smallest subtree containing
  them has $1/2 sum_i d(v_i, v_(i + 1))$ edges.]
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

=== Changing trees
#fact[Link-cut][Edges inserted and deleted online, with path
  folds, connectivity or a moving root: `tree/link-cut`. All
  changes known in advance and only connectivity or
  per-component questions: `offline/dynamic-connectivity` (a
  rollback DSU knows no paths).]
#fact[Edge, not vertices][When a subtree is cut and relinked
  elsewhere and all its vertices change alike (depth shift,
  subtree add), store the change once on the new edge or the
  subtree's root, and let root-path folds pick it up.]

#full[
=== Cues
- "For every subtree, a statistic of its colours": small to
  large. "For every vertex as the root": rerooting. "Count paths
  with a property", "nearest marked vertex": centroids.
  "$q$ queries of $k_i$ vertices, $sum k_i <= 2 dot 10^5$":
  virtual tree.
- A resource that only grows (rating, water level, time) and
  edges usable above a threshold: the reachable set is always a
  Kruskal reconstruction tree node; the process walks up it.
- "Each vertex becomes the XOR of itself and its children, $t$
  times": $(I + C)^t$ over $"GF"(2)$; the descendant $d$ levels
  below contributes iff $binom(t, d)$ is odd, iff $d$ is a submask
  of $t$ (Lucas). Split $t$ into high and low bits and lift over
  levels.
- $x mod a_i$ along a path changes $x$ only $O(log x)$ times (each
  change at least halves it): jump to the next $a_i <= x$ with a
  path-minimum descent (binary lifting, or a link-cut tree when
  edges change).
- $n$ up to $10^6$: flat arrays, BFS orders, fast input, an
  $O(1)$ LCA (`tree/lca-rmq`).

=== Pitfalls
- Sums of distances reach $n^2 \/ 2$: use 64-bit.
- Small to large: empty the whole subtree of a light child, not
  just the vertex; reset global maxima when the counts empty.
- Centroids: recompute subtree sizes inside the current piece;
  when counting pairs through a centroid, subtract the pairs
  inside one child part.
- Virtual tree: add the root as a key if the DP needs it; clear
  per-query arrays only on the virtual vertices.
- Kruskal reconstruction tree: same-component edges and
  self-loops create no node; a disconnected graph gives a forest.

=== Prüfer codes
Encode: repeatedly remove the smallest leaf and write down its
neighbour, $n - 2$ times. Decode: count $deg(v) = 1 +$ occurrences
of $v$; for each code entry $x$, join $x$ to the smallest $v$ with
$deg(v) = 1$, then decrease both degrees; the last two vertices
with degree 1 form the final edge. A min-heap gives
$O(n log n)$; walking a pointer over the smallest leaf gives
$O(n)$.
]
