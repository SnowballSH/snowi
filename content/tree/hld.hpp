/**
 * Title: Heavy-light decomposition
 * Description: Lays the tree out in an array so that every path is
 *   $O(log n)$ ranges of consecutive positions. `pos[v]` is a DFS
 *   order that enters the heavy child (largest subtree) first, so
 *   each chain from its top `head[v]` down is contiguous and the
 *   subtree of `v` is `[pos[v], pos[v] + size[v])`. `path(u, v)`
 *   returns half-open ranges covering the path, listed from `u` to
 *   `v`; the ranges from `u` up to the LCA run against position
 *   order (walk each right to left), the rest along it.
 *   `edges = true` leaves out the LCA, for values on edges stored
 *   at their lower endpoint. `climb(x, w, keep)` gives the (upward)
 *   ranges from `x` to its ancestor `w` (`w` itself if `keep`).
 *   Also `par` (-1 at the root), `depth`, `heavy` (-1 at leaves).
 * Precondition: `adj` holds both directions of the $n - 1$ edges of
 *   a tree on vertices `0..n-1`, $n >= 1$.
 * Time: $O(n)$ build (iterative); `path` returns $O(log n)$ ranges
 *   and `lca` takes $O(log n)$.
 * Usage:
 *   HLD h(adj); // rooted at 0; HLD h(adj, root)
 *   vector<ll> b(n); // ds/segtree (sum) over positions
 *   rep(v, 0, n) b[h.pos[v]] = a[v];
 *   SegTree st(b); // point update: st.set(h.pos[v], x)
 *   ll s = 0; // path sum
 *   for (auto [l, r] : h.path(u, v)) s += st.query(l, r);
 *   auto [l, r] = h.subtree(v); // st.query(l, r)
 *   // ds/lazy-segtree lz: add x to every vertex on a path
 *   for (auto [l, r] : h.path(u, v)) lz.apply(l, r, x);
 * Tweak: Edge values: put the value of edge `(par[v], v)` at
 *   `pos[v]` and pass `edges = true`; the subtree's edges are then
 *   `[pos[v] + 1, pos[v] + size[v])`. Any monoid works: the ranges
 *   go to whatever structure sits on the positions (Fenwick for
 *   sums, segment tree for max, lazy segment tree for path updates).
 *   Non-commutative `op`: with `w = lca(u, v)`, fold `path(u, w)`
 *   with `op` reversed, then `path(w, v, true)` as is.
 * Example:
 *   // Non-commutative: compose maps y -> b y + c along u -> v.
 *   struct M { ll b, c; };
 *   M then(M f, M g) { // f, then g
 *     return {g.b * f.b % MOD, (g.b * f.c + g.c) % MOD};
 *   }
 *   // In ds/segtree, nodes fold both ways (r: right to left):
 *   struct T { M f, r; };
 *   static constexpr T e{{1, 0}, {1, 0}};
 *   T op(T a, T b) { return {then(a.f, b.f), then(b.r, a.r)}; }
 *   // The leaf at pos[x] is {m, m} for the map m of vertex x.
 *   int w = h.lca(u, v);
 *   M f{1, 0};
 *   for (auto [l, r] : h.path(u, w))
 *     f = then(f, st.query(l, r).r);
 *   for (auto [l, r] : h.path(w, v, true))
 *     f = then(f, st.query(l, r).f);
 *   // (f.b * y + f.c) % MOD applies f_u first, f_v last
 * Notes: Order matters only for non-commutative operations.
 *   Positions grow downward along a chain, so `path(u, v)` is
 *   `climb(u, w)` (each range walked right to left) followed by
 *   `climb(v, w)` reversed (walked left to right), `w = lca(u, v)`.
 *   Hence the Example's split: `path(u, w)` is all upward,
 *   `path(w, v, true)` all downward without `w`. A second tree
 *   with `op`'s arguments swapped (a renamed copy of the struct)
 *   also works, but needs every update twice. Cues: "change the
 *   values on the path from `u` to `v`", "sum or max on a path
 *   with updates", "add to a subtree, query a path". Only subtree
 *   queries: plain DFS order suffices. Static values with path
 *   sums: prefix sums from the root and an LCA are simpler. A root
 *   path crosses at most $log_2 n$ light edges, so any path crosses
 *   at most $2 log_2 n$, and with a segment tree a path query costs
 *   $O(log^2 n)$. Subtree of `v` when the root moves to `r`: if `r`
 *   is a proper descendant of `v`, it is everything outside the
 *   subtree of `v`'s child `c` towards `r` (climb from `r`: `c` is
 *   the last head below `v`, or `heavy[v]` once on `v`'s chain).
 * Source: kactl HLD (CC0); maspy's library (CC0).
 * Status: stress-tested against explicit path enumeration
 */
#pragma once

struct HLD {
  vi par, heavy, head, pos, depth, size;
  HLD(const vector<vi>& adj, int root = 0) {
    int n = sz(adj), t = 0;
    par = heavy = vi(n, -1), head = pos = depth = vi(n);
    size = vi(n, 1);
    vi q{root}; // BFS order
    rep(i, 0, sz(q)) {
      int u = q[i];
      for (int v : adj[u])
        if (v != par[u])
          par[v] = u, depth[v] = depth[u] + 1, q.push_back(v);
    }
    for (int i = n; --i;) {
      int v = q[i], p = par[v];
      size[p] += size[v];
      if (heavy[p] < 0 || size[v] > size[heavy[p]])
        heavy[p] = v;
    }
    vi st{root};
    head[root] = root;
    while (sz(st)) {
      int v = st.back();
      st.pop_back();
      pos[v] = t++;
      for (int c : adj[v])
        if (c != par[v] && c != heavy[v])
          head[c] = c, st.push_back(c);
      if (heavy[v] >= 0)
        head[heavy[v]] = head[v], st.push_back(heavy[v]);
    }
  }
  vector<pii> path(int u, int v, bool edges = false) {
    int w = lca(u, v);
    vector<pii> a = climb(u, w, !edges), b = climb(v, w, false);
    a.insert(a.end(), b.rbegin(), b.rend());
    return a;
  }
  vector<pii> climb(int x, int w, bool keep) { // to ancestor w
    vector<pii> out;
    for (; head[x] != head[w]; x = par[head[x]])
      out.push_back({pos[head[x]], pos[x] + 1});
    if (x != w || keep)
      out.push_back({pos[w] + !keep, pos[x] + 1});
    return out;
  }
  pii subtree(int v) { return {pos[v], pos[v] + size[v]}; }
  int lca(int u, int v) {
    for (; head[u] != head[v]; v = par[head[v]])
      if (depth[head[u]] > depth[head[v]]) swap(u, v);
    return depth[u] < depth[v] ? u : v;
  }
};
