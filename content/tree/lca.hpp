/**
 * Title: LCA by binary lifting
 * Description: Lowest common ancestor, distance in edges, k-th
 *   ancestor, and `jump(u, v, k)`: the vertex `k` edges from `u` on
 *   the path to `v`. Also gives `depth`, the parent `up[0][v]` (the
 *   root's is itself) and a DFS preorder: `order[tin[v]] = v`.
 * Precondition: `adj` holds both directions of the $n - 1$ edges of
 *   a tree on vertices `0..n-1`, $n >= 1$. `k >= 0`.
 * Time: $O(n log n)$ time and memory to build (iterative), then
 *   $O(log n)$ per query.
 * Usage:
 *   LCA lc(adj); // rooted at 0; LCA lc(adj, root) otherwise
 *   lc.lca(u, v); lc.dist(u, v); lc.depth[v];
 *   lc.kth_ancestor(u, k); // -1 if k > depth[u]
 *   lc.jump(u, v, k); // -1 if k > dist(u, v); k = 0 gives u
 * Example:
 *   // Max edge weight on the path u - v (MST bottleneck queries).
 *   // w[v] = weight of the edge from v to its parent.
 *   vector<vi> mx(sz(lc.up), w); // max of 2^k edges above v
 *   rep(k, 1, sz(mx)) rep(v, 0, n) mx[k][v] =
 *       max(mx[k - 1][v], mx[k - 1][lc.up[k - 1][v]]);
 *   auto climb = [&](int v, int d) { // max of d edges above v
 *     int r = INT_MIN;
 *     rep(k, 0, sz(mx)) if (d >> k & 1)
 *       r = max(r, mx[k][v]), v = lc.up[k][v];
 *     return r;
 *   };
 *   int a = lc.lca(u, v), d = lc.depth[a];
 *   int best = max(climb(u, lc.depth[u] - d),
 *                  climb(v, lc.depth[v] - d));
 * Notes: Cues: "distance between two vertices", "k-th vertex on the
 *   path", "is `x` on the path `a - b`" (iff `dist(a, x) + dist(x,
 *   b) == dist(a, b)`), "middle of a path" (`jump(u, v, dist(u, v)
 *   / 2)`; on a diameter `u - v` it is a center). `x` is an
 *   ancestor of `v` iff `depth[x] <= depth[v]` and `kth_ancestor(v,
 *   depth[v] - depth[x]) == x`, or in $O(1)$ with subtree sizes:
 *   `tin[x] <= tin[v] < tin[x] + size[x]` (sizes: add each
 *   `order[i]` into its parent for `i` from `n - 1` down to 1).
 *   Paths `a - b` and `c - d` meet iff `lca(a, b)` lies on `c - d`
 *   or `lca(c, d)` lies on `a - b`. The same lifting finds the
 *   highest ancestor with a property that is monotone along the
 *   root path (value below `x`, Kruskal reconstruction trees): for
 *   `k` from the top, jump while `up[k][v]` still satisfies it.
 *   Memory is $n log_2 n$ ints (80 MB at $n = 10^6$): for $10^6$
 *   queries use tree/lca-rmq ($O(1)$ query); when memory is tight,
 *   `HLD::lca` ($O(n)$ memory). For a forest, join one vertex of
 *   each tree to a new root `n` and test reachability separately.
 * Source: Folklore binary lifting; cp-algorithms.
 * Status: stress-tested against naive parent walking
 */
#pragma once

struct LCA {
  vi depth, tin, order;
  vector<vi> up; // up[k][v]: 2^k-th ancestor or the root
  LCA(const vector<vi>& adj, int root = 0)
      : depth(sz(adj)), tin(sz(adj)) {
    int n = sz(adj);
    vi par(n, root), st{root};
    while (sz(st)) {
      int u = st.back();
      st.pop_back();
      tin[u] = sz(order), order.push_back(u);
      for (int v : adj[u])
        if (v != par[u])
          par[v] = u, depth[v] = depth[u] + 1, st.push_back(v);
    }
    up = {par};
    for (int k = 1; (1 << k) < n; k++) {
      up.emplace_back(n);
      rep(v, 0, n) up[k][v] = up[k - 1][up[k - 1][v]];
    }
  }
  int kth_ancestor(int u, int k) {
    if (k > depth[u]) return -1;
    rep(j, 0, sz(up))
      if (k >> j & 1) u = up[j][u];
    return u;
  }
  int lca(int u, int v) {
    if (depth[u] < depth[v]) swap(u, v);
    u = kth_ancestor(u, depth[u] - depth[v]);
    for (int j = sz(up); j--;)
      if (up[j][u] != up[j][v]) u = up[j][u], v = up[j][v];
    return u == v ? u : up[0][u];
  }
  int dist(int u, int v) {
    return depth[u] + depth[v] - 2 * depth[lca(u, v)];
  }
  int jump(int u, int v, int k) {
    int a = lca(u, v), du = depth[u] - depth[a];
    int d = du + depth[v] - depth[a];
    if (k > d) return -1;
    return k <= du ? kth_ancestor(u, k)
                   : kth_ancestor(v, d - k);
  }
};
