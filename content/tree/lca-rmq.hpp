/**
 * Title: LCA by Euler tour and RMQ
 * Description: Lowest common ancestor in $O(1)$ after an
 *   $O(n log n)$ build. With a DFS preorder `tin`, entry `i` of the
 *   tour holds `tin` of the parent of the vertex entered at time
 *   `i + 1`; for `u != v`, the LCA is the vertex whose `tin` is the
 *   minimum of entries `[min(tin), max(tin))`. Like tree/lca it
 *   has `depth`, `tin` and `order` (`order[tin[v]] = v`); it also
 *   has `par`, -1 at the root (tree/lca has `up[0]` instead,
 *   which maps the root to itself).
 * Precondition: `adj` holds both directions of the $n - 1$ edges of
 *   a tree on vertices `0..n-1`, $n >= 1$.
 * Time: $O(n log n)$ time and memory to build (iterative), $O(1)$
 *   per query.
 * Usage:
 *   RmqLCA lc(adj); // rooted at 0; RmqLCA lc(adj, root)
 *   lc.lca(u, v); lc.dist(u, v); lc.depth[v];
 * Tweak: Weighted distance: keep the weighted depth `wd` beside
 *   `depth` and return `wd[u] + wd[v] - 2 * wd[lca(u, v)]`.
 * Example:
 *   // Steiner tree: edges needed to connect a set of vertices
 *   // = half the closed tour through them in preorder.
 *   sort(all(s), [&](int a, int b) {
 *     return lc.tin[a] < lc.tin[b];
 *   });
 *   ll edges = 0;
 *   rep(i, 0, sz(s)) edges += lc.dist(s[i], s[(i + 1) % sz(s)]);
 *   edges /= 2;
 * Notes: Prefer it over binary lifting when queries dominate
 *   ($10^6$ LCA or distance queries) or inside another loop (virtual
 *   trees, sorting by distance). The tour has $n - 1$ entries, half
 *   the classic $2n - 1$ tour of depths. Memory is $n log_2 n$ ints:
 *   80 MB at $n = 10^6$. It answers only LCA and distance; for
 *   k-th ancestors and `jump` use tree/lca.
 * Source: kactl LCA (CC0), made iterative.
 * Status: stress-tested against naive parent walking
 */
#pragma once
#include "../ds/sparse-table.hpp"

struct RmqLCA {
  vi par, depth, tin, order;
  SparseTable<int> rmq;
  RmqLCA(const vector<vi>& adj, int root = 0)
      : par(sz(adj), -1), depth(sz(adj)), tin(sz(adj)),
        rmq(tour(adj, root)) {}
  vi tour(const vector<vi>& adj, int root) {
    vi st{root}, top;
    while (sz(st)) {
      int u = st.back();
      st.pop_back();
      if (u != root) top.push_back(tin[par[u]]);
      tin[u] = sz(order), order.push_back(u);
      for (int v : adj[u])
        if (v != par[u])
          par[v] = u, depth[v] = depth[u] + 1, st.push_back(v);
    }
    return top;
  }
  int lca(int u, int v) {
    if (u == v) return u;
    auto [a, b] = minmax(tin[u], tin[v]);
    return order[rmq.query(a, b)];
  }
  int dist(int u, int v) {
    return depth[u] + depth[v] - 2 * depth[lca(u, v)];
  }
};
