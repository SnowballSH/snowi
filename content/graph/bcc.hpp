/**
 * Title: Biconnected components
 * Description: Blocks of an undirected multigraph given as an edge
 *   list. `comp[e]` is the block of edge `e` (two edges share a
 *   block iff some simple cycle contains both); `bridge[e]` and
 *   `cut[v]` mark bridges and articulation points. `tree` is the
 *   block-cut forest: node `v < n` is vertex `v`, node `n + b` is
 *   block `b`, and each block is joined to all its vertices.
 * Precondition: Vertices 0-indexed. Self-loops are ignored
 *   (`comp[e] = -1`). Parallel edges form a cycle: they share a
 *   block and are never bridges. Isolated vertices are in no block.
 * Time: $O(n + m)$, iterative.
 * Usage:
 *   BCC b(n, edges); // vector<pii> edges
 *   b.comp[e]; b.bridge[e]; b.cut[v]; // nb blocks
 *   b.tree[n + b.comp[e]]; // the vertices of e's block
 * Tweak: 2-edge-connected components: join the endpoints of every
 *   non-bridge edge in a DSU. Vertex `v`'s blocks are its
 *   neighbours in `tree`; `v` is a cut vertex iff it has two or
 *   more.
 * Example:
 *   // 2-edge-connected components and the bridge tree.
 *   BCC b(n, edges);
 *   DSU d(n);
 *   rep(e, 0, m) if (!b.bridge[e])
 *     d.join(edges[e].first, edges[e].second);
 *   vector<vi> bt(n); // on DSU roots
 *   rep(e, 0, m) if (b.bridge[e]) {
 *     auto [u, v] = edges[e];
 *     int x = d.find(u), y = d.find(v);
 *     bt[x].push_back(y), bt[y].push_back(x);
 *   }
 * Notes: Use the block-cut tree for questions about removing one
 *   vertex: the vertices other than `u` and `v` whose removal
 *   separates them are exactly the vertex nodes strictly inside
 *   the tree path between `u` and `v` (all cut vertices; path
 *   queries with LCA). For removing one edge use bridges and the
 *   bridge tree (Example). A graph has a strongly connected
 *   orientation iff it is connected and bridgeless (orient along
 *   the DFS). Cactus graphs (each edge on at most one cycle):
 *   every block is a single edge or a cycle, so DP on the
 *   block-cut tree handles cycles like vertices. Blocks with 3 or
 *   more vertices are 2-vertex-connected: any two of their
 *   vertices lie on a common cycle. Parallel edges must be told
 *   apart by edge id, which is why the DFS skips only the parent
 *   edge, not the parent vertex.
 * Source: Hopcroft-Tarjan, iterative; kactl BiconnectedComponents
 *   (CC0).
 * Status: stress-tested against vertex/edge removal and cycle
 *   enumeration
 */
#pragma once

struct BCC {
  int nb = 0;
  vi comp;
  vector<bool> bridge, cut;
  vector<vi> tree;
  BCC(int n, const vector<pii>& edges)
      : comp(sz(edges), -1), bridge(sz(edges)), cut(n) {
    vector<vector<pii>> adj(n);
    rep(e, 0, sz(edges)) {
      auto [u, v] = edges[e];
      if (u != v)
        adj[u].push_back({v, e}), adj[v].push_back({u, e});
    }
    vi tin(n, -1), low(n), pe(n, -1), it(n), bid(n, -1);
    vi head, vs, cs; // block heads, vertex stack, DFS stack
    int t = 0;
    rep(r, 0, n) {
      if (tin[r] < 0) cs = {r}, tin[r] = low[r] = t++;
      while (sz(cs)) {
        int u = cs.back();
        if (it[u] < sz(adj[u])) {
          auto [v, e] = adj[u][it[u]++];
          if (tin[v] < 0) {
            tin[v] = low[v] = t++, pe[v] = e;
            cs.push_back(v), vs.push_back(v);
          } else if (e != pe[u]) low[u] = min(low[u], tin[v]);
          continue;
        }
        cs.pop_back();
        if (cs.empty()) break;
        int p = cs.back();
        low[p] = min(low[p], low[u]);
        if (low[u] < tin[p]) continue;
        // p cuts off u's subtree: pop it as a block with head p
        bridge[pe[u]] = low[u] > tin[p];
        int v;
        do v = vs.back(), vs.pop_back(), bid[v] = nb;
        while (v != u);
        head.push_back(p), nb++;
      }
    }
    tree.resize(n + nb);
    auto link = [&](int v, int b) {
      tree[v].push_back(n + b), tree[n + b].push_back(v);
    };
    rep(v, 0, n)
      if (bid[v] >= 0) link(v, bid[v]);
    rep(b, 0, nb) link(head[b], b);
    rep(v, 0, n) cut[v] = sz(tree[v]) > 1;
    // An edge belongs to the block of its deeper endpoint.
    rep(e, 0, sz(edges)) {
      auto [u, v] = edges[e];
      if (u != v) comp[e] = bid[tin[u] > tin[v] ? u : v];
    }
  }
};
