/**
 * Title: General matching
 * Description: Maximum matching in any undirected graph (Edmonds'
 *   blossom algorithm). Returns the size `k` and `mate`, where
 *   `mate[v]` is the partner of `v`, or -1 if `v` is unmatched.
 * Precondition: 0-indexed adjacency lists, each edge in both
 *   `g[u]` and `g[v]`. Self-loops and parallel edges are fine.
 * Time: $O(n^3)$, at most one search per vertex; far faster in
 *   practice: about 0.3 s for a random graph with $n = 10^5$,
 *   $m = 2 dot 10^5$, under 0.05 s for paths and random trees.
 * Usage:
 *   vector<vi> g(n);
 *   g[u].push_back(v), g[v].push_back(u);
 *   auto [k, mate] = general_matching(g);
 * Tweak: Bipartite graphs: `flow/hopcroft-karp` is faster.
 *   Weighted matching in a general graph needs the much longer
 *   weighted blossom algorithm; for $n <= 20$ use a bitmask DP.
 * Example:
 *   // Players alternately move a token to an unvisited
 *   // neighbour; whoever cannot move loses. The first player
 *   // wins from s iff every maximum matching covers s.
 *   auto h = g;
 *   h[s].clear();
 *   for (auto& adj : h) erase(adj, s);
 *   bool first_wins = general_matching(h).first < k;
 * Notes: Do not make a general graph bipartite by doubling the
 *   vertices: a triangle then gets 3 edges instead of 1 pair.
 *   Tutte-Berge: $2k = min_U (n + |U| - "odd"(G - U))$, where
 *   $"odd"(G - U)$ counts the components of odd size after
 *   deleting $U$; a perfect matching exists iff
 *   $"odd"(G - U) <= |U|$ for all $U$. The smallest edge cover
 *   (no isolated vertices) has $n - k$ edges, but Konig fails:
 *   the largest independent set is NP-hard here. Cues: "pair up"
 *   when any two items may pair (not two classes), dominoes on a
 *   graph with odd cycles, the token game above. Each vertex
 *   starts at most one search: a vertex with no augmenting path
 *   now never gets one later, and the whole tree of a failed
 *   search can be dropped for good (Edmonds). Hand-written
 *   versions that rebuild it are quadratic on a long path with
 *   many leaves at one end.
 * Source: Edmonds (1965); union-find blossom bases as in the
 *   common BFS formulation (Gabow 1976).
 * Status: stress-tested against bitmask DP and Tutte matrix rank
 */
#pragma once

pair<int, vi> general_matching(const vector<vi>& g) {
  int n = sz(g), k = 0, stamp = 0;
  vi mate(n, -1), pre(n), dsu(n), lab(n), seen(n), q;
  iota(all(dsu), 0); // lab: 0 unseen, 1 even, 2 odd, 3 dead
  auto base = [&](int x) {
    while (dsu[x] != x) x = dsu[x] = dsu[dsu[x]];
    return x;
  };
  auto lca = [&](int x, int y) {
    for (++stamp;; swap(x, y))
      if (x >= 0) {
        if (seen[x = base(x)] == stamp) return x;
        seen[x] = stamp, x = mate[x] < 0 ? -1 : pre[mate[x]];
      }
  };
  auto shrink = [&](int x, int y, int b) {
    for (; base(x) != b; x = pre[y]) {
      pre[x] = y, y = mate[x];
      if (lab[y] == 2) lab[y] = 1, q.push_back(y);
      if (dsu[x] == x) dsu[x] = b;
      if (dsu[y] == y) dsu[y] = b;
    }
  };
  auto augment = [&](int root) {
    for (int x : q) lab[x] = 0, dsu[x] = x; // undo last search
    q = {root}, lab[root] = 1;
    rep(i, 0, sz(q)) {
      int u = q[i];
      if (lab[u] == 1)
        for (int v : g[u]) {
          if (lab[v] == 1 && base(u) != base(v)) {
            int b = lca(u, v);
            shrink(u, v, b), shrink(v, u, b);
          } else if (!lab[v]) {
            lab[v] = 2, pre[v] = u, q.push_back(v);
            if (mate[v] < 0) {
              for (int x = v; x >= 0;) {
                int y = mate[pre[x]];
                mate[x] = pre[x], mate[pre[x]] = x, x = y;
              }
              return true;
            }
            lab[mate[v]] = 1, q.push_back(mate[v]);
          }
        }
    }
    for (int x : q) lab[x] = 3; // no path through them, ever
    q.clear();
    return false;
  };
  rep(v, 0, n)
    if (mate[v] < 0) k += augment(v);
  return {k, mate};
}
