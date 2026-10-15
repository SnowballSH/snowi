/**
 * Title: Offline dynamic connectivity
 * Description: Edges are added and removed in time order, and
 *   `query()` marks a moment to ask about the graph. `solve` then
 *   calls `answer(i)` for the queries `i = 0, 1, ...` in order,
 *   each time with the rollback DSU `d` holding exactly the edges
 *   present at query `i`. An edge present during queries `[l, r)`
 *   is stored in $O(log q)$ nodes of a segment tree over query
 *   times; a DFS joins a node's edges on entry and rolls them back
 *   on exit.
 * Precondition: Vertices `0..n-1`. `remove(u, v)` only while that
 *   edge is present; parallel edges and loops are fine. Call
 *   `solve` once, after all operations.
 * Time: $O(m log q log n)$ for $m$ added edges and $q$ queries,
 *   plus $O(m log m)$ for the map.
 * Usage:
 *   DynConn dc(n);
 *   dc.add(u, v); dc.remove(u, v); // in time order
 *   int i = dc.query(); // ask about the graph at this moment
 *   dc.solve([&](int i) {
 *     same[i] = dc.d.find(a[i]) == dc.d.find(b[i]);
 *     comps[i] = n - dc.d.time() / 2;
 *   });
 * Tweak: Component size: `dc.d.size(v)`. Bipartiteness, sums or
 *   counts per component: give the DSU that state, record it in
 *   `join` and restore it in `rollback`. Any structure with
 *   "insert" and "undo the last inserts" works instead of the DSU
 *   (XOR basis, knapsack array: save a copy per node). Prune: once
 *   a node's edges decide all queries below it (an odd cycle),
 *   answer them there and do not descend.
 * Example:
 *   // Ops "+ u v", "- u v" (1-indexed) and "?": print the number
 *   // of components at each "?".
 *   DynConn dc(n);
 *   rep(t, 0, k) {
 *     char op;
 *     int u = 0, v = 0;
 *     cin >> op;
 *     if (op != '?') cin >> u >> v, u--, v--;
 *     if (op == '+') dc.add(u, v);
 *     if (op == '-') dc.remove(u, v);
 *     if (op == '?') dc.query();
 *   }
 *   dc.solve([&](int) { cout << n - dc.d.time() / 2 << '\n'; });
 * Notes: Think of it when edges (or items) appear and disappear,
 *   everything is known in advance, and the structure can insert
 *   and undo but not delete: connectivity, bipartiteness (DSU with
 *   parity), sums or counts per component, XOR basis, knapsack.
 *   An edge is alive on the half-open range from its addition to
 *   its removal: a query between them sees it, one after the
 *   removal does not. Only deletions: reverse time, so that they
 *   become insertions into a plain DSU. Forced online (queries
 *   decoded with the last answer): this does not apply. Memory is
 *   $O(m log q)$ stored edges.
 * Source: Segment tree over time; kactl RollbackUF (CC0); UCup
 *   S4-S5 Training Packet 14.
 * Status: stress-tested against brute force
 */
#pragma once
#include "../ds/dsu-rollback.hpp"

struct DynConn {
  int q = 0; // queries so far
  RollbackDSU d;
  map<pii, vi> open; // edge -> start times of present copies
  vector<array<int, 4>> ev; // l, r, u, v: alive for [l, r)
  DynConn(int n) : d(n) {}
  void add(int u, int v) { open[minmax(u, v)].push_back(q); }
  void remove(int u, int v) {
    vi& s = open[minmax(u, v)];
    ev.push_back({s.back(), q, u, v}), s.pop_back();
  }
  int query() { return q++; }
  void solve(auto answer) {
    for (auto& [e, s] : open)
      for (int l : s) ev.push_back({l, q, e.first, e.second});
    int S = 1;
    while (S < q) S *= 2;
    vector<vector<pii>> at(2 * S);
    for (auto [l, r, u, v] : ev)
      for (l += S, r += S; l < r; l /= 2, r /= 2) {
        if (l & 1) at[l++].push_back({u, v});
        if (r & 1) at[--r].push_back({u, v});
      }
    auto dfs = [&](auto&& self, int x) -> void {
      int t = d.time();
      for (auto [u, v] : at[x]) d.join(u, v);
      if (x < S) self(self, 2 * x), self(self, 2 * x + 1);
      else if (x - S < q) answer(x - S);
      d.rollback(t);
    };
    dfs(dfs, 1);
  }
};
