/**
 * Title: Strongly connected components
 * Description: Tarjan's algorithm, iterative. Returns `comp` and
 *   the number of components `k`: `comp[u] == comp[v]` iff `u` and
 *   `v` reach each other. Components are numbered in reverse
 *   topological order: an edge `u -> v` between different
 *   components has `comp[u] > comp[v]`, so component 0 is a sink and
 *   `k - 1, ..., 0` is a topological order of the condensation.
 * Precondition: `adj[u]` lists the heads of the edges leaving `u`;
 *   vertices 0-indexed. Self-loops and parallel edges are fine.
 * Time: $O(n + m)$; no recursion, fine for $n = 10^6$.
 * Usage:
 *   auto [comp, k] = scc(adj);
 * Tweak: Topological numbering (sources first): `k - 1 - comp[v]`.
 *   Vertex lists: push each `v` into `g[comp[v]]` of a
 *   `vector<vi> g(k)`.
 * Example:
 *   // Condensation DAG, edges deduplicated.
 *   auto [comp, k] = scc(adj);
 *   vector<vi> dag(k);
 *   rep(u, 0, n) for (int v : adj[u])
 *     if (comp[u] != comp[v]) dag[comp[u]].push_back(comp[v]);
 *   for (vi& out : dag) {
 *     sort(all(out));
 *     out.erase(unique(all(out)), out.end());
 *   }
 *   // DP over it, sinks first: for (c = 0; c < k; c++) ...
 * Notes: Think SCC when the statement says "mutually reachable",
 *   "can return to", "cycles can be traversed many times" (collapse
 *   each component, then DP on the DAG: e.g. maximum total value
 *   collected on a walk), or for 2-SAT. Minimum number of edges to
 *   add to make the graph strongly connected: $max(s, t)$ for $s$
 *   sources and $t$ sinks of the condensation, 0 if $k = 1$. A
 *   vertex reachable from all others exists iff the condensation
 *   has exactly one sink. In the DP, process components in
 *   increasing id for "what can be reached from here" and in
 *   decreasing id for "what reaches here".
 * Source: kactl SCC (CC0), made iterative.
 * Status: stress-tested against brute-force reachability
 */
#pragma once

pair<vi, int> scc(const vector<vi>& adj) {
  int n = sz(adj), t = 0, k = 0;
  vi comp(n, -1), tin(n, -1), low(n), it(n), st, cs;
  auto visit = [&](int v) {
    cs.push_back(v), st.push_back(v), tin[v] = low[v] = t++;
  };
  rep(r, 0, n) {
    if (tin[r] < 0) visit(r);
    while (sz(cs)) {
      int u = cs.back();
      if (it[u] < sz(adj[u])) {
        int v = adj[u][it[u]++];
        if (tin[v] < 0) visit(v);
        else if (comp[v] < 0) low[u] = min(low[u], tin[v]);
        continue;
      }
      cs.pop_back();
      if (sz(cs)) low[cs.back()] = min(low[cs.back()], low[u]);
      if (low[u] < tin[u]) continue;
      int v;
      do v = st.back(), st.pop_back(), comp[v] = k;
      while (v != u);
      k++;
    }
  }
  return {comp, k};
}
