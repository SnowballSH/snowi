/**
 * Title: Topological sort
 * Description: Kahn's algorithm: orders the vertices so that every
 *   edge `u -> v` has `u` before `v`. Returns an empty vector if the
 *   graph has a cycle.
 * Precondition: `adj[u]` lists the heads of the edges leaving `u`;
 *   vertices 0-indexed. Parallel edges are fine; a self-loop is a
 *   cycle.
 * Time: $O(n + m)$.
 * Usage:
 *   vi order = topo_sort(adj); // empty: cycle (when n > 0)
 * Tweak: Lexicographically smallest order: replace the scan of
 *   `order` with a `priority_queue<int, vi, greater<int>>`,
 *   $O(n log n + m)$. To keep the acyclic part, drop the final
 *   check: the result then holds exactly the vertices not
 *   reachable from a cycle.
 * Example:
 *   // Longest path in a DAG with weights wadj[u] = {(v, w)}:
 *   vector<ll> best(n, 0); // longest path ending at v
 *   for (int u : topo_sort(adj))
 *     for (auto [v, w] : wadj[u])
 *       best[v] = max(best[v], best[u] + w);
 * Notes: Use it for "order tasks with prerequisites", for DP on a
 *   DAG (longest path, number of paths, reachability bitsets), and
 *   to detect cycles. The order is unique iff the queue never holds
 *   two vertices at once, i.e. consecutive vertices are joined by an
 *   edge. "Smallest label as early as possible" (vertex 0 first if
 *   it can be, then 1, ...) is not the lexicographically smallest
 *   order: reverse all edges, take the lexicographically largest
 *   order with a max-heap, and reverse the result. Constraints
 *   "`a` before `b`" from comparisons (sorting by partial
 *   information, dictionary order of an alien alphabet) are a
 *   topological sort.
 * Source: kactl TopoSort (CC0).
 * Status: stress-tested against brute force
 */
#pragma once

vi topo_sort(const vector<vi>& adj) {
  int n = sz(adj);
  vi indeg(n), order;
  for (auto& out : adj)
    for (int v : out) indeg[v]++;
  rep(u, 0, n)
    if (!indeg[u]) order.push_back(u);
  rep(i, 0, sz(order))
    for (int v : adj[order[i]])
      if (!--indeg[v]) order.push_back(v);
  if (sz(order) < n) order.clear();
  return order;
}
