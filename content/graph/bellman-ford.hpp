/**
 * Title: Bellman-Ford
 * Description: Shortest paths from `s` when weights may be negative.
 *   `dist[v]` is `LLONG_MAX` if `v` cannot be reached and
 *   `LLONG_MIN` if paths to `v` can be arbitrarily short (`v` is
 *   reachable from a negative cycle that `s` reaches). For finite
 *   `dist[v]`, `par[v]` is the vertex before `v` on a shortest path
 *   (-1 for `s`).
 * Precondition: Edges `(a, b, w)` are directed `a -> b`; vertices
 *   0-indexed. `ll` must hold $W n m$, where $W$ is the largest
 *   $|w|$.
 * Time: $O(n m)$; stops as soon as a round changes nothing.
 * Usage:
 *   vector<tuple<int, int, ll>> edges; // a -> b, weight w
 *   BellmanFord bf(n, edges, s);
 *   if (bf.dist[v] == LLONG_MIN) ... // no shortest path
 * Tweak: Print a negative cycle (any one if `dist` starts all 0):
 *   relax plainly for $n$ rounds (no `LLONG_MIN` marking); if round
 *   $n$ still relaxes some `b`, follow `par` from `b` $n$ times,
 *   then once around. Real weights (arbitrage): `double` for `ll`
 *   in `dist` and `edges` (the sentinels still work), and relax
 *   only if `da + w < dist[b] - 1e-9`, or rounding fakes negative
 *   cycles.
 * Example:
 *   // Difference constraints x[b] - x[a] <= w: edge a -> b
 *   // with weight w, plus the source n. Feasible iff no -inf;
 *   // x = dist is a solution (each x[v] <= 0 and maximal).
 *   rep(v, 0, n) edges.push_back({n, v, 0});
 *   BellmanFord bf(n + 1, edges, n);
 *   bool ok = count(all(bf.dist), LLONG_MIN) == 0;
 * Notes: Use it when weights are negative: arbitrage (weights
 *   $-log(r)$ in `double`, see Tweak), difference constraints
 *   (Example), "is there a cycle that gains money", or potentials
 *   for Johnson's all-pairs algorithm (then Dijkstra on $w(u, v) +
 *   d(u) - d(v) >= 0$). An undirected negative edge is already a
 *   negative cycle. SPFA (queue-based) has
 *   the same worst case and is often targeted by anti-tests.
 * Source: kactl BellmanFord (CC0), simplified (no edge sorting).
 * Status: stress-tested against walk-length dynamic programming
 */
#pragma once

struct BellmanFord {
  vector<ll> dist;
  vi par;
  BellmanFord(int n, const vector<tuple<int, int, ll>>& edges,
              int s)
      : dist(n, LLONG_MAX), par(n, -1) {
    dist[s] = 0;
    // Rounds >= n - 1 only mark: any change there means -inf.
    for (int i = 0, ch = 1; ch && i < 2 * n; i++) {
      ch = 0;
      for (auto [a, b, w] : edges) {
        auto da = dist[a];
        if (da == LLONG_MAX || dist[b] == LLONG_MIN) continue;
        if (da == LLONG_MIN || da + w < dist[b]) {
          dist[b] = i < n - 1 ? da + w : LLONG_MIN;
          par[b] = a, ch = 1;
        }
      }
    }
  }
};
