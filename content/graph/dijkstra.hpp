/**
 * Title: Dijkstra
 * Description: Shortest paths from a set of sources in a graph with
 *   non-negative weights. `dist[v]` is `LLONG_MAX` if `v` cannot be
 *   reached; `par[v]` is the vertex before `v` on a shortest path
 *   (-1 for sources and unreachable vertices). `path(t)` lists a
 *   shortest path from a source to `t`, empty if unreachable.
 * Precondition: `adj[u]` holds pairs `(v, w)` for edges `u -> v`;
 *   0-indexed, `w >= 0`, the sum of all weights fits in `ll`.
 * Time: $O((n + m) log m)$.
 * Usage:
 *   vector<vector<pair<int, ll>>> adj(n);
 *   adj[u].push_back({v, w});
 *   Dijkstra dj(adj, {s}); // or several sources
 *   ll d = dj.dist[t]; vi p = dj.path(t);
 * Example:
 *   // Dijkstra on states (vertex, coupons used): a coupon halves
 *   // one edge, at most k coupons in total.
 *   auto id = [&](int v, int c) { return v * (k + 1) + c; };
 *   vector<vector<pair<int, ll>>> g(n * (k + 1));
 *   for (auto [u, v, w] : edges) rep(c, 0, k + 1) {
 *     g[id(u, c)].push_back({id(v, c), w});
 *     if (c < k) g[id(u, c)].push_back({id(v, c + 1), w / 2});
 *   }
 *   Dijkstra dj(g, {id(s, 0)});
 *   ll best = LLONG_MAX;
 *   rep(c, 0, k + 1) best = min(best, dj.dist[id(t, c)]);
 * Notes: Think of it for "minimum total cost/time" with
 *   non-negative costs. Several starting points: pass them all as
 *   sources (multi-source), which also gives "distance to the
 *   nearest special vertex" for all vertices at once. Distances to
 *   one target: run
 *   on the reversed graph. Extra state (fuel left, parity, coupons
 *   used, time modulo $k$): run on pairs (vertex, state) as in the
 *   Example. Weights 0 and 1 only: 0-1 BFS with a deque (push front
 *   on weight 0, back on weight 1) is $O(n + m)$. Edge `(u, v, w)`
 *   lies on some shortest `s`-`t` walk iff
 *   $d_s (u) + w + d_t (v) = d_s (t)$ (run from `s` and, reversed,
 *   from `t`). Count shortest paths by DP over the tight edges
 *   ($d(u) + w = d(v)$) in topological order: increasing `dist` is
 *   one if all weights are positive, but 0-weight edges also need
 *   ties ordered along them (a 0-weight cycle gives infinitely many
 *   shortest walks). A negative edge breaks it silently: use
 *   Bellman-Ford or potentials (Johnson). Keep `dist` in `ll`.
 * Source: Classic; folklore contest implementation.
 * Status: stress-tested against Bellman-Ford style brute force
 */
#pragma once

struct Dijkstra {
  vector<ll> dist;
  vi par;
  Dijkstra(const auto& adj, const vi& srcs)
      : dist(sz(adj), LLONG_MAX), par(sz(adj), -1) {
    using P = pair<ll, int>;
    priority_queue<P, vector<P>, greater<P>> pq;
    for (int s : srcs) pq.push({dist[s] = 0, s});
    while (sz(pq)) {
      auto [d, u] = pq.top();
      pq.pop();
      if (d > dist[u]) continue;
      for (auto [v, w] : adj[u])
        if (d + w < dist[v])
          pq.push({dist[v] = d + w, v}), par[v] = u;
    }
  }
  vi path(int t) {
    vi p;
    if (dist[t] < LLONG_MAX)
      for (; t >= 0; t = par[t]) p.push_back(t);
    reverse(all(p));
    return p;
  }
};
