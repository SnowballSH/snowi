/**
 * Title: Min-cost max flow
 * Description: Sends as much flow as possible from `s` to `t` (at
 *   most `limit`), and among maximum flows one of minimum total
 *   cost. Returns `{flow, cost}`. `flow(id)` is the flow on an edge.
 *   Successive shortest paths: Bellman-Ford once for the
 *   potentials, then Dijkstra on reduced costs, stopping at `t`.
 * Precondition: Vertices 0-indexed, `s != t`, `cap >= 0`. Costs
 *   may be negative but no cycle of positive capacity may have
 *   negative total cost. The total cost and $8 n max|"cost"|$
 *   stay below $2^63$.
 * Time: $O(F m log m)$ for $F$ augmenting paths (at most the flow
 *   value), plus $O(n m)$ when costs are negative.
 * Usage:
 *   MinCostFlow mcf(n);
 *   int id = mcf.add_edge(u, v, cap, cost);
 *   auto [f, c] = mcf.max_flow(s, t); // or max_flow(s, t, k)
 *   mcf.flow(id);
 * Tweak: Maximum cost: negate costs (no positive cycles then).
 *   Cheapest flow of any size: put
 *   `if (pot[t] - pot[s] >= 0) break;` right after the `pot`
 *   update. That difference is the cost of the path about to be
 *   used (in the loop condition it would be stale). Exactly `k`
 *   units: call with `limit = k` and check that `f == k`.
 * Example:
 *   // k edge-disjoint paths 0 -> n-1 of minimum total length:
 *   MinCostFlow mcf(n);
 *   for (auto [u, v, w] : edges) mcf.add_edge(u, v, 1, w);
 *   auto [f, len] = mcf.max_flow(0, n - 1, k);
 *   if (f < k) ... // fewer than k disjoint paths exist
 * Notes: Model: capacities say how often something can be used,
 *   costs what each use costs. Cues: "minimum total cost to
 *   transport", "assign with capacities and prices", "choose
 *   $k$ disjoint paths or cycles", or a weighted matching too
 *   sparse for Hungarian. A convex cost per unit on one edge
 *   (1st unit costs 1, 2nd costs 3, ...) becomes parallel unit
 *   edges with increasing costs. Each augmenting path costs a
 *   Dijkstra: $F = 1000$ paths over $2.4 dot 10^4$ arcs take
 *   about 0.4 s in C++ (1 to 2 s in PyPy), growing like
 *   $F m log m$. Dijkstra stops once `t` is settled; unsettled
 *   vertices gain only $d_t$ in potential, which keeps reduced
 *   costs $>= 0$. Negative cycles are not supported. For a
 *   cheapest circulation with them, saturate every negative
 *   edge first (its residual reverse costs $> 0$) and route the
 *   resulting excesses from a super source to a super sink. On
 *   a path or tree the transport cost has a closed form: every
 *   edge carries $|$imbalance of one side$|$.
 * Source: Classic successive shortest paths with Johnson
 *   potentials; kactl MinCostMaxFlow (CC0); early stop at `t` as
 *   in the AtCoder Library (CC0).
 * Status: stress-tested against exhaustive flows and Bellman-Ford SSP
 */
#pragma once

struct MinCostFlow {
  struct E {
    int to, rev;
    ll cap, cost;
  };
  vector<vector<E>> g;
  vector<pii> pos; // edge id -> (u, index in g[u])
  vector<ll> pot;
  MinCostFlow(int n) : g(n), pot(n) {}
  int add_edge(int u, int v, ll cap, ll cost) {
    pos.push_back({u, sz(g[u])});
    g[u].push_back({v, sz(g[v]) + (u == v), cap, cost});
    g[v].push_back({u, sz(g[u]) - 1, 0, -cost});
    return sz(pos) - 1;
  }
  ll flow(int id) {
    E& x = g[pos[id].first][pos[id].second];
    return g[x.to][x.rev].cap;
  }
  pair<ll, ll> max_flow(int s, int t, ll limit = LLONG_MAX) {
    int n = sz(g);
    fill(all(pot), 0); // Bellman-Ford: valid potentials
    for (int i = 0, ch = 1; ch && i < n; i++) {
      ch = 0;
      rep(u, 0, n)
        for (E& x : g[u])
          if (x.cap && pot[u] + x.cost < pot[x.to])
            pot[x.to] = pot[u] + x.cost, ch = 1;
    }
    ll f = 0, c = 0;
    vector<ll> d(n);
    vector<pii> par(n); // (vertex, edge index) before v
    while (f < limit) {
      fill(all(d), LLONG_MAX), d[s] = 0;
      priority_queue<pair<ll, int>> pq; // (-dist, v)
      pq.push({0, s});
      while (sz(pq)) {
        auto [du, u] = pq.top();
        pq.pop();
        if (-du > d[u]) continue;
        if (u == t) break;
        rep(i, 0, sz(g[u])) {
          E& x = g[u][i];
          ll nd = d[u] + x.cost + pot[u] - pot[x.to];
          if (x.cap && nd < d[x.to])
            d[x.to] = nd, par[x.to] = {u, i},
            pq.push({-nd, x.to});
        }
      }
      if (d[t] == LLONG_MAX) break;
      rep(v, 0, n) pot[v] += min(d[v], d[t]);
      ll push = limit - f;
      for (int v = t; v != s; v = par[v].first)
        push = min(push, g[par[v].first][par[v].second].cap);
      for (int v = t; v != s; v = par[v].first) {
        E& x = g[par[v].first][par[v].second];
        x.cap -= push, g[v][x.rev].cap += push;
      }
      f += push, c += push * (pot[t] - pot[s]);
    }
    return {f, c};
  }
};
