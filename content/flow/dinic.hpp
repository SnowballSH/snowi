/**
 * Title: Dinic max flow
 * Description: Maximum flow from `s` to `t`. Afterwards `flow(id)`
 *   is the flow on edge `id` (negative if it runs backwards on an
 *   undirected edge), and `source_side(v)` tells whether `v` is on
 *   the source side of a minimum cut: reachable from `s` in the
 *   residual graph.
 * Precondition: Vertices 0-indexed, `s != t`, capacities `>= 0`.
 *   Each `cap + rev_cap` and the total capacity out of `s` stay
 *   below $2^63$.
 * Time: $O(n^2 m)$, far faster in practice. Unit capacities:
 *   $O(m sqrt(m))$, and $O(m sqrt(n))$ if every vertex but $s$,
 *   $t$ has a single edge in or a single edge out (matching,
 *   split vertices).
 * Usage:
 *   Dinic d(n);
 *   int id = d.add_edge(u, v, cap); // directed u -> v
 *   d.add_edge(a, b, c, c); // undirected, capacity c
 *   ll f = d.max_flow(s, t);
 *   d.flow(id); d.source_side(v);
 * Tweak: A second `max_flow` call keeps the flow already sent and
 *   adds to it. To start over on the same graph, set `x.cap =
 *   x.orig` for every edge `x` in `g`.
 * Example:
 *   // Fewest undirected roads to block so that 0 cannot reach
 *   // n - 1, and which ones (CSES Police Chase).
 *   Dinic d(n);
 *   for (auto [a, b] : roads) d.add_edge(a, b, 1, 1);
 *   cout << d.max_flow(0, n - 1) << "\n";
 *   for (auto [a, b] : roads)
 *     if (d.source_side(a) != d.source_side(b))
 *       cout << a + 1 << " " << b + 1 << "\n";
 * Notes: Read a minimum cut as the edges from source-side to
 *   sink-side vertices, not as "saturated edges": flow paths can
 *   saturate edges anywhere. An undirected edge is one call with
 *   `rev_cap = c`, not two calls. An "infinite" capacity such as
 *   $10^18$ must not overflow when a few are added. The recursion
 *   is as deep as the distance from `s` to `t`. Unit capacities
 *   without parallel edges also give $O(m n^(2/3))$; parallel unit
 *   edges can force $n - 2$ phases. Fast on matching
 *   and layered networks; on big grids with varied capacities it
 *   can take hundreds of phases. Think of max flow for "can every
 *   demand be met", "pair things up, each usable $k$ times",
 *   "fewest edges or vertices to delete to separate", "pick a
 *   subset where taking $x$ forces $y$": see the modelling note.
 *   With $n$ up to $10^6$ the intended solution is usually a
 *   greedy or a formula proved by the cut condition, not a flow.
 * Source: kactl Dinic (CC0); edge ids as in the AtCoder Library.
 * Status: stress-tested against cut enumeration and Edmonds-Karp
 */
#pragma once

struct Dinic {
  struct E {
    int to, rev;
    ll cap, orig;
  };
  vector<vector<E>> g;
  vector<pii> pos; // edge id -> (u, index in g[u])
  vi lvl, it;
  Dinic(int n) : g(n), lvl(n), it(n) {}
  int add_edge(int u, int v, ll cap, ll rev_cap = 0) {
    pos.push_back({u, sz(g[u])});
    g[u].push_back({v, sz(g[v]), cap, cap});
    g[v].push_back({u, sz(g[u]) - 1, rev_cap, rev_cap});
    return sz(pos) - 1;
  }
  ll flow(int id) {
    E& x = g[pos[id].first][pos[id].second];
    return x.orig - x.cap;
  }
  bool source_side(int v) { return lvl[v] >= 0; }
  bool bfs(int s, int t) {
    fill(all(lvl), -1), lvl[s] = 0;
    vi q{s};
    rep(i, 0, sz(q))
      for (E& x : g[q[i]]) {
        if (!x.cap || lvl[x.to] >= 0) continue;
        lvl[x.to] = lvl[q[i]] + 1, q.push_back(x.to);
        if (x.to == t) return true;
      }
    return false;
  }
  ll dfs(int u, int t, ll f) {
    if (u == t) return f;
    ll res = 0;
    for (int& i = it[u]; i < sz(g[u]); i++) {
      E& x = g[u][i];
      if (!x.cap || lvl[x.to] != lvl[u] + 1) continue;
      ll d = dfs(x.to, t, min(f - res, x.cap));
      x.cap -= d, g[x.to][x.rev].cap += d;
      if ((res += d) == f) return f;
    }
    return res;
  }
  ll max_flow(int s, int t) {
    ll f = 0;
    while (bfs(s, t))
      fill(all(it), 0), f += dfs(s, t, LLONG_MAX);
    return f;
  }
};
