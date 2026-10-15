/**
 * Title: Gomory-Hu tree
 * Description: Builds a tree on the same vertices such that the
 *   minimum cut between any `u` and `v` equals the lightest edge
 *   on their tree path; removing that edge splits the vertices
 *   into the two sides of such a cut. Returns `n - 1` edges
 *   `{u, v, w}`.
 * Precondition: Undirected edges `{u, v, c}`, 0-indexed,
 *   `c >= 0`; parallel edges and self-loops are fine. The same
 *   overflow bounds as `flow/dinic`.
 * Time: $n - 1$ max flows; a random graph with $n = 1000$,
 *   $m = 10^4$ takes about 0.2 s.
 * Usage:
 *   vector<array<ll, 3>> edges; // {u, v, capacity}
 *   auto tree = gomory_hu(n, edges); // {u, v, min cut}
 * Tweak: The tree is rooted at 0: `tree[v - 1]` holds `v`, its
 *   parent and the cut, ready for binary lifting of path minima.
 * Uses: flow/dinic
 * Example:
 *   // Sum of the minimum cut over all pairs: join tree edges
 *   // heaviest first; joining sets of sizes a and b gives a * b
 *   // pairs whose minimum cut is that edge (ds/dsu).
 *   sort(all(tree), [](auto& x, auto& y) { return x[2] > y[2]; });
 *   DSU dsu(n);
 *   ll total = 0;
 *   for (auto [u, v, c] : tree)
 *     total += c * dsu.size(u) * dsu.size(v), dsu.join(u, v);
 * Notes: Gusfield's method: no graph contraction, each step one
 *   max flow on the original graph between `s` and its current
 *   tree parent. Answers all $binom(n, 2)$ pairwise min cuts,
 *   that is max flows (with unit capacities: numbers of
 *   edge-disjoint paths), with $n - 1$ flows; the lightest tree
 *   edge is the global minimum cut. Cues: "for every pair", "$q$
 *   queries of the min cut between $u$ and $v$", "sum or count
 *   of pairs by connectivity", $n <= 3000$ (about 1 s for
 *   $m = 4500$). Undirected only: directed graphs have no such
 *   tree. Queries: min on a tree path (binary lifting) or the
 *   DSU order above, offline.
 * Source: Gomory and Hu (1961); Gusfield (1990), "Very simple
 *   methods for all pairs network flow analysis".
 * Status: stress-tested against cut enumeration and pairwise flows
 */
#pragma once
#include "dinic.hpp"

auto gomory_hu(int n, const vector<array<ll, 3>>& edges) {
  Dinic d(n);
  for (auto [u, v, c] : edges) d.add_edge(u, v, c, c);
  vi par(n);
  vector<ll> cut(n);
  rep(s, 1, n) {
    for (auto& adj : d.g)
      for (auto& x : adj) x.cap = x.orig;
    int t = par[s];
    cut[s] = d.max_flow(s, t);
    rep(i, 0, n)
      if (i != s && par[i] == t && d.source_side(i)) par[i] = s;
    if (d.source_side(par[t]))
      par[s] = par[t], par[t] = s, swap(cut[s], cut[t]);
  }
  vector<array<ll, 3>> tree;
  rep(v, 1, n) tree.push_back({v, par[v], cut[v]});
  return tree;
}
