/**
 * Title: Global minimum cut
 * Description: Stoer-Wagner: the cheapest way to split an
 *   undirected weighted graph into two non-empty parts. Returns
 *   the cut weight and the vertices of one part.
 * Precondition: Adjacency matrix with `w[u][v] = w[v][u] >= 0`
 *   (sum parallel edges; the diagonal is ignored), $n >= 2$.
 *   The edge weights sum to less than $2^63$.
 * Time: $O(n^3)$; about 0.5 s for $n = 1000$.
 * Usage:
 *   vector<vector<ll>> w(n, vector<ll>(n));
 *   w[u][v] += c, w[v][u] += c;
 *   auto [cut, side] = global_min_cut(w);
 * Tweak: Real weights: `double` for `ll`, `INFINITY` for
 *   `LLONG_MAX`. Sparse graphs: adjacency lists and a max-heap of
 *   `key` (lazy deletion) give $O(n m log n)$.
 * Example:
 *   // Fewest roads to close so that the towns fall apart, and
 *   // the towns on one side (edge connectivity).
 *   vector<vector<ll>> w(n, vector<ll>(n));
 *   for (auto [a, b] : roads) w[a][b]++, w[b][a]++;
 *   auto [cut, side] = global_min_cut(w); // 0: disconnected
 * Notes: Each phase adds vertices one by one, always the one
 *   most tightly connected to those already added; the last
 *   one, $t$, is separated from the one before, $s$, by exactly
 *   its `key`, which is a minimum $s$-$t$ cut. Then $s$ and $t$
 *   merge. Cues: "disconnect the graph" with no $s$ and $t$
 *   given, "fewest edges whose removal splits the network",
 *   "split into two non-empty groups, minimizing the weight
 *   between them". The same value is the minimum over $t$ of
 *   the $s$-$t$ max flow for any fixed $s$ ($n - 1$ flows), and
 *   the lightest edge of a Gomory-Hu tree. Directed graphs have
 *   no Stoer-Wagner: take the minimum of the flows from a fixed
 *   $s$ to every $t$ and back, $2(n - 1)$ flows.
 * Source: Stoer and Wagner (1997); kactl GlobalMinCut (CC0).
 * Status: stress-tested against enumeration of all cuts
 */
#pragma once

pair<ll, vi> global_min_cut(vector<vector<ll>> w) {
  int n = sz(w);
  pair<ll, vi> best{LLONG_MAX, {}};
  vector<vi> group(n); // original vertices merged into v
  rep(v, 0, n) group[v] = {v}, w[v][v] = 0;
  vi alive(n);
  iota(all(alive), 0);
  while (sz(alive) > 1) {
    vector<ll> key(n);
    vi rest = alive;
    int s = -1, t = -1;
    ll cut = 0;
    while (sz(rest)) {
      int i = 0;
      rep(j, 1, sz(rest))
        if (key[rest[j]] > key[rest[i]]) i = j;
      s = t, t = rest[i], cut = key[t];
      rest[i] = rest.back(), rest.pop_back();
      for (int v : rest) key[v] += w[t][v];
    }
    if (cut < best.first) best = {cut, group[t]};
    group[s].insert(group[s].end(), all(group[t]));
    for (int v : alive) w[s][v] += w[t][v], w[v][s] = w[s][v];
    alive.erase(find(all(alive), t));
  }
  return best;
}
