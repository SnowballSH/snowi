/**
 * Title: Hungarian algorithm
 * Description: Minimum-cost assignment of the `n` rows of `a` to
 *   distinct columns. Returns the cost and `match`, where
 *   `match[i]` is the column of row `i`.
 * Precondition: `n <= m` for an `n x m` matrix (transpose if not).
 *   Costs may be negative, with $(n + 4) max |a_(i j)| < 2^63$.
 * Time: $O(n^2 m)$.
 * Usage:
 *   vector<vector<ll>> a(n, vector<ll>(m));
 *   auto [cost, match] = hungarian(a);
 * Tweak: Rows that may stay unassigned: add $n$ dummy columns of
 *   cost 0. Real costs: replace every `ll` with `double` (the
 *   `LLONG_MAX` sentinels still work).
 * Example:
 *   // Maximum total score, each worker on a different job:
 *   for (auto& row : a) for (ll& x : row) x = -x;
 *   auto [c, match] = hungarian(a);
 *   ll best = -c;
 * Notes: This is successive shortest paths on the complete
 *   bipartite graph with dense $O(m)$ Dijkstra steps: rows are
 *   added one at a time, and potentials `u` (rows) and `v`
 *   (columns) keep $a_(i j) - u_i - v_j >= 0$, with equality on the
 *   matched pairs. Think of it for "assign workers to jobs", "match
 *   $n$ points to $n$ targets with minimum total distance", or a
 *   permutation minimizing $sum_i a_(i, p(i))$. For $n = 500$ it
 *   takes well under a second. Sparse or capacitated variants fit
 *   min-cost flow better.
 * Source: e-maxx (cp-algorithms) Hungarian with potentials; kactl
 *   (CC0).
 * Status: stress-tested against permutations and planted optima
 */
#pragma once

pair<ll, vi> hungarian(const vector<vector<ll>>& a) {
  int n = sz(a), m = n ? sz(a[0]) : 0;
  vector<ll> u(n), v(m + 1);   // column m is a virtual start
  vi p(m + 1, -1), way(m + 1); // p[j]: row matched to column j
  rep(i, 0, n) {
    vector<ll> dist(m, LLONG_MAX);
    vector<bool> done(m + 1);
    int j0 = m;
    p[m] = i;
    do { // Dijkstra step: grow the tree by the closest column
      done[j0] = true;
      int i0 = p[j0], j1 = 0;
      ll delta = LLONG_MAX;
      rep(j, 0, m)
        if (!done[j]) {
          ll cur = a[i0][j] - u[i0] - v[j];
          if (cur < dist[j]) dist[j] = cur, way[j] = j0;
          if (dist[j] < delta) delta = dist[j], j1 = j;
        }
      rep(j, 0, m + 1) {
        if (done[j]) u[p[j]] += delta, v[j] -= delta;
        else dist[j] -= delta;
      }
      j0 = j1;
    } while (p[j0] >= 0);
    for (int j1; j0 != m; j0 = j1) p[j0] = p[j1 = way[j0]];
  }
  vi match(n);
  rep(j, 0, m)
    if (p[j] >= 0) match[p[j]] = j;
  return {-v[m], match};
}
