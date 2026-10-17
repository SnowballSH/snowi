/**
 * Title: Floyd-Warshall
 * Description: All-pairs shortest paths; weights may be negative.
 *   On input `d[i][j]` is the weight of edge `i -> j`, `LLONG_MAX`
 *   if there is none. On output it is the distance: `LLONG_MAX` if
 *   `j` is unreachable from `i`, `LLONG_MIN` if a negative cycle on
 *   the way makes it arbitrarily small.
 * Precondition: `d` is $n times n$; keep the lightest of parallel
 *   edges. Distances without negative cycles stay below $2^62$ in
 *   absolute value.
 * Time: $O(n^3)$.
 * Usage:
 *   vector<vector<ll>> d(n, vector<ll>(n, LLONG_MAX));
 *   d[u][v] = min(d[u][v], w);
 *   floyd_warshall(d);
 * Tweak: Transitive closure: `bitset` rows and
 *   `if (r[i][k]) r[i] |= r[k]` in the same loop order (64 times
 *   faster).
 * Example:
 *   // Shortest cycle in a directed graph: run it with the diagonal
 *   // left at LLONG_MAX (drop the first loop), then
 *   ll best = LLONG_MAX;
 *   rep(i, 0, n) best = min(best, d[i][i]);
 *   // Path i -> j: store mid[i][j] = k whenever k improves d[i][j]
 *   // and recurse on (i, k), (k, j); -1 means a direct edge.
 *   // Add edge a -> b of weight w in O(n^2), if there is no
 *   // negative cycle before or after: d[b][a] == LLONG_MAX or
 *   // d[b][a] + w >= 0 (else rerun). Skip the sentinel.
 *   rep(i, 0, n) rep(j, 0, n)
 *     if (d[i][a] != LLONG_MAX && d[b][j] != LLONG_MAX)
 *       d[i][j] = min(d[i][j], d[i][a] + w + d[b][j]);
 * Notes: Use it for $n <= 500$ when many pairs are queried, for
 *   transitive closure, or when edges are added one at a time
 *   ($O(n^2)$ per edge while no negative cycle appears, Example).
 *   The loop over `k` must be the outermost. Processing vertices
 *   `k` in a chosen order (sorted by some value) answers "shortest
 *   path using only intermediate vertices with value $<= x$"
 *   offline. Without the clamp at $-2^62$, negative cycles make
 *   values grow exponentially and overflow.
 * Source: kactl FloydWarshall (CC0).
 * Status: stress-tested against walk-length dynamic programming
 */
#pragma once

void floyd_warshall(vector<vector<ll>>& d) {
  int n = sz(d);
  const ll inf = LLONG_MAX, lo = LLONG_MIN / 2;
  rep(i, 0, n) d[i][i] = min(d[i][i], 0LL);
  rep(k, 0, n)
    rep(i, 0, n)
      rep(j, 0, n)
        if (d[i][k] != inf && d[k][j] != inf)
          d[i][j] = min(d[i][j], max(d[i][k] + d[k][j], lo));
  rep(k, 0, n)
    rep(i, 0, n)
      rep(j, 0, n)
        if (d[k][k] < 0 && d[i][k] != inf && d[k][j] != inf)
          d[i][j] = LLONG_MIN;
}
