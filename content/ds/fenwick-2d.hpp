/**
 * Title: 2D Fenwick tree
 * Description: Point add and rectangle sums on a dense `n` x `m`
 *   grid. `query(x, y)` is the sum over `[0, x) x [0, y)`;
 *   `rect(x1, y1, x2, y2)` the sum over `[x1, x2) x [y1, y2)`.
 * Precondition: 0-indexed, `0 <= x <= n`, `0 <= y <= m`; `add` at
 *   `x = n` or `y = m` does nothing. Memory is $n m$ values, so
 *   $n m <= 10^7$ or so.
 * Time: $O(log n log m)$ per operation.
 * Usage:
 *   Fenwick2D fw(n, m);
 *   fw.add(x, y, 5);
 *   ll s = fw.rect(x1, y1, x2, y2);
 * Tweak: Rectangle add with point query: add `v` at `(x1, y1)` and
 *   `(x2, y2)`, `-v` at `(x1, y2)` and `(x2, y1)`; the value at
 *   `(x, y)` is `query(x + 1, y + 1)`. Other abelian groups: change
 *   `ll` and the `+=`, `-`, `+` (xor: all become `^`).
 * Example:
 *   // Count points in [x1, x2) x [y1, y2) after insertions:
 *   Fenwick2D cnt(n, m);
 *   for (auto [x, y] : pts) cnt.add(x, y, 1);
 *   ll inside = cnt.rect(x1, y1, x2, y2);
 * Notes: For sparse points with big coordinates, go offline: sort
 *   the distinct `x` of all updates, and for each Fenwick node `i`
 *   collect (sorted, unique) the `y` of updates whose `x` passes
 *   through `i`. Each node then holds a 1D Fenwick tree over its
 *   own `y` list, searched by `lower_bound`: $O(log^2 n)$ per
 *   operation and $O(k log k)$ memory for `k` updates. If there are
 *   no updates between queries, a sweep over `x` with a 1D Fenwick
 *   tree over `y` is simpler and faster. Cues: "count points in a
 *   rectangle", "sum of a submatrix with point changes".
 * Source: kactl FenwickTree2d (CC0) for the offline variant.
 * Status: stress-tested against brute force
 */
#pragma once

struct Fenwick2D {
  int n, m;
  vector<vector<ll>> t;
  Fenwick2D(int rows, int cols)
      : n(rows), m(cols), t(rows, vector<ll>(cols)) {}
  void add(int x, int y, ll delta) {
    for (int i = x; i < n; i |= i + 1)
      for (int j = y; j < m; j |= j + 1) t[i][j] += delta;
  }
  ll query(int x, int y) { // sum of [0, x) x [0, y)
    ll s = 0;
    for (int i = x; i > 0; i &= i - 1)
      for (int j = y; j > 0; j &= j - 1) s += t[i - 1][j - 1];
    return s;
  }
  ll rect(int x1, int y1, int x2, int y2) {
    return query(x2, y2) - query(x1, y2) - query(x2, y1) +
           query(x1, y1);
  }
};
