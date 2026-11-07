/**
 * Title: Sparse table
 * Description: Static range minimum.
 * Precondition: 0-indexed, `l < r`. `op` is associative and
 *   idempotent, `op(x, x) = x`: min, max, gcd, bitwise and/or.
 * Time: $O(n log n)$ time and memory to build, $O(1)$ per query.
 * Usage:
 *   SparseTable st(a); // vector<ll> a
 *   ll m = st.query(l, r); // min of [l, r)
 * Example:
 *   // Leftmost position of the maximum of [l, r):
 *   vector<pair<ll, int>> b;
 *   rep(i, 0, sz(a)) b.push_back({-a[i], i});
 *   SparseTable st(b);
 *   int pos = st.query(l, r).second;
 * Notes: Use it when the array never changes and there are many
 *   queries: $10^6$ range minimum queries, LCA as the shallowest
 *   vertex between two Euler tour positions (pairs `{depth, v}`),
 *   or the LCP of two suffixes from the LCP array. A non-idempotent
 *   `op` (sum, product) double counts the overlap; use prefix sums
 *   instead. Memory is $n log n$ values: for $n = 10^6$ and `ll`
 *   that is 160 MB, so use `int` when values fit.
 * Source: kactl RMQ (CC0).
 * Status: stress-tested against brute force
 */
#pragma once

template<class T> struct SparseTable {
  static T op(T a, T b) { return min(a, b); }
  vector<vector<T>> t;
  SparseTable(const vector<T>& a) : t(1, a) {
    for (int k = 1, j = 1; 2 * k <= sz(a); k *= 2, j++) {
      t.emplace_back(sz(a) - 2 * k + 1);
      rep(i, 0, sz(t[j]))
        t[j][i] = op(t[j - 1][i], t[j - 1][i + k]);
    }
  }
  T query(int l, int r) {
    int j = __lg(r - l);
    return op(t[j][l], t[j][r - (1 << j)]);
  }
};
