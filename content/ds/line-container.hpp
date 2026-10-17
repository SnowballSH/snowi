/**
 * Title: Line container
 * Description: Upper envelope of lines $k x + m$ added in any order;
 *   `query(x)` returns the maximum value at `x`.
 * Precondition: `query` needs at least one line. $|k|, |m| < 2^62$
 *   and every queried value $k x + m$ fits in `ll`.
 * Time: $O(log n)$ amortized per operation.
 * Usage:
 *   LineContainer lc;
 *   lc.add(k, m); ll y = lc.query(x); // max of k x + m
 * Tweak: Real values: `ll` becomes `double`, `inf = HUGE_VAL` and
 *   `div(a, b) = a / b`.
 * Example:
 *   // dp[j] = max over i < j of dp[i] + a[i] * b[j], with the
 *   // b[j] in any order (a monotone stack needs sorted b).
 *   LineContainer lc;
 *   lc.add(a[0], dp[0]);
 *   rep(j, 1, n) {
 *     dp[j] = lc.query(b[j]);
 *     lc.add(a[j], dp[j]);
 *   }
 * Notes: Does the job of the Li Chao tree (`ds/li-chao`) with
 *   online queries and unbounded $x$, but takes full lines only
 *   and cannot undo. PyPy has no balanced tree: use the Li Chao
 *   tree there (read the queries first and compress them). Lines
 *   are kept sorted by slope; `isect` finds `p` by floor division,
 *   so equal slopes and integer breakpoints are handled exactly. A
 *   new line that is never the maximum is erased at once, and the
 *   neighbours it hides are erased while it is added, so each line
 *   is erased at most once. Signals: a DP term $a_i b_j$ queried
 *   in any order; maximizing $sum a - lambda sum b$ for many
 *   $lambda$ (Dinkelbach, parametric search), where each item is a
 *   line in $lambda$; the farthest point of a set in direction
 *   $(x, 1)$, a dot product.
 * Source: kactl LineContainer (CC0), by Simon Lindholm.
 * Status: stress-tested against brute force
 */
#pragma once

struct HullLine {
  mutable ll k, m, p;
  bool operator<(const HullLine& o) const { return k < o.k; }
  bool operator<(ll x) const { return p < x; }
};

struct LineContainer : multiset<HullLine, less<>> {
  static constexpr ll inf = LLONG_MAX;
  ll div(ll a, ll b) { return a / b - ((a ^ b) < 0 && a % b); }
  bool isect(iterator x, iterator y) {
    if (y == end()) return x->p = inf, 0;
    if (x->k == y->k) x->p = x->m > y->m ? inf : -inf;
    else x->p = div(y->m - x->m, x->k - y->k);
    return x->p >= y->p;
  }
  void add(ll k, ll m) {
    auto z = insert({k, m, 0}), y = z++, x = y;
    while (isect(y, z)) z = erase(z);
    if (x != begin() && isect(--x, y)) isect(x, y = erase(y));
    while ((y = x) != begin() && (--x)->p >= y->p)
      isect(x, erase(y));
  }
  ll query(ll x) {
    auto l = *lower_bound(x);
    return l.k * x + l.m;
  }
};
