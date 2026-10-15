/**
 * Title: Ternary search
 * Description: `ternary(lo, hi, f)` returns the smallest `x` in
 *   `[lo, hi]` that minimises an integer function `f`, by binary
 *   search for the first `x` with `f(x) <= f(x + 1)`.
 *   `golden(a, b, f)` returns a minimiser of a real function on
 *   `[a, b]` by golden-section search.
 * Precondition: `ternary`: `lo <= hi`; `f` strictly decreases up to
 *   its first minimum and never decreases after it. Every convex
 *   `f` qualifies (a flat bottom is fine); a flat stretch before
 *   the minimum breaks it. `golden`: `f` strictly decreases, then
 *   strictly increases, except on a flat bottom.
 * Time: `ternary`: $2 ceil(log_2 n)$ calls of `f` for $n$ values.
 *   `golden`: 102 calls; each shrinks the interval by 0.618.
 * Usage:
 *   ll x = ternary(0, n - 1, [&](ll i) { return cost(i); });
 *   double t = golden(0, 1e9, [&](double x) { return g(x); });
 * Tweak: In both, `f` may return any type with `<` and `<=` (`ll`,
 *   `double`, `pair`). Maximum: negate `f`. Use `long double` for
 *   more digits, fewer iterations when `f` is slow.
 * Example:
 *   // Points move along a line, x_i(t) = p[i] + v[i] t. The
 *   // spread max - min is convex in t (a max of lines minus a
 *   // min of lines): find the time when they are closest.
 *   auto spread = [&](double t) {
 *     double lo = 1e18, hi = -1e18;
 *     rep(i, 0, n) {
 *       double x = p[i] + v[i] * t;
 *       lo = min(lo, x), hi = max(hi, x);
 *     }
 *     return hi - lo;
 *   };
 *   double t = golden(0, 1e9, spread), best = spread(t);
 * Notes: Unimodality does not survive sums but convexity does:
 *   sums and maxima of convex functions are convex ($sum_i |x -
 *   a_i|$, the farthest of several points, the spread of moving
 *   points). If $F(x, y)$ is convex in both together,
 *   $min_y F(x, y)$ is convex in $x$: nest two searches. A product
 *   of monotone functions such as $i dot f(i)$ is often not
 *   unimodal: print `f` on small cases before trusting it. Near a
 *   smooth minimum `f` is flat, so `golden` finds the minimiser
 *   only to about $10^(-8)$ relative error, but the minimum value
 *   to full precision: print `f(x)`, not a formula in `x`.
 *   On integers prefer `ternary`; rounding `golden`'s answer can
 *   miss the best integer. Golden section costs one call per step
 *   (factor 0.618), plain ternary search two (factor 0.667 per
 *   step, 0.816 per call). Cues: "choose a time, position or
 *   parameter minimising the maximum (or the total) cost", "the
 *   answer first improves, then gets worse".
 * Source: kactl (CC0).
 * Status: stress-tested against brute force
 */
#pragma once

template<class F> ll ternary(ll lo, ll hi, F f) {
  while (lo < hi) {
    ll mid = lo + (hi - lo) / 2;
    if (f(mid) <= f(mid + 1)) hi = mid;
    else lo = mid + 1;
  }
  return lo;
}
template<class F> double golden(double a, double b, F f) {
  double r = (sqrt(5.0) - 1) / 2, x1 = b - r * (b - a),
         x2 = a + r * (b - a);
  auto f1 = f(x1), f2 = f(x2);
  rep(it, 0, 100) {
    if (f1 < f2) {
      b = x2, x2 = x1, f2 = f1;
      x1 = b - r * (b - a), f1 = f(x1);
    } else {
      a = x1, x1 = x2, f1 = f2;
      x2 = a + r * (b - a), f2 = f(x2);
    }
  }
  return (a + b) / 2;
}
