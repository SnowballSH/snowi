/**
 * Title: Slope trick
 * Description: A convex piecewise-linear function $f$ with
 *   integer slopes, stored as its minimum `min_f` and its
 *   breakpoints: `L` (max-heap) left of the minimum, where the
 *   slope drops by 1 per point going left, and `R` (min-heap)
 *   right of it. The minimum is reached on `[top_l(), top_r()]`
 *   (`INF` means unbounded). Shifts are lazy.
 * Precondition: Breakpoints, shifts and values stay below
 *   $10^18$ in absolute value. `slide` needs `a <= b`.
 * Time: $O(log n)$ per hinge, $O(1)$ for the rest (clearing is
 *   amortised), $O(n log n)$ for `eval`.
 * Usage:
 *   SlopeTrick st; // f = 0
 *   st.add_abs(a); // f += |x - a|; also add_x_minus_a,
 *   // add_a_minus_x (one-sided hinges), add_const
 *   st.prefix_min(); // f(x) = min of f(y) over y <= x
 *   st.slide(a, b); // f(x) = min of f(y), x - b <= y <= x - a
 *   ll lo = st.top_l(), hi = st.top_r(), v = st.min_f;
 * Tweak: Weight $c dot |x - a|$: call `c` times, or store
 *   `(point, count)` pairs and split a pair when only part of it
 *   moves. Start pinned at $x_0$ ($f = 0$ there, $infinity$
 *   elsewhere): call `add_abs(x0)` `K` times, `K` above any slope
 *   added later. Domain $x >= 0$: add `K` times
 *   `add_a_minus_x(0)`, or show the optimum never leaves it.
 * Example:
 *   // Fewest +-1 steps making a non-decreasing (CF 13C):
 *   // f_i(x) = |x - a_i| + min over y <= x of f_(i-1)(y).
 *   SlopeTrick st;
 *   for (ll x : a) st.prefix_min(), st.add_abs(x);
 *   ll steps = st.min_f;
 *   // An optimal b: store opt[i] = st.top_l() after step i,
 *   // then b[n-1] = opt[n-1], b[i] = min(b[i + 1], opt[i]).
 *   // Strictly increasing: run on a[i] - i. Neighbours within
 *   // d, |b[i] - b[i-1]| <= d: slide(-d, d) in place of
 *   // prefix_min(). Rising by 0 to d: slide(0, d) instead.
 * Notes: Think of it when a DP has one numeric state `x` (a
 *   value, position or time) and costs like $|x - a|$ or
 *   $max(0, x - a)$, with transitions "min over $y <= x$", "min
 *   over a window around $x$", "shift", or "add two functions".
 *   Signals: make an array monotone with fewest $plus.minus 1$
 *   steps, move at speed at most 1 and pay distances at given
 *   times (`slide(-t, t)`), equalise root-to-leaf path lengths
 *   on a tree (merge children's heaps small-to-large). Adding a
 *   right hinge pays `top_l - a` if `a` is left of the minimum,
 *   then moves the new largest left point to `R`; the left hinge
 *   mirrors it. `min` of two different convex functions or a
 *   concave term break convexity. Track `min_f`; `eval` is for
 *   the final answer only. A concave function of a count works
 *   the same way with roles swapped: store its marginal gains in
 *   a heap (buy low, sell high: for each price `p`, if the
 *   cheapest kept price `q < p`, gain `p - q`, pop `q`, push `p`;
 *   always push `p`).
 * Source: Packet 13 (UCup); "Slope trick" by maspy.
 * Status: stress-tested against brute force on a grid
 */
#pragma once

struct SlopeTrick {
  static constexpr ll INF = LLONG_MAX / 4;
  ll min_f = 0, add_l = 0, add_r = 0;
  priority_queue<ll> L;
  priority_queue<ll, vector<ll>, greater<ll>> R;
  ll top_l() { return L.empty() ? -INF : L.top() + add_l; }
  ll top_r() { return R.empty() ? INF : R.top() + add_r; }
  void push_l(ll a) { L.push(a - add_l); }
  void push_r(ll a) { R.push(a - add_r); }
  void add_const(ll c) { min_f += c; }
  void add_x_minus_a(ll a) { // f += max(0, x - a)
    min_f += max(0LL, top_l() - a);
    push_l(a), push_r(top_l()), L.pop();
  }
  void add_a_minus_x(ll a) { // f += max(0, a - x)
    min_f += max(0LL, a - top_r());
    push_r(a), push_l(top_r()), R.pop();
  }
  void add_abs(ll a) { add_x_minus_a(a), add_a_minus_x(a); }
  void prefix_min() { R = {}; }
  void suffix_min() { L = {}; }
  void slide(ll a, ll b) { add_l += a, add_r += b; }
  ll eval(ll x) {
    ll res = min_f;
    for (auto l = L; sz(l); l.pop())
      res += max(0LL, l.top() + add_l - x);
    for (auto r = R; sz(r); r.pop())
      res += max(0LL, x - r.top() - add_r);
    return res;
  }
};
