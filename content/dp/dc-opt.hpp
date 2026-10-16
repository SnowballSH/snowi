/**
 * Title: Divide and conquer optimisation
 * Description: One layer of a DP,
 *   $"new"[i] = min_(j < i) "old"[j] + "cost"(j, i)$ for every
 *   `i`, when the smallest optimal `j` is non-decreasing in `i`.
 * Precondition: `cost` satisfies the quadrangle inequality (DP
 *   techniques), or the smallest optimal `j` is known to be
 *   non-decreasing in `i` (a monotone choice among tied `j` is not
 *   enough). Unreachable states are `LLONG_MAX / 4`, also used
 *   for `new[0]` and as a cap; `old[j] + cost(j, i)` fits in `ll`.
 * Time: $O(n log n)$ calls to `cost` per layer.
 * Usage:
 *   // dp[i]: best split of the first i elements; cost(j, i)
 *   // is the cost of one group [j, i).
 *   vector<ll> dp(n + 1, LLONG_MAX / 4);
 *   dp[0] = 0;
 *   rep(k, 0, K) dp = dc_layer(dp, cost); // dp[n]: K groups
 * Tweak: Only some `i` needed, or `j <= i` allowed: change the
 *   root call and the scan bound `m`. A cost without an $O(1)$
 *   formula (equal pairs in `[j, i)`): move a Mo-style window to
 *   `[j, i)` inside `cost`, $O(n log n)$ moves per layer. All `K`
 *   layers in $O(n^2)$: Knuth's bound
 *   $"opt"_(k - 1)[i] <= "opt"_k [i] <= "opt"_k [i + 1]$, `i`
 *   filled downwards.
 * Example:
 *   // Split a into K groups, minimising the sum of the squared
 *   // group sums; pre[i] = a[0] + ... + a[i - 1], a[i] >= 0.
 *   auto cost = [&](int j, int i) {
 *     ll s = pre[i] - pre[j];
 *     return s * s;
 *   };
 *   vector<ll> dp(n + 1, LLONG_MAX / 4);
 *   dp[0] = 0;
 *   rep(k, 0, K) dp = dc_layer(dp, cost);
 *   ll ans = dp[n];
 * Notes: The quadrangle inequality says crossing intervals cost
 *   at most nested ones. It holds for $f(S_i - S_j)$ with $f$
 *   convex and $S$ non-decreasing (squared sums), for "pairs
 *   inside the group" counts such as inversions or equal pairs,
 *   and for sums of such costs. When unsure, brute-force `opt` on
 *   small cases and check that it never decreases. Window costs
 *   (Tweak) also fit distinct values for a maximum. If the cost
 *   has the quadrangle inequality, the answer for exactly `K`
 *   groups is convex in `K`, so the Aliens trick (`dp/wqs`) with a
 *   fast inner DP avoids the factor `K`. This does not apply to
 *   $"dp"[i] = min_(j < i) "dp"[j] + "cost"(j, i)$ on a single
 *   array, where `dp[j]` is still unknown: use CDQ divide and
 *   conquer on `i`, or a deque of candidates with binary search.
 * Source: Classic; cp-algorithms "Divide and Conquer DP".
 * Status: stress-tested against an $O(n^2)$ DP
 */
#pragma once

template<class F>
vector<ll> dc_layer(const vector<ll>& old, F cost) {
  int n = sz(old);
  vector<ll> res(n, LLONG_MAX / 4);
  // Solve i in [l, r) knowing the optimal j lies in [lo, hi].
  auto rec = [&](auto& self, int l, int r, int lo, int hi) {
    if (l >= r) return;
    int m = (l + r) / 2, best = lo;
    ll val = LLONG_MAX;
    rep(j, lo, min(m, hi + 1)) {
      ll v = old[j] + cost(j, m);
      if (v < val) val = v, best = j;
    }
    res[m] = min(res[m], val);
    self(self, l, m, lo, best);
    self(self, m + 1, r, best, hi);
  };
  rec(rec, 1, n, 0, n - 1);
  return res;
}
