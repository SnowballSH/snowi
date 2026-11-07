/**
 * Title: Knuth optimisation
 * Description: The interval DP
 *   $"dp"[l][r] = w(l, r) + min_(l < k < r) "dp"[l][k] + "dp"[k][r]$
 *   for all $0 <= l < r <= n$, with $"dp"[i][i + 1] = 0$, in
 *   $O(n^2)$ instead of $O(n^3)$: the smallest optimal split
 *   satisfies $"opt"[l][r - 1] <= "opt"[l][r] <= "opt"[l + 1][r]$.
 * Precondition: `w` satisfies the quadrangle inequality (DP
 *   techniques) and is monotone, $w(b, c) <= w(a, d)$ for
 *   $a <= b <= c <= d$. Sums of values over the
 *   interval, with values $>= 0$, satisfy both.
 * Time: $O(n^2)$ time and memory (12 bytes per pair: 300 MB at
 *   $n = 5000$).
 * Usage:
 *   // merge adjacent piles; a merge costs the sum of both
 *   auto dp = knuth(n, [&](int l, int r) {
 *     return pre[r] - pre[l];
 *   });
 *   ll ans = dp[0][n];
 * Tweak: Leaf values: set `dp[i][i + 1]` before the loops; they
 *   do not change `opt`.
 * Example:
 *   // Optimal binary search tree on keys 0..n-1 searched f[i]
 *   // times; cost = sum of f[i] * depth (root depth 1). Root k
 *   // splits keys [l, k) and [k + 1, r); to fit the form, let
 *   // dp[l][r] cover keys l..r-2 (dp[l][l + 1] is empty).
 *   auto dp = knuth(n + 1, [&](int l, int r) {
 *     return pre[r - 1] - pre[l]; // pre: prefix sums of f
 *   });
 *   ll ans = dp[0][n + 1];
 * Notes: Signals: merging neighbours where a merge costs the
 *   total size (slimes, file merging, optimal alphabetic code),
 *   cutting a stick or string where a cut costs the current
 *   length (cut positions as boundaries), optimal BSTs. When
 *   unsure, test both conditions on all quadruples of small
 *   random inputs. If any two piles may merge (not only
 *   neighbours), it is Huffman instead: a greedy with a heap.
 *   Merging `K > 2` piles at once breaks the form; use an
 *   $O(n^3 / K)$ DP over `(l, r)` alone: `[l, r)` ends as
 *   $(r - l - 1) mod (K - 1) + 1$ piles, and a split's left part
 *   as one. The same `opt` bounds speed up layer DPs
 *   (`dp/dc-opt`).
 * Source: D. Knuth (1971), F. Yao (1980).
 * Status: stress-tested against the $O(n^3)$ DP
 */
#pragma once

template<class W> vector<vector<ll>> knuth(int n, W w) {
  vector<vector<ll>> dp(n + 1, vector<ll>(n + 1));
  vector<vi> opt(n + 1, vi(n + 1));
  rep(i, 0, n) opt[i][i + 1] = i;
  for (int l = n - 2; l >= 0; l--)
    rep(r, l + 2, n + 1) {
      ll best = LLONG_MAX;
      rep(k, max(l + 1, opt[l][r - 1]), opt[l + 1][r] + 1)
        if (ll v = dp[l][k] + dp[k][r]; v < best)
          best = v, opt[l][r] = k;
      dp[l][r] = best + w(l, r);
    }
  return dp;
}
