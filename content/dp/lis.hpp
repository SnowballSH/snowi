/**
 * Title: Longest increasing subsequence
 * Description: Indices of one longest strictly increasing
 *   subsequence of `a`, in increasing order. While it runs,
 *   `k + 1` is the length of the longest one ending at `i`.
 * Precondition: Values only need `<`. `a` may be empty.
 * Time: $O(n log n)$.
 * Usage:
 *   vi idx = lis(a); // a[idx[0]] < a[idx[1]] < ...
 *   int len = sz(idx);
 * Tweak: Longest non-decreasing: `upper_bound` instead of
 *   `lower_bound`. Decreasing: negate the values, or search with
 *   `greater<T>()` as a fourth argument.
 * Example:
 *   // Most envelopes nested strictly in both width and height:
 *   // sort by width up, height down for equal widths, so two
 *   // equal widths never chain; then LIS on the heights.
 *   sort(all(env), [](pii p, pii q) {
 *     return p.first != q.first ? p.first < q.first
 *                               : p.second > q.second;
 *   });
 *   vi h;
 *   for (auto [w, ht] : env) h.push_back(ht);
 *   int best = sz(lis(h));
 * Notes: `tail[k]` is the smallest last value of an increasing
 *   subsequence of length `k + 1` seen so far; it stays sorted.
 *   Cues: fewest deletions to make an array strictly increasing
 *   is `n - LIS` (non-decreasing: the Tweak); fewest changes to
 *   make it strictly increasing with any integers is `n` minus
 *   the longest non-decreasing subsequence of `a[i] - i` (if
 *   values must stay `>= 1`, first drop every `i` with
 *   `a[i] <= i`). Dilworth: the fewest non-increasing
 *   subsequences covering `a` equals the LIS
 *   length. Pairs ordered in two coordinates: sort by one, LIS on
 *   the other, with the tie order chosen as in the Example. `i`
 *   lies in some LIS iff $"end"[i] + "start"[i] - 1 = L$ (`start`
 *   from the reversed, negated array); in every LIS iff also no
 *   other such `j` has the same `end[j]`. Counting LIS: a Fenwick
 *   tree over values storing (length, count) pairs.
 * Source: kactl LIS (CC0).
 * Status: stress-tested against an $O(n^2)$ DP
 */
#pragma once

template<class T> vi lis(const vector<T>& a) {
  vector<T> tail;
  vi at, prev(sz(a)); // at[k]: index of tail[k]
  rep(i, 0, sz(a)) {
    int k = lower_bound(all(tail), a[i]) - begin(tail);
    if (k == sz(tail)) tail.push_back(a[i]), at.push_back(i);
    tail[k] = a[i], at[k] = i;
    prev[i] = k ? at[k - 1] : -1;
  }
  vi res;
  for (int i = sz(at) ? at.back() : -1; i >= 0; i = prev[i])
    res.push_back(i);
  reverse(all(res));
  return res;
}
