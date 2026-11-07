/**
 * Title: CDQ divide and conquer
 * Description: 3D dominance counting: `cnt[i]` is the number of
 *   points `j != i` with $a_j <= a_i$, $b_j <= b_i$ and
 *   $c_j <= c_i$; equal points count each other.
 * Precondition: Any `int` coordinates; `c` is compressed inside.
 * Time: $O(n log^2 n)$.
 * Usage:
 *   vi cnt = dominance(pts); // pts: vector<array<int, 3>>
 * Tweak: Strict in `b`: `<` in the sweep; in `c`:
 *   `query(c[idx[i]])`. Then equal points never count each other, so
 *   drop the final loop. Strict in `a` only: make it `b` or `c`.
 *   Weights: add and undo `wt` instead of 1, with `ll` counts;
 *   copies of a point get `cnt[last] + wt[last] - wt[i]`. DP over
 *   the order (`j` must be final before it feeds `i`): recurse
 *   left, do the cross step, then recurse right, sorting copies of
 *   both halves by `b` instead of merging.
 * Example:
 *   // CDQ on time: operation t inserts point (X[t], Y[t]) or,
 *   // if ask[t], counts inserted points <= (X[t], Y[t]). Keep
 *   // idx = 0..q-1 (time is a), merge by Y, Fenwick over X
 *   // (compressed). Only inserts of [l, m) add (and are undone),
 *   // only questions of [m, r) read; no copy fix.
 *   rep(i, m, r) {
 *     int t = idx[i];
 *     for (; k < m && Y[idx[k]] <= Y[t]; k++)
 *       if (!ask[idx[k]]) fw.add(X[idx[k]], 1);
 *     if (ask[t]) res[t] += fw.query(X[t] + 1);
 *   }
 * Notes: Think of it for "count pairs `j < i` with two more
 *   conditions", a 3D partial order, interleaved insertions and
 *   rectangle queries given offline (split a rectangle into four
 *   prefix corners with signs $plus.minus 1$), or
 *   `dp[i] = 1 + max dp[j]` over `j < i` with $x_j <= x_i$ and
 *   $y_j <= y_i$ (prefix-max Fenwick, reset touched cells by
 *   assignment). Undo the Fenwick inserts; clearing it per call
 *   costs $O(n)$ each time. Two conditions only: sort and sweep, no
 *   CDQ. Four: CDQ inside CDQ, $O(n log^3 n)$, or bitsets
 *   ($n^2 / 64$). Other cross steps work too: a merge of two convex
 *   hulls, a two-pointer count, a convolution. A static rectangle
 *   count needs no CDQ: sweep $x$ with a Fenwick tree over $y$.
 * Source: Chen Danqi (CDQ); UCup S4-S5 Training Packet 14.
 * Status: stress-tested against brute force
 */
#pragma once
#include "../ds/fenwick.hpp"

vi dominance(const vector<array<int, 3>>& p) {
  int n = sz(p);
  vi ord(n), c(n), cnt(n), cs;
  iota(all(ord), 0);
  sort(all(ord), [&](int i, int j) { return p[i] < p[j]; });
  for (auto& x : p) cs.push_back(x[2]);
  sort(all(cs));
  rep(i, 0, n)
    c[i] = int(lower_bound(all(cs), p[i][2]) - cs.begin());
  vi idx = ord;
  Fenwick fw(n);
  auto byb = [&](int i, int j) { return p[i][1] < p[j][1]; };
  auto rec = [&](auto&& self, int l, int r) -> void {
    if (r - l < 2) return;
    int m = (l + r) / 2, k = l;
    self(self, l, m), self(self, m, r); // now sorted by b
    rep(i, m, r) { // cross step: [l, m) counts into [m, r)
      for (; k < m && p[idx[k]][1] <= p[idx[i]][1]; k++)
        fw.add(c[idx[k]], 1);
      cnt[idx[i]] += int(fw.query(c[idx[i]] + 1));
    }
    rep(j, l, k) fw.add(c[idx[j]], -1);
    inplace_merge(idx.begin() + l, idx.begin() + m,
                  idx.begin() + r, byb);
  };
  rec(rec, 0, n);
  for (int i = n - 1; i-- > 0;) // copies count as the last one
    if (p[ord[i]] == p[ord[i + 1]])
      cnt[ord[i]] = cnt[ord[i + 1]];
  return cnt;
}
