/**
 * Title: Wavelet matrix
 * Description: Static array of non-negative integers answering
 *   order statistics on any subarray: `kth(l, r, k)` is the
 *   `k`-th smallest value (0-indexed) of `a[l, r)`,
 *   `count_less(l, r, x)` counts its values `< x`,
 *   `count(l, r, lo, hi)` its values in `[lo, hi)`, and
 *   `rank(x, r)` the copies of `x` in `a[0, r)`. Level `h` stores
 *   bit `h` of every value as a bit vector with popcount prefix
 *   sums, then moves the values whose bit `h` is 0 to the front
 *   (stably): that is the order level `h - 1` sees.
 * Precondition: 0-indexed, half-open ranges, $0 <= a_i < 2^62$.
 *   `kth` needs `0 <= k < r - l`; any `ll` `x` is fine elsewhere.
 * Time: $O(n log V)$ build and $O(log V)$ per query, where $V$ is
 *   the largest value plus one. Memory: $(3/16) n log V$ bytes.
 * Usage:
 *   WaveletMatrix wm(a); // vector<ll> a
 *   ll v = wm.kth(l, r, k); // k-th smallest of a[l, r)
 *   int c = wm.count_less(l, r, x), d = wm.count(l, r, lo, hi);
 *   int copies = wm.rank(x, r); // of x in a[0, r)
 * Tweak: Negative or huge values: compress to indices in `vals`
 *   first (also makes $log V = log n$); map `kth` back through
 *   `vals`. Sum of the `k` smallest in `a[l, r)`: per level `h`,
 *   also keep prefix sums of the values in their order after level
 *   `h`'s partition (zeros first). In the `kth` walk, add the sum
 *   over `[l0, r0)` whenever it goes right; at the end add
 *   `k * res`, the copies of the result among the `k` smallest.
 * Example:
 *   WaveletMatrix wm({5, 1, 4, 1, 3, 9});
 *   wm.kth(1, 5, 2);          // 3: [1, 4, 1, 3] sorted is 1 1 3 4
 *   wm.count_less(1, 5, 4);   // 3
 *   wm.count(0, 6, 3, 6);     // 3: the values 5, 4, 3
 *   wm.rank(1, 4);            // 2
 *   // Largest value < x in a[l, r), if c > 0:
 *   int c = wm.count_less(l, r, x);
 *   ll prev = wm.kth(l, r, c - 1);
 *   // Smallest value >= x: kth(l, r, c) if c < r - l.
 *   // Median of a[l, r): kth(l, r, (r - l) / 2).
 * Notes: Think of it for static "k-th smallest in a range",
 *   "median of a subarray", "how many values in `[x, y)` lie in
 *   `a[l, r)`" and "closest value to `x` in a range", online and
 *   without updates. It needs far less memory than a persistent
 *   segment tree (bits instead of 16-byte nodes) and runs faster.
 *   Level `h` sends the values whose bit `h` is 0 to the front
 *   and the others behind all `zeros[h]` zeros, so a range maps
 *   to one range on either side by counting zeros, and a query
 *   walks one root-to-leaf path. The order at level `h` is not
 *   sorted by the higher bits as a number (bit `h + 1` is the
 *   main key), so data kept per level, like the sums in Tweak,
 *   must follow the partitions. On a tree, build it on the Euler
 *   tour for subtree queries. With point updates use offline
 *   methods (parallel binary search, sorted sweeps) or sqrt
 *   decomposition instead.
 * Source: Claude, Navarro, Ordonez, "The wavelet matrix" (2012).
 * Status: stress-tested against brute force
 */
#pragma once

struct WaveletMatrix {
  int n, lg = 0;
  vector<vector<uint64_t>> bits; // bits[h]: bit h, 64 per word
  vector<vi> pre; // pre[h][j]: ones in the first j words
  vi zeros;       // zeros[h]: zeros at level h
  WaveletMatrix(vector<ll> a) : n(sz(a)) {
    ll mx = 0;
    for (ll x : a) mx = max(mx, x);
    while (mx >> lg) lg++;
    bits.assign(lg, vector<uint64_t>(n / 64 + 1));
    pre.assign(lg, vi(n / 64 + 2));
    zeros.resize(lg);
    for (int h = lg; h--;) {
      auto& w = bits[h];
      vi& p = pre[h];
      rep(i, 0, n)
        if (a[i] >> h & 1) w[i / 64] |= 1ULL << i % 64;
      rep(j, 0, sz(w)) p[j + 1] = p[j] + popcount(w[j]);
      zeros[h] = n - p.back();
      stable_partition(all(a),
                       [&](ll x) { return ~x >> h & 1; });
    }
  }
  int rank0(int h, int i) { // zeros among the first i bits
    uint64_t w = bits[h][i / 64] & ((1ULL << i % 64) - 1);
    return i - pre[h][i / 64] - popcount(w);
  }
  ll kth(int l, int r, int k) {
    ll res = 0;
    for (int h = lg; h--;) {
      int l0 = rank0(h, l), r0 = rank0(h, r);
      if (k < r0 - l0) l = l0, r = r0;
      else {
        k -= r0 - l0, res |= 1LL << h;
        l += zeros[h] - l0, r += zeros[h] - r0;
      }
    }
    return res;
  }
  int count_less(int l, int r, ll x) {
    if (x >> lg) return x < 0 ? 0 : r - l;
    int res = 0;
    for (int h = lg; h--;) {
      int l0 = rank0(h, l), r0 = rank0(h, r);
      if (x >> h & 1) {
        res += r0 - l0;
        l += zeros[h] - l0, r += zeros[h] - r0;
      } else l = l0, r = r0;
    }
    return res;
  }
  int count(int l, int r, ll lo, ll hi) { // 0 if lo >= hi
    return max(0, count_less(l, r, hi) - count_less(l, r, lo));
  }
  int rank(ll x, int r) { // x >> lg: absent; no x + 1 overflow
    return x >> lg ? 0 : count(0, r, x, x + 1);
  }
};
