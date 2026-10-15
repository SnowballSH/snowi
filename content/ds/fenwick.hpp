/**
 * Title: Fenwick tree
 * Description: Point add and prefix sums. `lower_bound(s)` returns
 *   the smallest `i` whose prefix sum over `[0, i]` is at least `s`:
 *   `n` if none is, `-1` if `s <= 0`.
 * Precondition: 0-indexed. `lower_bound` needs all values `>= 0`.
 * Time: $O(log n)$ per operation.
 * Usage:
 *   Fenwick fw(n);
 *   fw.add(i, 5);
 *   ll s = fw.query(r) - fw.query(l); // sum of [l, r)
 * Tweak: Range add, point read: `add(l, x)`, `add(r, -x)`, value at
 *   `i` is `query(i + 1)`. Change `ll` for other abelian groups.
 * Example:
 *   // Range add + range sum with two trees b1, b2:
 *   // add x on [l, r):
 *   b1.add(l, x), b1.add(r, -x);
 *   b2.add(l, x * l), b2.add(r, -x * r);
 *   // sum of [0, p):
 *   b1.query(p) * p - b2.query(p);
 * Notes: Build in $O(n)$: set `t[i] = a[i]`, then for each `i` add
 *   `t[i]` to `t[i | (i + 1)]` when that index is `< n`.
 * Source: kactl FenwickTree (CC0).
 * Status: stress-tested against an array
 */
#pragma once

struct Fenwick {
  vector<ll> t;
  Fenwick(int n) : t(n) {}
  void add(int i, ll delta) {
    for (; i < sz(t); i |= i + 1) t[i] += delta;
  }
  ll query(int r) { // sum of [0, r)
    ll s = 0;
    for (; r > 0; r &= r - 1) s += t[r - 1];
    return s;
  }
  int lower_bound(ll s) {
    if (s <= 0) return -1;
    int pos = 0;
    for (int pw = 1 << 25; pw; pw >>= 1)
      if (pos + pw <= sz(t) && t[pos + pw - 1] < s)
        pos += pw, s -= t[pos - 1];
    return pos;
  }
};
