/**
 * Title: Segment tree
 * Description: Bottom-up point assign, range query over a monoid.
 *   The largest `r` with `pred(query(l, r))` true is
 *   `max_right(l, pred)`, and the smallest such `l` is
 *   `min_left(r, pred)`.
 * Precondition: 0-indexed, half-open ranges. `op` need not be
 *   commutative (left operand is the left part). `pred(e)` is true
 *   and `pred` is monotone: once false for a range, false for
 *   every range containing it.
 * Time: $O(n)$ build, $O(log n)$ per operation.
 * Usage:
 *   SegTree st(a); // or SegTree st(n): all leaves e
 *   st.set(i, 5); st.get(i); st.query(l, r); st.query_all();
 *   // position of the k-th one (k >= 1) in a 0/1 array, or n:
 *   int i = st.max_right(0, [&](ll s) { return s < k; });
 * Tweak: Max: `e = LLONG_MIN`; then `max_right(l, v < x)` is the
 *   first index `>= l` with `a[i] >= x`. A struct node needs `e` as
 *   a `constexpr` aggregate.
 * Example:
 *   // Point set range composite: f(x) = (b x + c) % MOD,
 *   // op(f, g) = "f, then g" (non-commutative).
 *   struct T { ll b, c; };
 *   static constexpr T e{1, 0};
 *   T op(T f, T g) {
 *     return {g.b * f.b % MOD, (g.b * f.c + g.c) % MOD};
 *   }
 *   // f = st.query(l, r) applies maps l..r-1 in order:
 *   // T f = st.query(l, r); ll y = (f.b * x + f.c) % MOD;
 * Notes: Think of it when a problem has point updates and range
 *   questions about any associative summary, or asks for "the first
 *   position after `l` where ...": `max_right` walks down the tree
 *   in $O(log n)$ instead of binary searching with queries. On
 *   sums, `s <= x` is monotone only when all values are `>= 0`. The
 *   query keeps a left and a right accumulator, so order is
 *   preserved: swapping them is invisible for max or sum but wrong
 *   for matrices or maps.
 *   Leaves past `n` hold `e`, so `e` must be a true identity (for
 *   (max,+) matrices: 0 on the diagonal, $-infinity$ elsewhere).
 *   With many test cases, build per test with that test's `n`.
 * Source: AtCoder Library segtree (CC0); kactl (CC0).
 * Status: stress-tested against brute force
 */
#pragma once

struct SegTree {
  typedef ll T;
  static constexpr T e = 0;
  T op(T a, T b) { return a + b; }
  int n, size = 1;
  vector<T> t;
  SegTree(int len) : SegTree(vector<T>(len, e)) {}
  SegTree(const vector<T>& a) : n(sz(a)) {
    while (size < n) size *= 2;
    t.assign(2 * size, e);
    rep(i, 0, n) t[size + i] = a[i];
    for (int i = size; --i;) pull(i);
  }
  void pull(int i) { t[i] = op(t[2 * i], t[2 * i + 1]); }
  void set(int i, T v) {
    for (t[i += size] = v; i /= 2;) pull(i);
  }
  T get(int i) { return t[i + size]; }
  T query_all() { return t[1]; }
  T query(int l, int r) {
    T a = e, b = e;
    for (l += size, r += size; l < r; l /= 2, r /= 2) {
      if (l & 1) a = op(a, t[l++]);
      if (r & 1) b = op(t[--r], b);
    }
    return op(a, b);
  }
  int max_right(int l, auto pred) {
    if (l == n) return n;
    T acc = e;
    l += size;
    do {
      while (l % 2 == 0) l /= 2;
      if (!pred(op(acc, t[l]))) {
        while (l < size) {
          l *= 2;
          if (pred(op(acc, t[l]))) acc = op(acc, t[l++]);
        }
        return l - size;
      }
      acc = op(acc, t[l++]);
    } while ((l & -l) != l);
    return n;
  }
  int min_left(int r, auto pred) {
    if (r == 0) return 0;
    T acc = e;
    r += size;
    do {
      r--;
      while (r > 1 && r % 2) r /= 2;
      if (!pred(op(t[r], acc))) {
        while (r < size) {
          r = 2 * r + 1;
          if (pred(op(t[r], acc))) acc = op(t[r--], acc);
        }
        return r + 1 - size;
      }
      acc = op(t[r], acc);
    } while ((r & -r) != r);
    return 0;
  }
};
