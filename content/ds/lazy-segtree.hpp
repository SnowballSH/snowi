/**
 * Title: Lazy segment tree
 * Description: Range apply and range fold. Nodes `S` form a monoid
 *   under `op` with identity `e`; tags `F` act on nodes through
 *   `mapping(f, s)` and compose by `composition(f, g)` = "`g`, then
 *   `f`", with identity `id`. `max_right` and `min_left` work as in
 *   the segment tree. Default: range add, range min.
 * Precondition: 0-indexed, half-open ranges. `mapping(f, op(a, b))
 *   = op(mapping(f, a), mapping(f, b))` and `mapping(id, s) = s`.
 *   `pred(e)` is true and `pred` is monotone.
 * Time: $O(n)$ build, $O(log n)$ per operation.
 * Usage:
 *   LazySegTree st(a); // vector<ll> a(n)
 *   st.apply(l, r, 5); st.apply(i, -1); st.query(l, r);
 *   st.set(i, 7); st.get(i); st.query_all();
 *   int r = st.max_right(l, [&](ll v) { return v > x; });
 * Tweak: Change only the block above `n`. Max instead of min:
 *   `e = LLONG_MIN`, `max`. Range sum needs the length:
 *   `S = {sum, len}` with leaves `{a_i, 1}` and `mapping` returning
 *   `{s.sum + f * s.len, s.len}`. Range assign: `F` with a sentinel
 *   `id` for "none", `composition` returns `f == id ? g : f` and
 *   `mapping` returns `f == id ? s : f` (`push` also passes `id`).
 * Example:
 *   // Range affine a_i = b a_i + c, range sum, mod MOD:
 *   struct S { ll sum, len; };
 *   struct F { ll b, c; };
 *   static constexpr S e{0, 0};
 *   static constexpr F id{1, 0};
 *   S op(S x, S y) {
 *     return {(x.sum + y.sum) % MOD, x.len + y.len};
 *   }
 *   S mapping(F f, S s) {
 *     return {(f.b * s.sum + f.c * s.len) % MOD, s.len};
 *   }
 *   F composition(F f, F g) { // g, then f
 *     return {f.b * g.b % MOD, (f.b * g.c + f.c) % MOD};
 *   }
 *   // Build from {a_i, 1}; st.apply(l, r, {b, c});
 *   // st.query(l, r).sum
 * Notes: Check the three laws on paper before coding: tags compose,
 *   act on a whole node in $O(1)$, and distribute over `op`. Sums
 *   need `len`; products of tags must stay in a small family. The
 *   new tag goes on the outside: `composition(new, old)`. If a tag
 *   cannot update a node's summary without looking inside (chmin
 *   with sum, `%=` with sum), use segment tree beats. `max_right`
 *   finds e.g. the first index after `l` with value `<= x` under
 *   range adds: `pred = v > x` on the min. Leaves past `n` hold `e`
 *   and never receive a tag other than `id`, so `e = LLONG_MAX`
 *   cannot overflow.
 * Source: AtCoder Library lazy_segtree (CC0).
 * Status: stress-tested against brute force
 */
#pragma once

struct LazySegTree {
  typedef ll S; // node: min of the range
  typedef ll F; // tag: add to every element
  static constexpr S e = LLONG_MAX;
  static constexpr F id = 0;
  S op(S a, S b) { return min(a, b); }
  S mapping(F f, S s) { return s + f; }
  F composition(F f, F g) { return f + g; }
  int n, lg = 0, size;
  vector<S> d;
  vector<F> lz;
  LazySegTree(const vector<S>& a) : n(sz(a)) {
    while (1 << lg < n) lg++;
    size = 1 << lg;
    d.assign(2 * size, e), lz.assign(size, id);
    rep(i, 0, n) d[size + i] = a[i];
    for (int i = size; --i;) pull(i);
  }
  void pull(int k) { d[k] = op(d[2 * k], d[2 * k + 1]); }
  void put(int k, F f) {
    d[k] = mapping(f, d[k]);
    if (k < size) lz[k] = composition(f, lz[k]);
  }
  void push(int k) {
    put(2 * k, lz[k]), put(2 * k + 1, lz[k]), lz[k] = id;
  }
  // ancestors of leaves l, r - 1 not fully inside [l, r)
  void push_ends(int l, int r) {
    for (int i = lg; i; i--) {
      if (l % (1 << i)) push(l >> i);
      if (r % (1 << i)) push((r - 1) >> i);
    }
  }
  void pull_ends(int l, int r) {
    rep(i, 1, lg + 1) {
      if (l % (1 << i)) pull(l >> i);
      if (r % (1 << i)) pull((r - 1) >> i);
    }
  }
  void set(int i, S v) {
    i += size, push_ends(i, i + 1);
    d[i] = v, pull_ends(i, i + 1);
  }
  S get(int i) {
    i += size, push_ends(i, i + 1);
    return d[i];
  }
  S query_all() { return d[1]; }
  S query(int l, int r) {
    S a = e, b = e;
    l += size, r += size, push_ends(l, r);
    for (; l < r; l /= 2, r /= 2) {
      if (l & 1) a = op(a, d[l++]);
      if (r & 1) b = op(d[--r], b);
    }
    return op(a, b);
  }
  void apply(int l, int r, F f) {
    l += size, r += size, push_ends(l, r);
    for (int a = l, b = r; a < b; a /= 2, b /= 2) {
      if (a & 1) put(a++, f);
      if (b & 1) put(--b, f);
    }
    pull_ends(l, r);
  }
  void apply(int i, F f) { apply(i, i + 1, f); }
  int max_right(int l, auto pred) {
    if (l == n) return n;
    S acc = e;
    l += size, push_ends(l, l + 1);
    do {
      while (l % 2 == 0) l /= 2;
      if (!pred(op(acc, d[l]))) {
        while (l < size) {
          push(l), l *= 2;
          if (pred(op(acc, d[l]))) acc = op(acc, d[l++]);
        }
        return l - size;
      }
      acc = op(acc, d[l++]);
    } while ((l & -l) != l);
    return n;
  }
  int min_left(int r, auto pred) {
    if (r == 0) return 0;
    S acc = e;
    r += size, push_ends(r - 1, r);
    do {
      r--;
      while (r > 1 && r % 2) r /= 2;
      if (!pred(op(d[r], acc))) {
        while (r < size) {
          push(r), r = 2 * r + 1;
          if (pred(op(d[r], acc))) acc = op(d[r--], acc);
        }
        return r + 1 - size;
      }
      acc = op(d[r], acc);
    } while ((r & -r) != r);
    return 0;
  }
};
