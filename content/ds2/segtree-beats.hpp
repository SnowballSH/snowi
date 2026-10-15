/**
 * Title: Segment tree beats
 * Description: Range chmin ($a_i = min(a_i, x)$), range chmax
 *   ($a_i = max(a_i, x)$) and range add, with range sum, max and
 *   min. A node keeps its max, strict second max and count of the
 *   max, and the same for the min. Chmin stops at a node where only
 *   the max class changes and fixes the sum from the count.
 * Precondition: 0-indexed, half-open ranges. Values, `x` and range
 *   sums stay below $4 dot 10^18$ in absolute value.
 * Time: $O((n + q) log^2 n)$ amortized with range add, and
 *   $O((n + q) log n)$ with only chmin, chmax and queries. Memory:
 *   $4 n$ nodes of 64 bytes.
 * Usage:
 *   Beats st(a); // vector<ll> a
 *   st.chmin(l, r, x); st.chmax(l, r, x); st.add(l, r, x);
 *   auto nd = st.query(l, r); // nd.sum, nd.mx, nd.mn
 * Tweak: Assign `x` on `[l, r)`: `chmin(l, r, x)`, `chmax(l, r,
 *   x)`. For only chmin and sum, drop the min fields and
 *   `put_max`. More node data (say a count of odd values) goes
 *   into `Node` and `merge`, and every `put_` must update it in
 *   $O(1)$ from the fields it changes.
 * Example:
 *   Beats st(vector<ll>{5, 1, 4, 1, 3});
 *   st.chmin(0, 5, 3); // 3 1 3 1 3
 *   st.add(1, 4, 2);   // 3 3 5 3 3
 *   st.chmax(0, 2, 4); // 4 4 5 3 3
 *   auto nd = st.query(0, 5); // sum 19, mx 5, mn 3
 *   // nd.cmx = 1 copy of the max, nd.mx2 = 4 (strict 2nd max)
 *   st.chmin(2, 4, 0), st.chmax(2, 4, 0); // 4 4 0 0 3
 * Notes: Think of it when a range update clamps values ("set all
 *   above `x` to `x`", "everybody gets at least `x`") and sums or
 *   counts are asked, which a lazy tag cannot do. Without adds,
 *   the bound charges each extra recursion to the number of
 *   distinct values in a node: it drops when the max class merges
 *   with the second, and a chmin or chmax raises it by at most one
 *   in each of the $O(log n)$ partly covered nodes. One range add
 *   can raise it by half a node's length, so with adds Ji's proof
 *   needs a finer potential (over "tag classes") and gets
 *   $O(log^2 n)$; random tests visit about $4 log_2 n$ nodes per
 *   operation, like a plain segment tree. The same "recurse while
 *   something changes" recursion handles `a_i %= x`, `a_i /= x`,
 *   `a_i = floor(sqrt(a_i))` and `a_i &= x` on values $>= 0$: stop
 *   where nothing changes (`mx < x` for `%=`) and update leaves;
 *   each value shrinks $O(log V)$ times. Pitfalls: the second max
 *   is strict and $-infinity$ at a leaf, adds must not shift the
 *   sentinels, and a node with one or two distinct values links
 *   the max and min fields (`put_min`, `put_max`). Members named
 *   `max` or `min` would hide `std::max` inside the struct.
 * Source: Ji Ruyi, segment tree beats (2016).
 * Status: stress-tested against brute force
 */
#pragma once

struct Beats {
  static constexpr ll INF = LLONG_MAX; // never shifted
  struct Node {
    ll sum = 0, lz = 0;
    ll mx = -INF, mx2 = -INF, mn = INF, mn2 = INF;
    int cmx = 0, cmn = 0, len = 0;
  };
  int n;
  vector<Node> t;
  Beats(const vector<ll>& a) : n(sz(a)), t(4 * n + 4) {
    if (n) build(1, 0, n, a);
  }
  static Node merge(const Node& a, const Node& b) {
    Node c;
    c.sum = a.sum + b.sum, c.len = a.len + b.len;
    c.mx = max(a.mx, b.mx), c.mn = min(a.mn, b.mn);
    c.cmx = (a.mx == c.mx) * a.cmx + (b.mx == c.mx) * b.cmx;
    c.cmn = (a.mn == c.mn) * a.cmn + (b.mn == c.mn) * b.cmn;
    c.mx2 = max(a.mx == c.mx ? a.mx2 : a.mx,
                b.mx == c.mx ? b.mx2 : b.mx);
    c.mn2 = min(a.mn == c.mn ? a.mn2 : a.mn,
                b.mn == c.mn ? b.mn2 : b.mn);
    return c;
  }
  void build(int k, int lo, int hi, const vector<ll>& a) {
    if (hi - lo == 1) {
      t[k].sum = t[k].mx = t[k].mn = a[lo];
      t[k].cmx = t[k].cmn = t[k].len = 1;
      return;
    }
    int mid = (lo + hi) / 2;
    build(2 * k, lo, mid, a), build(2 * k + 1, mid, hi, a);
    t[k] = merge(t[2 * k], t[2 * k + 1]);
  }
  void put_add(int k, ll x) {
    Node& c = t[k];
    c.sum += x * c.len, c.mx += x, c.mn += x, c.lz += x;
    if (c.mx2 != -INF) c.mx2 += x;
    if (c.mn2 != INF) c.mn2 += x;
  }
  void put_min(int k, ll x) { // needs mx2 < x < mx
    Node& c = t[k];
    c.sum += (x - c.mx) * c.cmx;
    if (c.mn == c.mx) c.mn = x;
    else if (c.mn2 == c.mx) c.mn2 = x;
    c.mx = x;
  }
  void put_max(int k, ll x) { // needs mn < x < mn2
    Node& c = t[k];
    c.sum += (x - c.mn) * c.cmn;
    if (c.mx == c.mn) c.mx = x;
    else if (c.mx2 == c.mn) c.mx2 = x;
    c.mn = x;
  }
  void push(int k) {
    for (int c : {2 * k, 2 * k + 1}) {
      if (t[k].lz) put_add(c, t[k].lz);
      if (t[c].mx > t[k].mx) put_min(c, t[k].mx);
      if (t[c].mn < t[k].mn) put_max(c, t[k].mn);
    }
    t[k].lz = 0;
  }
  // op 0: chmin, 1: chmax, 2: add
  void upd(int k, int lo, int hi, int l, int r, ll x, int op) {
    if (r <= lo || hi <= l || (op == 0 && t[k].mx <= x) ||
        (op == 1 && t[k].mn >= x))
      return;
    if (l <= lo && hi <= r &&
        (op == 2 || (op == 0 && t[k].mx2 < x) ||
         (op == 1 && t[k].mn2 > x))) {
      if (op == 0) put_min(k, x);
      else if (op == 1) put_max(k, x);
      else put_add(k, x);
      return;
    }
    push(k);
    int mid = (lo + hi) / 2;
    upd(2 * k, lo, mid, l, r, x, op);
    upd(2 * k + 1, mid, hi, l, r, x, op);
    t[k] = merge(t[2 * k], t[2 * k + 1]);
  }
  void chmin(int l, int r, ll x) { upd(1, 0, n, l, r, x, 0); }
  void chmax(int l, int r, ll x) { upd(1, 0, n, l, r, x, 1); }
  void add(int l, int r, ll x) { upd(1, 0, n, l, r, x, 2); }
  Node query(int k, int lo, int hi, int l, int r) {
    if (r <= lo || hi <= l) return Node();
    if (l <= lo && hi <= r) return t[k];
    push(k);
    int mid = (lo + hi) / 2;
    return merge(query(2 * k, lo, mid, l, r),
                 query(2 * k + 1, mid, hi, l, r));
  }
  Node query(int l, int r) { return query(1, 0, n, l, r); }
};
