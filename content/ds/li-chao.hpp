/**
 * Title: Li Chao tree
 * Description: Minimum of lines $k x + m$ at the points of a fixed
 *   sorted list `xs`. `add_segment` adds a line that exists only
 *   for `x` in `[l, r)`. `query(x)` returns `INF` if no line covers
 *   `x`.
 * Precondition: `xs` is sorted and non-empty; `query(x)` needs `x`
 *   in `xs`; `l` and `r` may be any values. Every `k * x` and
 *   `k * x + m` with `x` in `xs` fits in `ll` and is below `INF`
 *   (e.g. $|k|, |x| <= 10^9$ and $|m| <= 10^18$).
 * Time: $O(log n)$ per line or query, $O(log^2 n)$ per segment.
 * Usage:
 *   LiChao lc(xs); // vector<ll>, the query points
 *   lc.add_line(k, m); lc.add_segment(k, m, l, r);
 *   ll y = lc.query(x);
 * Tweak: Online `x` in a large range: create children on demand
 *   over the integers `[lo, hi)` and use `mid` itself instead of
 *   `xs[mid]`.
 * Example:
 *   // dp[j] = min over i < j of dp[i] + (h[j] - h[i])^2 + C
 *   LiChao lc(sorted_unique_h);
 *   lc.add_line(-2 * h[0], h[0] * h[0]); // dp[0] = 0
 *   rep(j, 1, n) {
 *     dp[j] = lc.query(h[j]) + h[j] * h[j] + C;
 *     lc.add_line(-2 * h[j], dp[j] + h[j] * h[j]);
 *   }
 * Notes: Think of it when a DP transition expands to
 *   $f(i) g(j) + a(i) + b(j)$: each `i` is a line with slope $f(i)$
 *   and intercept $a(i)$, queried at $x = g(j)$. Signals: squared
 *   differences $(h_j - h_i)^2$, products of prefix sums, "for
 *   every query `c`, the best value minus `c` times size". Slopes
 *   and queries may come in any order, unlike the monotone hull,
 *   and segments handle lines valid on a range of `x` only. Add
 *   line `j` after computing `dp[j]`, never before. Insertion
 *   evaluates lines at points that are never queried, so keep
 *   `xs` to the real query points to avoid overflow. Each node
 *   keeps the line that wins at its middle point; the loser can
 *   win on one side only and moves down into that child. A query
 *   takes the minimum along the leaf-to-root path. On a tree,
 *   root-path queries need only rollback: `add_line` changes one
 *   node per level, so in a DFS log `(v, K[v], M[v])` before each
 *   change and undo the log on leaving the vertex. Subtree
 *   queries merge Li Chao trees.
 * Source: Li Chao tree as in cp-algorithms; Packet 06 (UCup).
 * Status: stress-tested against brute force
 */
#pragma once

struct LiChao {
  static constexpr ll INF = LLONG_MAX;
  int size = 1;
  vector<ll> xs, K, M; // node v holds the line K[v] x + M[v]
  LiChao(const vector<ll>& x) : xs(x) {
    while (size < sz(x)) size *= 2;
    xs.resize(size, x.back());
    K.assign(2 * size, 0), M.assign(2 * size, INF);
  }
  int idx(ll x) { return lower_bound(all(xs), x) - begin(xs); }
  // Node v covers indices [lo, hi); the line exists on [l, r).
  void add(ll k, ll m, int l, int r, int v, int lo, int hi) {
    if (r <= lo || hi <= l) return;
    int mid = (lo + hi) / 2;
    if (l > lo || r < hi) {
      add(k, m, l, r, 2 * v, lo, mid);
      return add(k, m, l, r, 2 * v + 1, mid, hi);
    }
    auto below = [&](int i) {
      return k * xs[i] + m < K[v] * xs[i] + M[v];
    };
    bool a = below(lo);
    if (a == below(hi - 1)) { // one line wins everywhere
      if (a) K[v] = k, M[v] = m;
      return;
    }
    if (below(mid)) swap(k, K[v]), swap(m, M[v]), a = !a;
    if (a) add(k, m, l, r, 2 * v, lo, mid);
    else add(k, m, l, r, 2 * v + 1, mid, hi);
  }
  void add_line(ll k, ll m) { add(k, m, 0, size, 1, 0, size); }
  void add_segment(ll k, ll m, ll l, ll r) { // x in [l, r)
    add(k, m, idx(l), idx(r), 1, 0, size);
  }
  ll query(ll x) {
    ll y = INF;
    for (int v = idx(x) + size; v; v /= 2)
      y = min(y, K[v] * x + M[v]);
    return y;
  }
};
