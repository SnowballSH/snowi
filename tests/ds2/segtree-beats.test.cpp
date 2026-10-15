#include "../testlib.hpp"
#include "../../content/ds2/segtree-beats.hpp"

const ll INF = Beats::INF;

// Checks every field of a query node against the subarray b.
void check_node(const Beats::Node& nd, vector<ll> b) {
  CHECK_EQ(nd.sum, accumulate(all(b), 0LL));
  CHECK_EQ(nd.len, sz(b));
  if (b.empty()) {
    CHECK_EQ(nd.mx, -INF);
    CHECK_EQ(nd.mn, INF);
    return;
  }
  sort(all(b));
  ll mx = b.back(), mn = b[0];
  CHECK_EQ(nd.mx, mx);
  CHECK_EQ(nd.mn, mn);
  CHECK_EQ(nd.cmx, (int)count(all(b), mx));
  CHECK_EQ(nd.cmn, (int)count(all(b), mn));
  auto below = lower_bound(all(b), mx);
  CHECK_EQ(nd.mx2, below == b.begin() ? -INF : *prev(below));
  auto above = upper_bound(all(b), mn);
  CHECK_EQ(nd.mn2, above == b.end() ? INF : *above);
}

// Random operations on [l, r); x within [-c, c], adds within +-d.
void stress(int rounds, int maxn, int ops, ll c, ll d) {
  rep(it, 0, iters(rounds)) {
    int n = (int)rnd(0, maxn);
    vector<ll> a = rnd_vec(n, -c, c);
    Beats st(a);
    rep(op, 0, ops) {
      int l = (int)rnd(0, n), r = (int)rnd(0, n), kind = (int)rnd(0, 4);
      if (l > r) swap(l, r);
      ll x = rnd(-c, c);
      if (kind == 0) {
        st.chmin(l, r, x);
        rep(i, l, r) a[i] = min(a[i], x);
      } else if (kind == 1) {
        st.chmax(l, r, x);
        rep(i, l, r) a[i] = max(a[i], x);
      } else if (kind == 2) {
        x = rnd(-d, d);
        st.add(l, r, x);
        rep(i, l, r) a[i] += x;
      } else if (kind == 3) { // Tweak: assign
        st.chmin(l, r, x), st.chmax(l, r, x);
        rep(i, l, r) a[i] = x;
      } else {
        check_node(st.query(l, r), vector<ll>(a.begin() + l, a.begin() + r));
      }
    }
    rep(i, 0, n) CHECK_EQ(st.query(i, i + 1).sum, a[i]);
    check_node(st.query(0, n), a);
  }
}

int main() {
  { // the Example
    Beats st(vector<ll>{5, 1, 4, 1, 3});
    st.chmin(0, 5, 3), st.add(1, 4, 2), st.chmax(0, 2, 4);
    auto nd = st.query(0, 5);
    CHECK_EQ(nd.sum, 19LL);
    CHECK_EQ(nd.mx, 5LL);
    CHECK_EQ(nd.mn, 3LL);
    CHECK_EQ(nd.cmx, 1);
    CHECK_EQ(nd.mx2, 4LL);
    st.chmin(2, 4, 0), st.chmax(2, 4, 0);
    check_node(st.query(0, 5), {4, 4, 0, 0, 3});
  }
  stress(3000, 12, 50, 3, 2); // few distinct values
  stress(2000, 16, 50, 10, 10);
  stress(300, 64, 200, 1000, 300);
  // extremes: |values| near 2e17, sums up to 1.6e18
  stress(500, 8, 50, 200'000'000'000'000'000, 1'000'000'000'000'000);
  // Medium against brute force, then a large run for speed.
  elapsed();
  stress(1, 2000, 20'000, 1'000'000'000, 1'000'000);
  int n = 100'000;
  vector<ll> a = rnd_vec(n, 0, 1'000'000'000);
  Beats st(a);
  ll total = 0;
  rep(op, 0, n) {
    int l = (int)rnd(0, n - 1), r = (int)rnd(l + 1, n), kind = (int)rnd(0, 3);
    ll x = rnd(0, 1'000'000'000);
    if (kind == 0) st.chmin(l, r, x);
    else if (kind == 1) st.chmax(l, r, x);
    else if (kind == 2) st.add(l, r, x - 500'000'000);
    else total += st.query(l, r).sum;
  }
  ll s = 0;
  rep(i, 0, n) s += st.query(i, i + 1).sum;
  CHECK_EQ(st.query(0, n).sum, s);
  CHECK(total != 0);
  fprintf(stderr, "  large inputs: %.2fs\n", elapsed());
  puts("OK ds2/segtree-beats");
}
