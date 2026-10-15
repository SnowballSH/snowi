#include "../testlib.hpp"
#include "../../content/ds/sparse-table.hpp"

int main() {
  rep(it, 0, iters(3000)) {
    int n = (int)rnd(0, 40);
    bool extreme = rnd(0, 1);
    vector<ll> a =
        extreme ? rnd_vec(n, LLONG_MIN, LLONG_MAX) : rnd_vec(n, -5, 5);
    SparseTable st(a);
    vector<pii> b(n);
    rep(i, 0, n) b[i] = {(int)rnd(0, 3), i};
    SparseTable sb(b);
    rep(l, 0, n) {
      ll m = LLONG_MAX;
      pii mb = {INT_MAX, 0};
      rep(r, l + 1, n + 1) {
        m = min(m, a[r - 1]), mb = min(mb, b[r - 1]);
        CHECK_EQ(st.query(l, r), m);
        CHECK(sb.query(l, r) == mb);
      }
    }
  }
  // Large: 1e5 queries against block minima of size 256.
  int n = 100'000, B = 256;
  vector<ll> a = rnd_vec(n, -1e18, 1e18);
  vector<ll> block(n / B + 1, LLONG_MAX);
  rep(i, 0, n) block[i / B] = min(block[i / B], a[i]);
  SparseTable st(a);
  rep(it, 0, n) {
    int l = (int)rnd(0, n - 1), r = (int)rnd(l + 1, n);
    ll m = LLONG_MAX;
    int i = l;
    for (; i < r && i % B; i++) m = min(m, a[i]);
    for (; i + B <= r; i += B) m = min(m, block[i / B]);
    for (; i < r; i++) m = min(m, a[i]);
    CHECK_EQ(st.query(l, r), m);
  }
  puts("OK ds/sparse-table");
}
