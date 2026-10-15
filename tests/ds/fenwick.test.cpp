#include "../testlib.hpp"
#include "../../content/ds/fenwick.hpp"

int brute_lower_bound(const vector<ll>& a, ll s) {
  if (s <= 0) return -1;
  ll pre = 0;
  rep(i, 0, sz(a))
    if ((pre += a[i]) >= s) return i;
  return sz(a);
}

int main() {
  rep(it, 0, iters(3000)) {
    int n = (int)rnd(1, 20);
    vector<ll> a(n);
    Fenwick fw(n);
    rep(op, 0, 30) {
      int i = (int)rnd(0, n - 1);
      ll x = rnd(0, 10);
      a[i] += x, fw.add(i, x);
      int l = (int)rnd(0, n), r = (int)rnd(l, n);
      CHECK_EQ(fw.query(r) - fw.query(l),
               accumulate(a.begin() + l, a.begin() + r, 0LL));
      ll s = rnd(-2, accumulate(all(a), 0LL) + 2);
      CHECK_EQ(fw.lower_bound(s), brute_lower_bound(a, s));
    }
  }
  int n = 1 << 17;
  Fenwick fw(n);
  rep(i, 0, n) fw.add(i, 1);
  CHECK_EQ(fw.query(n), (ll)n);
  CHECK_EQ(fw.lower_bound(n), n - 1);
  CHECK_EQ(fw.lower_bound(n + 1), n);
  puts("OK ds/fenwick");
}
