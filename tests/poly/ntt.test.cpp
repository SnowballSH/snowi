#include "../testlib.hpp"
#include "../../content/poly/ntt.hpp"

vl brute_conv(const vl& a, const vl& b) {
  if (a.empty() || b.empty()) return {};
  vl c(sz(a) + sz(b) - 1);
  rep(i, 0, sz(a))
    rep(j, 0, sz(b)) c[i + j] = (c[i + j] + a[i] * b[j]) % MOD;
  return c;
}

// ntt(a)[k] = a(w^k) with w = 3^((MOD - 1) / n), natural order.
vl brute_dft(const vl& a) {
  int n = sz(a);
  ll w = modpow(3, (MOD - 1) / n);
  vl res(n);
  rep(k, 0, n) {
    ll x = modpow(w, k), p = 1;
    for (ll c : a) res[k] = (res[k] + c * p) % MOD, p = p * x % MOD;
  }
  return res;
}

int main() {
  rep(it, 0, iters(1500)) {
    int n = (int)rnd(0, it % 3 ? 50 : 120);
    int m = (int)rnd(0, it % 3 ? 50 : 120);
    ll hi = it % 5 ? MOD - 1 : 3;
    vl a = rnd_vec(n, 0, hi), b = rnd_vec(m, 0, hi);
    if (it % 7 == 0) fill(all(a), MOD - 1), fill(all(b), MOD - 1);
    CHECK(conv(a, b) == brute_conv(a, b));
    // Correlation from the header: sum_i a_i b_(i + d) at
    // c[sz(a) - 1 + d] of conv(reversed a, b).
    vl c = conv(vl(a.rbegin(), a.rend()), b);
    rep(d, 0, n ? m - n + 1 : 0) {
      ll want = 0;
      rep(i, 0, n) want = (want + a[i] * b[i + d]) % MOD;
      CHECK_EQ(c[n - 1 + d], want);
    }
  }
  rep(lg, 0, 8) {
    vl a = rnd_vec(1 << lg, 0, MOD - 1), f = a;
    ntt(f);
    CHECK(f == brute_dft(a));
    int n = sz(a);
    ntt(f), reverse(f.begin() + 1, f.end());
    for (ll& x : f) x = x * modinv(n) % MOD;
    CHECK(f == a);
  }
  // Large: all-ones gives c[k] = min(k + 1, n, m, n + m - 1 - k).
  int n = 1 << 16, m = (1 << 16) + 3;
  vl c = conv(vl(n, 1), vl(m, 1));
  CHECK_EQ(sz(c), n + m - 1);
  rep(k, 0, sz(c)) CHECK_EQ(c[k], (ll)min({k + 1, n, m, n + m - 1 - k}));
  vl a = rnd_vec(n, 0, MOD - 1), b = rnd_vec(n, 0, MOD - 1);
  c = conv(a, b);
  rep(t, 0, 20) {
    int k = (int)rnd(0, 2 * n - 2);
    ll want = 0;
    rep(i, max(0, k - n + 1), min(k, n - 1) + 1)
      want = (want + a[i] * b[k - i]) % MOD;
    CHECK_EQ(c[k], want);
  }
  puts("OK poly/ntt");
}
