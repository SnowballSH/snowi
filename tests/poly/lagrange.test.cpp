#include "../testlib.hpp"
#include "../../content/poly/lagrange.hpp"

ll horner(const vl& c, ll x) {
  x = (x % MOD + MOD) % MOD;
  ll r = 0;
  for (int k = sz(c); k--;) r = (r * x + c[k]) % MOD;
  return r;
}

int main() {
  rep(it, 0, iters(3000)) {
    int n = (int)rnd(1, it % 4 ? 12 : 40);
    vl c = rnd_vec(n, 0, it % 5 ? MOD - 1 : 3);
    if (it % 9 == 0) fill(all(c), MOD - 1);
    vl y(n);
    rep(i, 0, n) y[i] = horner(c, i);
    // Nodes themselves, small, huge and negative x.
    rep(t, 0, 4) {
      ll x = t == 0   ? rnd(0, n - 1)
             : t == 1 ? rnd(0, MOD - 1)
             : t == 2 ? rnd(MOD, (ll)4e18)
                      : -rnd(0, (ll)4e18);
      CHECK_EQ(lagrange(y, x), horner(c, x));
    }
    // A degree below n through n distinct points, some at 0, MOD - 1.
    set<ll> seen;
    vl xs, ys;
    if (it % 3 == 0) xs.push_back(0), seen.insert(0);
    if (it % 5 == 0) xs.push_back(MOD - 1), seen.insert(MOD - 1);
    xs.resize(min(sz(xs), n));
    while (sz(xs) < n) {
      ll x = it % 2 ? rnd(0, MOD - 1) : rnd(0, 3 * n);
      if (seen.insert(x).second) xs.push_back(x);
    }
    shuffle(all(xs), test_rng);
    for (ll x : xs) ys.push_back(horner(c, x));
    CHECK(interpolate(xs, ys) == c);
  }
  CHECK(interpolate({}, {}).empty());
  CHECK_EQ(lagrange({7}, 123456789012LL), 7LL);
  // Usage line: sum of squares.
  CHECK_EQ(lagrange({0, 1, 5, 14}, 10), 385LL);
  // Example: sum_(i = 1)^n i^k against the brute sum.
  rep(k, 0, 8) {
    vl y(k + 2);
    rep(i, 1, k + 2) y[i] = (y[i - 1] + modpow(i, k)) % MOD;
    ll brute = 0;
    rep(m, 0, 60) {
      if (m) brute = (brute + modpow(m, k)) % MOD;
      CHECK_EQ(lagrange(y, m), brute);
    }
  }
  // Large: f(x) = (x + 5)^d with d = 2e5; interpolate n = 1500.
  int d = 200000;
  vl y(d + 1);
  rep(i, 0, d + 1) y[i] = modpow(i + 5, d);
  rep(t, 0, 3) {
    ll x = rnd(0, (ll)1e18);
    CHECK_EQ(lagrange(y, x), modpow((x + 5) % MOD, d));
  }
  int n = 1500;
  vl c = rnd_vec(n, 0, MOD - 1), xs(n), ys(n);
  rep(i, 0, n) xs[i] = (ll)i * i % MOD, ys[i] = horner(c, xs[i]);
  CHECK(interpolate(xs, ys) == c);
  puts("OK poly/lagrange");
}
