#include "../testlib.hpp"
#include "../../content/combo/binom.hpp"

const int P = 80; // Pascal's triangle size

ll pw(ll b, ll e, ll m) {
  ll r = 1 % m;
  for (b %= m; e; e >>= 1, b = b * b % m)
    if (e & 1) r = r * b % m;
  return r;
}

// C(n, k) mod a prime p from the product formula, with the factors
// of p counted separately (Legendre), for small k and huge n.
ll product_binom(ll n, ll k, ll p) {
  if (k < 0 || k > n) return 0;
  k = min(k, n - k);
  ll num = 1, den = 1, v = 0;
  rep(i, 0, k) {
    ll a = n - i, b = i + 1;
    for (; a % p == 0; a /= p) v++;
    for (; b % p == 0; b /= p) v--;
    num = num * (a % p) % p, den = den * (b % p) % p;
  }
  return v > 0 ? 0 : num * pw(den, p - 2, p) % p;
}

vector<vector<ll>> pascal(ll m) {
  vector<vector<ll>> c(P, vector<ll>(P));
  rep(n, 0, P) {
    c[n][0] = 1 % m;
    rep(k, 1, n + 1) c[n][k] = (c[n - 1][k - 1] + c[n - 1][k]) % m;
  }
  return c;
}

// T is the storage type: Binom<int> must match Binom<ll>.
template<class T> void test_small() {
  auto pas = pascal(MOD);
  vector<ll> cat(P / 2, 0);
  cat[0] = 1;
  rep(m, 1, P / 2)
    rep(i, 0, m) cat[m] = (cat[m] + cat[i] * cat[m - 1 - i]) % MOD;
  for (int size : {0, 1, 2, 5, P - 1}) {
    Binom<T> bc(size);
    rep(n, -3, size + 1)
      rep(k, -3, size + 3) {
        bool in = 0 <= k && k <= n;
        CHECK_EQ(bc.C(n, k), in ? pas[n][k] : 0);
        ll want = in ? 1 : 0;
        if (in)
          rep(i, 0, k) want = want * (n - i) % MOD;
        CHECK_EQ(bc.perm(n, k), want);
      }
    rep(k, 0, size / 2 + 1) CHECK_EQ(bc.catalan(k), cat[k]);
    rep(i, 1, size + 1) CHECK_EQ(bc.inv(i) * i % MOD, 1LL);
    rep(i, 0, size + 1) CHECK_EQ(1LL * bc.fact[i] * bc.ifact[i] % MOD, 1LL);
    CHECK_EQ(bc.multinom({}), 1LL);
    rep(it, 0, iters(300)) {
      vi ks;
      int s = 0;
      while (true) {
        int k = (int)rnd(0, 6);
        if (s + k > size || rnd(0, 4) == 0) break;
        ks.push_back(k), s += k;
      }
      ll want = 1; // choose the slots of each part in turn
      s = 0;
      for (int k : ks) s += k, want = want * pas[s][k] % MOD;
      CHECK_EQ(bc.multinom(ks), want);
    }
  }
}

template<class T> void test_large(int n) {
  Binom<T> bc(n);
  rep(it, 0, iters(100'000)) {
    int a = (int)rnd(1, n), k = (int)rnd(0, min(a, 25));
    if (rnd(0, 1)) k = a - k;
    CHECK_EQ(bc.C(a, k), product_binom(a, k, MOD));
    CHECK_EQ(bc.C(a, k), (bc.C(a - 1, k - 1) + bc.C(a - 1, k)) % MOD);
  }
  // Catalan via the ballot formula C(2k, k) / (k + 1).
  rep(it, 0, iters(1000)) {
    int k = (int)rnd(0, n / 2);
    CHECK_EQ(bc.catalan(k), bc.C(2 * k, k) * modinv(k + 1) % MOD);
  }
}

void test_lucas() {
  for (ll p : {2, 3, 5, 7, 13, 97, 1009}) {
    Lucas lu(p);
    auto pas = pascal(p);
    rep(i, 0, p) CHECK_EQ(lu.fact[i] * lu.ifact[i] % p, 1LL);
    rep(n, -2, P)
      rep(k, -2, P + 2) {
        bool in = 0 <= k && k <= n;
        CHECK_EQ(lu.C(n, k), in ? pas[n][k] : 0);
      }
    rep(it, 0, iters(3000)) {
      ll n = rnd(0, (ll)1e18), k = rnd(0, 40);
      if (rnd(0, 1)) k = n - k;
      if (rnd(0, 3) == 0) n = rnd(0, 3000), k = rnd(-5, n + 5);
      CHECK_EQ(lu.C(n, k), product_binom(n, k, p));
    }
  }
  ll p = 1'000'003;
  Lucas lu(p);
  rep(it, 0, iters(20'000)) {
    ll n = rnd(0, (ll)1e18), k = rnd(0, 30);
    if (rnd(0, 1)) k = n - k;
    CHECK_EQ(lu.C(n, k), product_binom(n, k, p));
  }
  // C(p - 1, k) = (-1)^k and C(p, k) = 0 for 0 < k < p.
  rep(it, 0, iters(1000)) {
    ll k = rnd(1, p - 1);
    CHECK_EQ(lu.C(p - 1, k), k % 2 ? p - 1 : 1);
    CHECK_EQ(lu.C(p, k), 0LL);
  }
}

// Regression: int storage once overflowed int * int products
// (C(1e7, 1) came out negative).
void test_int_storage() {
  int n = 10'000'000;
  Binom<int> bc(n);
  CHECK_EQ(bc.C(n, 1), (ll)n);
  CHECK_EQ(bc.perm(n, 2), (ll)n * (n - 1) % MOD);
  CHECK_EQ(bc.inv(n) * n % MOD, 1LL);
  CHECK_EQ(bc.multinom({n / 2, n - n / 2}), bc.C(n, n / 2));
}

int main() {
  Binom bc(5); // the default storage is ll
  static_assert(is_same_v<decltype(bc.fact), vector<ll>>);
  test_small<ll>();
  test_small<int>();
  test_large<ll>(1'000'000);
  test_large<int>(2'000'000);
  test_int_storage();
  test_lucas();
  puts("OK combo/binom");
}
