#include "../testlib.hpp"
#include "../../content/numtheory/modint.hpp"
#include "../../content/numtheory/euclid.hpp"

ll norm(ll v, ll m = MOD) { return (v % m + m) % m; }

// The runtime-modulus Tweak applied to the snippet, as documented.
namespace dyn {
ll md = 1;
struct Mint {
  using M = Mint;
  ll x;
  Mint(ll v = 0) : x((v %= md) < 0 ? v + md : v) {}
  M& operator+=(M o) {
    if ((x += o.x) >= md) x -= md;
    return *this;
  }
  M& operator-=(M o) {
    if ((x -= o.x) < 0) x += md;
    return *this;
  }
  M& operator*=(M o) { return *this = x * o.x; }
  M& operator/=(M o) { return *this *= o.inv(); }
  friend M operator+(M a, M b) { return a += b; }
  friend M operator-(M a, M b) { return a -= b; }
  friend M operator*(M a, M b) { return a *= b; }
  friend M operator/(M a, M b) { return a /= b; }
  M operator-() const { return M() - *this; }
  M operator^(ll e) const {
    if (e < 0) return inv() ^ -e;
    M r = 1, b = *this;
    for (; e; e >>= 1, b *= b)
      if (e & 1) r *= b;
    return r;
  }
  M inv() const { return inv_mod(x, md); }
  bool operator==(const M&) const = default;
};
// Built before main sets md: must not divide by zero.
Mint early = 7;
} // namespace dyn

void check_dynamic(ll m) {
  dyn::md = m;
  rep(it, 0, iters(3000)) {
    const ll big = 1'000'000'000'000;
    ll u = rnd(-big, big), v = rnd(-big, big);
    ll a = norm(u, m), b = norm(v, m);
    dyn::Mint x = u, y = v;
    CHECK_EQ((x + y).x, (a + b) % m);
    CHECK_EQ((x - y).x, (a - b + m) % m);
    CHECK_EQ((x * y).x, a * b % m);
    if (gcd(b, m) != 1) continue;
    CHECK_EQ((x / y * y).x, a);
    CHECK_EQ(((y ^ -3) * (y ^ 3)).x, 1 % m);
  }
}

ll slow_pow(ll b, ll e) {
  ll r = 1;
  rep(i, 0, e) r = r * b % MOD;
  return r;
}

constexpr auto fact = [] {
  array<Mint, 21> f{1};
  rep(i, 1, 21) f[i] = f[i - 1] * i;
  return f;
}();
static_assert(fact[5] == 120 && fact[0].x == 1);
static_assert((Mint(2) ^ 10) == 1024 && Mint(-1).x == MOD - 1);
static_assert(Mint(3) / 3 == 1 && -Mint(0) == 0);
static_assert((Mint(2) ^ -1) * 2 == 1 && (Mint(3) ^ -2) * 9 == 1);

ll pick() {
  int kind = (int)rnd(0, 3);
  if (kind == 0) return rnd(-5, 5);
  if (kind == 1) return MOD - rnd(1, 3);
  if (kind == 2) return rnd(LLONG_MIN, LLONG_MAX);
  return rnd(0, MOD - 1);
}

int main() {
  rep(it, 0, iters(20000)) {
    ll u = pick(), v = pick(), a = norm(u), b = norm(v);
    Mint x = u, y = v;
    CHECK_EQ(x.x, a);
    CHECK_EQ((x + y).x, (a + b) % MOD);
    CHECK_EQ((x - y).x, (a - b + MOD) % MOD);
    CHECK_EQ((x * y).x, a * b % MOD);
    CHECK_EQ((-x).x, (MOD - a) % MOD);
    CHECK_EQ(x == y, a == b);
    Mint z = x;
    z += y, z -= y, z *= y;
    CHECK_EQ(z.x, a * b % MOD);
    if (b) {
      CHECK_EQ((x / y * y).x, a);
      CHECK_EQ((y.inv() * y).x, 1LL);
    }
    ll e = rnd(0, 40);
    CHECK_EQ((x ^ e).x, slow_pow(a, e));
    ll e1 = rnd(0, LLONG_MAX / 2), e2 = rnd(0, LLONG_MAX / 2);
    CHECK_EQ((x ^ e1) * (x ^ e2), x ^ (e1 + e2));
    if (a) {
      CHECK_EQ(x ^ -1, x.inv());
      CHECK_EQ((x ^ -e) * (x ^ e), Mint(1));
      CHECK_EQ((x ^ -e1) * (x ^ e1), Mint(1));
      CHECK_EQ(x ^ -LLONG_MAX, x.inv() ^ LLONG_MAX);
    }
  }
  CHECK_EQ(dyn::early.x, 0LL); // reduced modulo the default md = 1
  for (ll m : {1LL, 2LL, 10LL, 12LL, 360LL, 1'000'000'000LL, MOD})
    check_dynamic(m);
  CHECK_EQ(Mint(0).inv(), Mint(0));
  CHECK_EQ(Mint(0) ^ 0, Mint(1));
  CHECK_EQ(Mint(LLONG_MIN).x, norm(LLONG_MIN));
  CHECK_EQ(fact[20].x, (ll)(2432902008176640000LL % MOD));
  stringstream ss;
  ss << Mint(-2) << ' ' << Mint(7) * 3;
  CHECK_EQ(ss.str(), to_string(MOD - 2) + " 21");
  Mint a = 5, b = (a ^ 3) / 2 - 1; // Usage
  CHECK(b * 2 == 123 && b * b.inv() == 1);
  Mint s = 0;
  for (int k : {1, 2, 3}) s += Mint(k).inv();
  CHECK(s * 6 == 11);

  elapsed();
  Mint acc = 1, sum = 0;
  rep(i, 1, 200'001) acc *= i, sum += acc / i;
  Mint want = 0, f = 1;
  rep(i, 1, 200'001) want += f, f *= i;
  CHECK_EQ(sum, want);
  fprintf(stderr, "  2e5 mul + div: %.2fs\n", elapsed());
  puts("OK numtheory/modint");
}
