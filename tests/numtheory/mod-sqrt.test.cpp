#include "../testlib.hpp"
#include "../../content/numtheory/mod-sqrt.hpp"

typedef unsigned __int128 u128;

ll sq(ll x, ll p) { return (ll)((u128)x * x % p); }

bool is_prime_slow(ll n) {
  if (n < 2) return false;
  for (ll d = 2; d * d <= n; d++)
    if (n % d == 0) return false;
  return true;
}

// A large prime p: squares must give back a root of the same square,
// and a non-square (Euler's criterion) times a square gives -1.
void check_prime(ll p) {
  ll y = rnd(0, p - 1), r = mod_sqrt(sq(y, p), p);
  CHECK(r == min(y, p - y) || (y == 0 && r == 0));
  ll non = 2;
  while (powmod(non, (p - 1) / 2, p) != p - 1) non++;
  if (y) CHECK_EQ(mod_sqrt(mulmod(non, sq(y, p), p), p), -1LL);
  CHECK_EQ(mod_sqrt(sq(y, p) - p, p), r); // negative input
}

int main() {
  // Exhaustive for small primes: the smaller root, or -1.
  rep(p, 2, 1200) {
    if (!is_prime_slow(p)) continue;
    vi root(p, -1);
    for (int x = p - 1; x >= 0; x--) root[x * x % p] = x;
    rep(a, 0, p) CHECK_EQ(mod_sqrt(a, p), (ll)root[a]);
    CHECK_EQ(mod_sqrt(p + 1, p), 1LL);
    CHECK_EQ(mod_sqrt(5 * p, p), 0LL);
    rep(a, 1, 2 * p + 1) // negative inputs
      CHECK_EQ(mod_sqrt(-a, p), (ll)root[(2 * p - a) % p]);
  }
  // Primes p = k 2^s + 1 with large s: the slowest case.
  for (ll p : {998244353LL, 469762049LL, 2013265921LL, 167772161LL,
               4179340454199820289LL, (1LL << 61) - 1, 1000000007LL, 3LL, 5LL,
               17LL, 65537LL, 9223372036854775783LL})
    rep(it, 0, 300) check_prime(p);
  rep(it, 0, iters(2000)) {
    ll p = rnd(2, rnd(0, 1) ? 1000000 : 3000000000LL);
    while (!is_prime_slow(p)) p++;
    check_prime(p);
  }
  CHECK_EQ(mod_sqrt(2, 7), 3LL);
  CHECK_EQ(mod_sqrt(3, 7), -1LL);
  CHECK_EQ(mod_sqrt(0, 5), 0LL);
  CHECK_EQ(mod_sqrt(2, 998244353), 116195171LL);
  CHECK_EQ(mod_sqrt(5, 2), 1LL);
  CHECK_EQ(mod_sqrt(LLONG_MAX, 2), 1LL);
  CHECK_EQ(mod_sqrt(-1, 13), 5LL);
  CHECK_EQ(mod_sqrt(-3, 7), 2LL);
  CHECK_EQ(mod_sqrt(-4, 5), 1LL);
  CHECK_EQ(mod_sqrt(LLONG_MIN, 2), 0LL);
  CHECK_EQ(mod_sqrt(LLONG_MIN + 1, 2), 1LL);
  ll big = 9223372036854775783LL; // the largest prime below 2^63
  CHECK_EQ(mod_sqrt(LLONG_MIN, big), mod_sqrt(LLONG_MIN + big, big));

  elapsed();
  ll p = 998244353, found = 0;
  rep(i, 0, 100000) {
    ll a = rnd(0, p - 1), r = mod_sqrt(a, p);
    if (r < 0) continue;
    found++;
    CHECK_EQ(sq(r, p), a);
  }
  CHECK(found > 49000 && found < 51000);
  fprintf(stderr, "  1e5 roots mod 998244353: %.2fs\n", elapsed());
  puts("OK numtheory/mod-sqrt");
}
