#include "../testlib.hpp"
#include "../../content/numtheory/factor.hpp"

typedef unsigned long long ull;
typedef unsigned __int128 u128;

// Independent oracle: Miller-Rabin with 24 random bases.
bool oracle_prime(ull n) {
  if (n < 4) return n >= 2;
  if (n % 2 == 0) return false;
  ull d = n - 1;
  int s = 0;
  while (d % 2 == 0) d /= 2, s++;
  rep(it, 0, 24) {
    ull a = (ull)rnd(2, (ll)min<ull>(n - 2, LLONG_MAX)), x = 1;
    for (ull b = a, e = d; e; e /= 2, b = (u128)b * b % n)
      if (e % 2) x = (u128)x * b % n;
    bool ok = x == 1 || x == n - 1;
    for (int i = 1; i < s && !ok; i++) x = (u128)x * x % n, ok = x == n - 1;
    if (!ok) return false;
  }
  return true;
}

vector<ll> trial_factor(ll x) {
  vector<ll> f;
  for (ll d = 2; d * d <= x; d++)
    for (; x % d == 0; x /= d) f.push_back(d);
  if (x > 1) f.push_back(x);
  return f;
}

ll random_prime(ll lo, ll hi) {
  for (;;)
    if (ll p = rnd(lo, hi); oracle_prime(p)) return p;
}

ll slow_pow(ll b, ll e, ll m) {
  ll r = 1 % m;
  rep(i, 0, e) r = (ll)((u128)r * b % m);
  return r;
}

void check_factor(ll n) {
  vector<ll> f = factor(n);
  CHECK(is_sorted(all(f)));
  ll prod = 1;
  for (ll p : f) {
    CHECK(oracle_prime(p) && is_prime(p));
    CHECK((u128)prod * p <= (u128)n);
    prod *= p;
  }
  CHECK_EQ(prod, n);
}

int main() {
  // mulmod and powmod, reused by other snippets: any 1 <= m < 2^63.
  CHECK_EQ(powmod(5, 0, 1), 0LL);
  rep(m, 1, 40)
    rep(b, 0, 3 * m) {
      rep(a, 0, 3 * m) CHECK_EQ(mulmod(a, b, m), a * b % m);
      rep(e, 0, 12) CHECK_EQ(powmod(b, e, m), slow_pow(b, e, m));
    }
  rep(it, 0, iters(20000)) {
    ll m = rnd(1, LLONG_MAX >> rnd(0, 62)), b = rnd(0, LLONG_MAX);
    ll e1 = rnd(0, LLONG_MAX / 2), e2 = rnd(0, LLONG_MAX / 2);
    ll x = powmod(b, e1, m), y = powmod(b, e2, m);
    CHECK(0 <= x && x < m && 0 <= y && y < m);
    CHECK_EQ(powmod(b, e1 + e2, m), mulmod(x, y, m));
    CHECK_EQ(mulmod(b, 1, m), b % m);
  }

  vector<bool> composite(200'001);
  rep(i, 2, 200'001) {
    CHECK_EQ(is_prime(i), !composite[i]);
    for (int j = 2 * i; j <= 200'000; j += i) composite[j] = true;
  }
  CHECK(is_prime(1'000'000'007) &&
        factor(360) == vector<ll>({2, 2, 2, 3, 3, 5}));
  CHECK(!is_prime(0) && !is_prime(1) && !is_prime(-7));
  CHECK(!is_prime(LLONG_MIN));
  // Strong pseudoprimes to several small bases, Carmichael numbers.
  for (ll n : {2047LL, 3215031751LL, 4759123141LL, 1122004669633LL,
               3474749660383LL, 341550071728321LL, 3825123056546413051LL, 561LL,
               41041LL, 825265LL, 321197185LL, 5394826801LL}) {
    CHECK(!is_prime(n));
    check_factor(n);
  }
  // Strong pseudoprimes to all bases but one, in order of the missing
  // base: dropping any base of A must fail here.
  for (ll n :
       {1921077350011LL, 1411807385341LL, 443538368977861LL, 4341937413061LL,
        5517315475561LL, 3933464309633LL, 107528788110061LL})
    CHECK(!is_prime(n));
  rep(it, 0, iters(20000)) {
    ll n = rnd(1, LLONG_MAX >> rnd(0, 62));
    CHECK_EQ(is_prime(n), oracle_prime(n));
  }
  rep(it, 0, iters(300)) {
    ll p = random_prime(2, LLONG_MAX >> rnd(0, 61));
    CHECK(is_prime(p));
    if (p <= INT_MAX) CHECK(!is_prime(p * p));
  }

  rep(n, 1, 3000) CHECK(factor(n) == trial_factor(n));
  rep(it, 0, iters(1000)) {
    ll n = rnd(1, (ll)1e9);
    CHECK(factor(n) == trial_factor(n));
  }
  rep(it, 0, iters(2000)) check_factor(rnd(1, LLONG_MAX >> rnd(0, 62)));
  for (ll n :
       {LLONG_MAX, LLONG_MAX - 1, LLONG_MAX - 24, 1LL << 62,
        4052555153018976267LL, // 3^39
        (1LL << 61) - 1, 999999999999999989LL, (ll)1e18, 897612484786617600LL})
    check_factor(n);
  ll m31 = INT_MAX; // 2^31 - 1 is prime
  CHECK(factor(m31 * m31) == vector<ll>({m31, m31}));
  rep(it, 0, iters(30)) { // prime squares and cubes
    ll p = random_prime(2e9, 3e9), q = random_prime(1e6, 2e6);
    CHECK(factor(p * p) == vector<ll>({p, p}));
    CHECK(factor(q * q * q) == vector<ll>({q, q, q}));
  }

  // Example: all divisors from the factorisation.
  auto divisors = [](ll n) {
    vector<ll> divs = {1};
    ll last = 0;
    for (int from = 0, to; ll p : factor(n)) {
      if (p != last) from = 0, last = p;
      to = sz(divs);
      rep(i, from, to) divs.push_back(divs[i] * p);
      from = to;
    }
    return divs;
  };
  rep(n, 1, 2000) {
    vector<ll> d = divisors(n), want;
    rep(k, 1, n + 1)
      if (n % k == 0) want.push_back(k);
    sort(all(d));
    CHECK(d == want);
  }
  CHECK_EQ(sz(divisors(897612484786617600LL)), 103680);

  elapsed();
  rep(it, 0, 60) {
    ll p = random_prime(5e8, 1e9), q = random_prime(1e9, 2e9);
    CHECK(factor(p * q) == vector<ll>({p, q}));
  }
  double t = elapsed();
  fprintf(stderr, "  60 semiprimes near 1e18: %.2fs\n", t);
  // Worst case allowed: two primes near sqrt(2^63) = 3.04e9.
  rep(it, 0, 30) {
    ll p = random_prime(3e9, 3037000499), q = random_prime(3e9, 3e9 + 3e7);
    CHECK(factor(p * q) == vector<ll>({min(p, q), max(p, q)}));
  }
  t = elapsed() - t;
  fprintf(stderr, "  30 semiprimes near 2^63: %.2fs\n", t);
  puts("OK numtheory/factor");
}
