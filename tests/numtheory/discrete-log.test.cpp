#include "../testlib.hpp"
#include "../../content/numtheory/discrete-log.hpp"
#include "../../content/numtheory/euclid.hpp"

// a^x mod m takes at most m distinct values before it cycles, so the
// first x hitting b, if any, is below m.
ll brute_log(ll a, ll b, ll m) {
  ll cur = 1 % m;
  rep(x, 0, m + 1) {
    if (cur == b % m) return x;
    cur = cur * a % m;
  }
  return -1;
}

ll pow_mod(ll b, ll e, ll m) {
  ll r = 1 % m;
  for (b %= m; e; e >>= 1, b = b * b % m)
    if (e & 1) r = r * b % m;
  return r;
}

bool is_prime_slow(ll n) {
  if (n < 2) return false;
  for (ll d = 2; d * d <= n; d++)
    if (n % d == 0) return false;
  return true;
}

// For a large m: the answer must solve the equation, and no smaller
// x from a known solution may be skipped.
void check_large(ll a, ll x, ll m) {
  ll b = pow_mod(a, x, m), got = discrete_log(a, b, m);
  CHECK(got >= 0 && got <= x);
  CHECK_EQ(pow_mod(a, got, m), b);
}

int main() {
  rep(m, 1, 70)
    rep(a, 0, m)
      rep(b, 0, m) CHECK_EQ(discrete_log(a, b, m), brute_log(a, b, m));
  rep(it, 0, iters(3000)) {
    ll m = rnd(1, rnd(0, 1) ? 3000 : 300);
    if (rnd(0, 2) == 0) m = 1LL << rnd(0, 11); // prime powers
    else if (rnd(0, 2) == 0) m = rnd(1, 30) * rnd(1, 30) * 4;
    ll a = rnd(0, m - 1), b = rnd(0, m - 1);
    if (rnd(0, 1)) b = pow_mod(a, rnd(0, 2 * m), m);
    CHECK_EQ(discrete_log(a, b, m), brute_log(a, b, m));
  }
  // Unreduced inputs.
  CHECK_EQ(discrete_log(17 + 3, 17 * 5 + 13, 17), 4LL);
  CHECK_EQ(discrete_log(3, 13, 17), 4LL);
  CHECK_EQ(discrete_log(2, 0, 8), 3LL);
  CHECK_EQ(discrete_log(2, 1, 8), 0LL);
  CHECK_EQ(discrete_log(2, 3, 8), -1LL);
  CHECK_EQ(discrete_log(2, 8, 12), 3LL);
  CHECK_EQ(discrete_log(0, 0, 1), 0LL);
  CHECK_EQ(discrete_log(0, 1, 7), 0LL);
  CHECK_EQ(discrete_log(0, 0, 7), 1LL);
  CHECK_EQ(discrete_log(5, 0, 7), -1LL);
  // Negative inputs are reduced modulo m first.
  rep(m, 1, 40)
    rep(a, -2 * m, 0)
      rep(b, -2 * m, 1) {
        ll ra = (a % m + m) % m, rb = (b % m + m) % m;
        CHECK_EQ(discrete_log(a, b, m), brute_log(ra, rb, m));
      }
  CHECK_EQ(discrete_log(-1, 2, 3), 1LL);
  CHECK_EQ(discrete_log(-2, 2, 4), 1LL);
  CHECK_EQ(discrete_log(-14, -4, 17), 4LL);
  CHECK_EQ(discrete_log(LLONG_MIN, 0, 2), 1LL);
  // Example: the order of a unit.
  rep(m, 1, 200)
    rep(a, 1, m) {
      if (gcd(a, m) != 1) continue;
      ll ord = 1 + discrete_log(a, inv_mod(a, m), m), want = 1;
      while (pow_mod(a, want, m) != 1 % m) want++;
      CHECK_EQ(ord, want);
    }
  // Large moduli, up to the documented 3e9: primes (a^x is a long
  // cycle), smooth numbers, powers of two and random composites.
  for (ll m : {2999999929LL, 3000000000LL, 2147483648LL, 1000000007LL,
               998244353LL, 999999937LL * 3, 735134400LL}) {
    rep(it, 0, 3) check_large(rnd(0, m - 1), rnd(0, m), m);
    ll a = rnd(0, m - 1), b = rnd(0, m - 1);
    CHECK_EQ(discrete_log(a - m, b - 2 * m, m), discrete_log(a, b, m));
    CHECK_EQ(discrete_log(m - 1, 1, m), 0LL);
  }
  CHECK_EQ(discrete_log(5, 1, 1000000007), 0LL);
  CHECK_EQ(discrete_log(5, 5, 1000000007), 1LL);
  CHECK_EQ(discrete_log(5, 1000000006, 1000000007), 500000003LL);
  // 5 is a primitive root of 1e9 + 7: every x < p - 1 is the answer.
  rep(it, 0, 3) {
    ll p = 1000000007, x = rnd(0, p - 2);
    CHECK_EQ(discrete_log(5, pow_mod(5, x, p), p), x);
  }
  // No solution: a quadratic non-residue is not a power of a square.
  CHECK_EQ(discrete_log(25, 5, 1000000007), -1LL);

  elapsed();
  int queries = 0;
  rep(it, 0, 60) {
    ll m = rnd(5e8, 1e9);
    if (it % 2) // primes: the slowest case, no early exit
      while (!is_prime_slow(m)) m++;
    check_large(rnd(1, m - 1), rnd(0, m - 2), m), queries++;
  }
  fprintf(stderr, "  %d logs mod ~1e9: %.2fs\n", queries, elapsed());
  puts("OK numtheory/discrete-log");
}
