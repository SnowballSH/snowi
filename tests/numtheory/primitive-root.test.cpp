#include "../testlib.hpp"
#include "../../content/numtheory/primitive-root.hpp"

bool is_prime_slow(ll n) {
  if (n < 2) return false;
  for (ll d = 2; d * d <= n; d++)
    if (n % d == 0) return false;
  return true;
}

int brute_order(ll a, ll p) {
  ll cur = a % p;
  int e = 1;
  for (; cur != 1; e++) cur = cur * a % p;
  return e;
}

vector<ll> distinct_primes(ll n) { // trial division
  vector<ll> f;
  for (ll d = 2; d * d <= n; d++)
    if (n % d == 0) {
      f.push_back(d);
      while (n % d == 0) n /= d;
    }
  if (n > 1) f.push_back(n);
  return f;
}

vector<ll> rho_primes(ll n) { // factor is tested on its own
  vector<ll> f = factor(n);
  f.erase(unique(all(f)), f.end());
  return f;
}

bool generates(ll g, ll p, const vector<ll>& qs) {
  for (ll q : qs)
    if (powmod(g, (p - 1) / q, p) == 1) return false;
  return true;
}

// Large p: the answer generates and no smaller g does, judged with
// the prime factors qs of p - 1.
void check_large(ll p, const vector<ll>& qs) {
  ll g = primitive_root(p);
  CHECK(generates(g, p, qs));
  rep(h, 2, g) CHECK(!generates(h, p, qs));
}

ll order_example(ll a, ll p) {
  vector<ll> f = factor(p - 1);
  ll e = p - 1;
  for (ll q : f)
    if (powmod(a, e / q, p) == 1) e /= q;
  return e;
}

int main() {
  rep(p, 2, 6000) {
    if (!is_prime_slow(p)) continue;
    int g = 1;
    while (brute_order(g, p) != p - 1) g++;
    CHECK_EQ(primitive_root(p), (ll)g);
    if (p < 600)
      rep(a, 1, p) CHECK_EQ(order_example(a, p), (ll)brute_order(a, p));
  }
  CHECK_EQ(primitive_root(998244353), 3LL);
  CHECK_EQ(primitive_root(1'000'000'007), 5LL);
  CHECK_EQ(primitive_root((1LL << 61) - 1), 37LL);
  CHECK_EQ(primitive_root(1'000'000'009), 13LL);
  CHECK_EQ(primitive_root(90441961), 113LL); // record below 1e8
  // k-th root of unity from the example.
  ll p = 998244353, w = powmod(3, (p - 1) >> 20, p);
  CHECK(powmod(w, 1 << 20, p) == 1 && powmod(w, 1 << 19, p) == p - 1);
  // Up to 1e12: independent trial-division factorization of p - 1.
  rep(it, 0, iters(40)) {
    ll q = rnd(2, (ll)1e12);
    while (!is_prime(q)) q++;
    check_large(q, distinct_primes(q - 1));
  }
  // Up to 2^63: primes of the form k 2^s + 1 and random ones.
  for (ll q : {4179340454199820289LL, 9223372036854775783LL,
               1000000000000000003LL, 2305843009213693951LL})
    check_large(q, rho_primes(q - 1));
  elapsed();
  rep(it, 0, 300) {
    ll q = rnd(1, LLONG_MAX);
    while (!is_prime(q)) q = rnd(1, LLONG_MAX);
    check_large(q, rho_primes(q - 1));
  }
  fprintf(stderr, "  300 random 63-bit primes: %.2fs\n", elapsed());
  puts("OK numtheory/primitive-root");
}
