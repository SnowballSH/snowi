#include "../testlib.hpp"
#include "../../content/numtheory/sieve.hpp"

vi trial_factor(int x) {
  vi f;
  for (int d = 2; d * d <= x; d++)
    for (; x % d == 0; x /= d) f.push_back(d);
  if (x > 1) f.push_back(x);
  return f;
}

// phi and mu straight from their definitions.
int brute_phi(int x) {
  int c = 0;
  rep(k, 1, x + 1) c += gcd(k, x) == 1;
  return c;
}
int brute_mu(int x) {
  vi f = trial_factor(x);
  if (adjacent_find(all(f)) != f.end()) return 0;
  return sz(f) % 2 ? -1 : 1;
}

void check_sieve(int n) {
  Sieve s(n);
  CHECK_EQ(sz(s.spf), n + 1);
  vi primes;
  rep(x, 2, n + 1)
    if (sz(trial_factor(x)) == 1) primes.push_back(x);
  CHECK(s.primes == primes);
  CHECK_EQ(s.spf[0] + s.spf[1] + s.phi[0] + s.mu[0], 0);
  rep(x, 1, n + 1) {
    vi f = trial_factor(x);
    CHECK(s.factorize(x) == f);
    if (x > 1) CHECK_EQ(s.spf[x], f[0]);
    CHECK_EQ(s.phi[x], brute_phi(x));
    CHECK_EQ(s.mu[x], brute_mu(x));
  }
}

int main() {
  rep(n, 1, 40) check_sieve(n);
  check_sieve(1000), check_sieve(1024), check_sieve(997);
  Sieve u(1e3); // Usage
  CHECK(u.factorize(360) == vi({2, 2, 2, 3, 3, 5}));

  // Example: coprime pairs (i, j) in [1, N] x [1, M].
  Sieve s(300);
  rep(it, 0, iters(300)) {
    int N = (int)rnd(1, 300), M = (int)rnd(1, 60);
    ll cnt = 0, want = 0;
    rep(d, 1, min(N, M) + 1) cnt += s.mu[d] * ll(N / d) * (M / d);
    rep(i, 1, N + 1)
      rep(j, 1, M + 1) want += gcd(i, j) == 1;
    CHECK_EQ(cnt, want);
  }

  elapsed();
  int n = 2'000'000;
  Sieve big(n);
  CHECK_EQ(sz(big.primes), 148933);
  // sum_d mu(d) (n/d) = 1 and sum_d phi(d) (n/d) = n (n+1) / 2
  ll mu_sum = 0, phi_sum = 0;
  rep(d, 1, n + 1) {
    mu_sum += big.mu[d] * ll(n / d);
    phi_sum += big.phi[d] * ll(n / d);
  }
  CHECK_EQ(mu_sum, 1LL);
  CHECK_EQ(phi_sum, ll(n) * (n + 1) / 2);
  CHECK(big.factorize(1999993) == vi{1999993});
  CHECK(big.factorize(1 << 20) == vi(20, 2));
  fprintf(stderr, "  sieve 2e6: %.2fs\n", elapsed());
  puts("OK numtheory/sieve");
}
