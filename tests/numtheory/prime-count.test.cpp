#include "../testlib.hpp"
#include "../../content/numtheory/prime-count.hpp"

const int LIM = 2'000'000;
vi pi_table; // pi_table[v] = number of primes <= v, by Eratosthenes

// Every value the structure promises: v <= sqrt(n) and v = n / i.
void check_all(ll n) {
  PrimeCount pc(n);
  CHECK(pc.r * pc.r <= n && (pc.r + 1) * (pc.r + 1) > n);
  rep(v, 0, pc.r + 1) CHECK_EQ(pc(v), (ll)pi_table[v]);
  for (ll i = 1; i <= n; i = n / (n / i) + 1)
    CHECK_EQ(pc(n / i), (ll)pi_table[n / i]);
}

int main() {
  vector<bool> comp(LIM + 1);
  pi_table.assign(LIM + 1, 0);
  rep(i, 2, LIM + 1) {
    pi_table[i] = pi_table[i - 1] + !comp[i];
    if (!comp[i])
      for (ll j = (ll)i * i; j <= LIM; j += i) comp[j] = true;
  }
  rep(n, 0, 3000) check_all(n);
  rep(it, 0, iters(300)) check_all(rnd(0, LIM));
  for (ll n : {LIM, LIM - 1, 1 << 20, 999'999, 1'000'000, 1'000'001})
    check_all(n);
  rep(r, 2, 1415) { // squares and their neighbours: sqrt rounding
    for (ll n : {(ll)r * r - 1, (ll)r * r, (ll)r * r + 1})
      if (n <= LIM) check_all(n);
  }

  // Example: numbers with exactly 4 divisors (p^3 or p q, p < q).
  vi divs(20001);
  for (int d = 1; d <= 20000; d++)
    for (int k = d; k <= 20000; k += d) divs[k]++;
  int brute = 0;
  rep(n, 1, 20001) {
    brute += divs[n] == 4;
    if (n % 97 && n > 200) continue;
    PrimeCount pc(n);
    ll ans = 0;
    for (ll p = 2; p <= pc.r; p++)
      if (pc(p) > pc(p - 1)) ans += (p * p * p <= n) + pc(n / p) - pc(p);
    CHECK_EQ(ans, (ll)brute);
  }

  elapsed();
  PrimeCount big(10'000'000'000LL);
  CHECK_EQ(big(10'000'000'000LL), 455052511LL);
  CHECK_EQ(big(10'000'000'000LL / 10), 50847534LL); // pi(1e9)
  CHECK_EQ(big(100'000), 9592LL);
  CHECK_EQ(PrimeCount(1LL << 32)(1LL << 32), 203280221LL);
  double t = elapsed();
  fprintf(stderr, "  pi(1e10), pi(2^32): %.2fs\n", t);
  PrimeCount huge(100'000'000'000LL);
  CHECK_EQ(huge(100'000'000'000LL), 4118054813LL);
  fprintf(stderr, "  pi(1e11): %.2fs\n", elapsed() - t);
  puts("OK numtheory/prime-count");
}
