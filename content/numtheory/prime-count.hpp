/**
 * Title: Prime counting
 * Description: Lucy_Hedgehog's sieve: `pc(v)` is $pi(v)$, the number
 *   of primes $<= v$, for every $v = floor(n\/i)$ and every
 *   $v <= sqrt(n)$.
 * Precondition: $n >= 0$. Other values of `v` give wrong answers.
 * Time: $O(n^(3\/4) \/ log n)$, memory $O(sqrt(n))$: 0.05 s for
 *   $n = 10^10$, 0.3 s for $10^11$, 1.4 s for $10^12$.
 * Usage:
 *   PrimeCount pc(n);
 *   ll a = pc(n), b = pc(n / 7);
 * Tweak: Sum of primes: use `__int128` (or values mod a prime),
 *   start from $i(i + 1)\/2 - 1$ and $v(v + 1)\/2 - 1$, and subtract
 *   $p (S(floor(v\/p)) - S(p - 1))$. Any completely multiplicative
 *   weight works the same way, e.g. one table per power $p^k$.
 * Example:
 *   // n <= 1e11 with exactly 4 divisors: p^3 or p q with p < q
 *   ll ans = 0;
 *   for (ll p = 2; p <= pc.r; p++)
 *     if (pc(p) > pc(p - 1))
 *       ans += (p * p * p <= n) + pc(n / p) - pc(p);
 * Notes: $pi(10^10) = 455052511$, $pi(10^11) = 4118054813$,
 *   $pi(10^12) = 37607912018$. After sieving with the primes below
 *   `p`, `S(v)` counts the numbers in $[2, v]$ that are prime or
 *   have no prime factor below `p`, so only $v >= p^2$ change.
 *   Cues: one query about primes with $n$ up to $10^11$, sums like
 *   $sum_p f(floor(n\/p))$, numbers with few prime factors (as in
 *   the example).
 *   Min_25 and Du's sieves: _Number theory facts_.
 * Source: Lucy_Hedgehog's post in the Project Euler problem 10
 *   forum; min_25's blog.
 * Status: stress-tested against a sieve
 */
#pragma once

struct PrimeCount {
  ll n, r;
  vector<ll> lo, hi; // lo[v] = pi(v), hi[i] = pi(n / i)
  PrimeCount(ll n_) : n(n_), r(sqrtl(n_)) {
    while (r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    lo.resize(r + 1), hi.resize(r + 1);
    rep(i, 1, r + 1) lo[i] = i - 1, hi[i] = n / i - 1;
    for (ll p = 2; p <= r; p++) {
      if (lo[p] == lo[p - 1]) continue; // p is not prime
      ll cnt = lo[p - 1], q = p * p, m = min(r, n / q);
      for (ll i = 1, d = p; i <= m; i++, d += p)
        hi[i] -= (d <= r ? hi[d] : lo[n / d]) - cnt;
      for (ll v = r; v >= q; v--) lo[v] -= lo[v / p] - cnt;
    }
  }
  ll operator()(ll v) { return v <= r ? lo[v] : hi[n / v]; }
};
