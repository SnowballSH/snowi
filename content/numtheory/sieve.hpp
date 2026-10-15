/**
 * Title: Linear sieve
 * Description: For all $i <= n$: smallest prime factor `spf[i]`,
 *   Euler's totient `phi[i]` and Möbius `mu[i]`, plus the list of
 *   primes. `factorize(x)` lists the prime factors of $x <= n$ in
 *   increasing order, with multiplicity.
 * Precondition: $n >= 1$. Index 0 holds zeros, `spf[1] = 0`.
 * Time: $O(n)$, 0.2 s for $n = 10^7$; `factorize` $O(log x)$.
 *   Memory: 12 bytes per number, 120 MB for $n = 10^7$; drop the
 *   tables you do not need.
 * Usage:
 *   Sieve s(1e7);
 *   s.primes[k]; s.spf[x]; s.phi[x]; s.mu[x];
 *   vi f = s.factorize(360); // {2, 2, 2, 3, 3, 5}
 * Tweak: Another multiplicative $f$: set $f(p)$ for a new prime,
 *   $f(i p) = f(i) f(p)$ when $p < op("spf")(i)$; when
 *   $p = op("spf")(i)$ you need $f$ of one more power of $p$, so
 *   also store the power of the smallest prime in `i`.
 * Example:
 *   // pairs 1 <= i <= N, 1 <= j <= M with gcd(i, j) = 1:
 *   ll cnt = 0; // sieve up to min(N, M)
 *   rep(d, 1, min(N, M) + 1) cnt += s.mu[d] * ll(N / d) * (M / d);
 * Notes: $phi(n)$ counts $1 <= k <= n$ coprime to `n`; $mu$ is the
 *   Möbius function: $sum_(d | n) mu(d) = [n = 1]$ and
 *   $sum_(d | n) phi(d) = n$. A gcd condition in a sum turns into
 *   $[gcd = 1] = sum_(d | gcd) mu(d)$. There are 664579 primes below
 *   $10^7$ and 78498 below $10^6$. For only the primes, a
 *   `vector<bool>` Eratosthenes is about 3 times faster (0.08 s for
 *   $n = 10^7$) and needs 1 bit per number.
 * Source: Linear (Euler) sieve, UCup Packet 12.
 * Status: stress-tested against trial division
 */
#pragma once

struct Sieve {
  vi primes, spf, phi, mu;
  Sieve(int n) : spf(n + 1), phi(n + 1), mu(n + 1) {
    phi[1] = mu[1] = 1;
    rep(i, 2, n + 1) {
      if (!spf[i]) spf[i] = i, phi[i] = i - 1, mu[i] = -1;
      if (spf[i] == i) primes.push_back(i);
      for (int p : primes) {
        if (p > spf[i] || (ll)i * p > n) break;
        spf[i * p] = p;
        phi[i * p] = phi[i] * (p == spf[i] ? p : p - 1);
        mu[i * p] = p == spf[i] ? 0 : -mu[i];
      }
    }
  }
  vi factorize(int x) {
    vi f;
    for (; x > 1; x /= spf[x]) f.push_back(spf[x]);
    return f;
  }
};
