"""
Title: Linear sieve
Description: For all $i <= n$: smallest prime factor `spf[i]`,
  Euler's totient `phi[i]` and Möbius `mu[i]`, plus the list of
  primes. `factorize(x)` lists the prime factors of $x <= n$ in
  increasing order, with multiplicity.
Precondition: $n >= 1$. Index 0 holds zeros, `spf[1] = 0`.
Time: $O(n)$; `factorize` $O(log x)$. PyPy: 0.1 s for $n = 10^6$;
  about 1 s and 340 MB for $n = 10^7$.
Usage:
  s = Sieve(10**6)
  s.primes[k]; s.spf[x]; s.phi[x]; s.mu[x]
  f = s.factorize(360)  # [2, 2, 2, 3, 3, 5]
Tweak: Drop `phi` or `mu` when unused: each list costs 8 bytes per
  number. Another multiplicative $f$: set $f(p)$ for a new prime,
  $f(i p) = f(i) f(p)$ when $p < op("spf")(i)$; when
  $p = op("spf")(i)$ you need $f$ of one more power of $p$, so also
  store the power of the smallest prime in `i`.
Example:
  >>> s = Sieve(30)
  >>> s.primes[:6], s.phi[12], s.mu[30], s.factorize(24)
  ([2, 3, 5, 7, 11, 13], 4, -1, [2, 2, 2, 3])
  >>> # only primes, n = 1e7 in 0.3 s: bytearray Eratosthenes
  >>> from math import isqrt
  >>> n = 50
  >>> is_p = bytearray([1]) * (n + 1); is_p[0] = is_p[1] = 0
  >>> for i in range(2, isqrt(n) + 1):
  ...   if is_p[i]:
  ...     is_p[i * i :: i] = bytes(len(range(i * i, n + 1, i)))
  >>> [i for i in range(n + 1) if is_p[i]][-3:]
  [41, 43, 47]
Notes: $phi(n)$ counts $1 <= k <= n$ coprime to `n`; $mu$ is the
  Möbius function: $sum_(d | n) mu(d) = [n = 1]$ and
  $sum_(d | n) phi(d) = n$. A gcd condition in a sum turns into
  $[gcd = 1] = sum_(d | gcd) mu(d)$: the pairs $i <= N$, $j <= M$
  with $gcd(i, j) = 1$ number $sum_d mu(d) floor(N/d) floor(M/d)$.
  There are 664579 primes below $10^7$ and 78498 below $10^6$.
Source: Linear (Euler) sieve, UCup Packet 12.
Status: stress-tested against trial division
"""


class Sieve:
  def __init__(self, n):
    self.spf = spf = [0] * (n + 1)
    self.phi = phi = [0] * (n + 1)
    self.mu = mu = [0] * (n + 1)
    self.primes = primes = []
    phi[1] = mu[1] = 1
    for i in range(2, n + 1):
      if not spf[i]:
        spf[i], phi[i], mu[i] = i, i - 1, -1
        primes.append(i)
      s, f, m = spf[i], phi[i], -mu[i]
      top = min(s, n // i)
      for p in primes:
        if p > top:
          break
        spf[i * p] = p
        if p == s:
          phi[i * p] = f * p
        else:
          phi[i * p], mu[i * p] = f * (p - 1), m

  def factorize(self, x):
    spf, f = self.spf, []
    while x > 1:
      f.append(spf[x])
      x //= spf[x]
    return f
