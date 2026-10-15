"""
Title: Prime counting
Description: Lucy_Hedgehog's sieve. After an $O(n^(3\/4))$ build,
  `pc(v)` is $pi(v)$, the number of primes $<= v$, for $v = n$ and
  every $v = floor(n\/i)$, and for every $v <= sqrt(n)$.
Precondition: $n >= 0$. Other values of `v` give wrong answers.
Time: $O(n^(3\/4) \/ log n)$, memory $O(sqrt(n))$. PyPy: 0.2 s for
  $n = 10^10$, 0.9 s for $10^11$.
Usage:
  pc = PrimeCount(n)
  a, b = pc(n), pc(n // 7)
Tweak: Sum of primes: start from $i(i + 1)\/2 - 1$ and
  $v(v + 1)\/2 - 1$, and subtract $p (S(floor(v\/p)) - S(p - 1))$;
  reduce modulo a prime if allowed, since big ints are slow. Any
  completely multiplicative weight works the same way, e.g. one
  table per power $p^k$.
Example:
  >>> pc = PrimeCount(10**6)
  >>> pc(10**6), pc(10**6 // 3), pc(1000)
  (78498, 28665, 168)
  >>> n, ans = 100, 0  # exactly 4 divisors: p^3 or p q, p < q
  >>> pc = PrimeCount(n)
  >>> for p in range(2, pc.r + 1):
  ...   if pc(p) > pc(p - 1):
  ...     ans += (p**3 <= n) + pc(n // p) - pc(p)
  >>> ans
  32
Notes: $pi(10^10) = 455052511$, $pi(10^11) = 4118054813$,
  $pi(10^12) = 37607912018$. After sieving with the primes below
  `p`, `S(v)` counts the numbers in $[2, v]$ that are prime or
  have no prime factor below `p`, so only $v >= p^2$ change.
  Cues: one query about primes with $n$ up to $10^11$, sums like
  $sum_p f(floor(n\/p))$, numbers with few prime factors (as in
  the example). Min_25 sieve, for a multiplicative `f` with $f(p)$
  a polynomial in `p`: get $G(v) = sum_(p <= v) f(p)$ from one
  table per monomial; then $S(v, j)$, the sum of $f(i)$ over
  $2 <= i <= v$ with no prime factor below $p_j$, is
  $G(v) - G(p_(j - 1)) + sum_(k >= j, p_k^2 <= v)
  sum_(e >= 1, p_k^(e + 1) <= v) (f(p_k^e) S(floor(v\/p_k^e),
  k + 1) + f(p_k^(e + 1)))$, and $sum_(i <= n) f(i) = S(n, 1) +
  f(1)$. Plain recursion is fine: its depth is $O(log n)$. For
  $Phi(n) = sum_(i <= n) phi(i)$ and $M(n) = sum_(i <= n) mu(i)$,
  Du's sieve is shorter: $Phi(n) = n(n + 1)\/2 - sum_(d = 2)^n
  Phi(floor(n\/d))$ and $M(n) = 1 - sum_(d = 2)^n M(floor(n\/d))$,
  summed over floor blocks and memoized on the values
  $floor(n\/i)$, with a sieve for the arguments up to $n^(2\/3)$:
  $O(n^(2\/3))$ in all.
Source: Lucy_Hedgehog's post in the Project Euler problem 10
  forum; min_25's blog.
Status: stress-tested against a sieve
"""

from math import isqrt


class PrimeCount:
  def __init__(self, n):
    r = isqrt(n)
    # lo[v] = pi(v), hi[i] = pi(n // i)
    lo = [max(i - 1, 0) for i in range(r + 1)]
    hi = [0] + [n // i - 1 for i in range(1, r + 1)]
    for p in range(2, r + 1):
      if lo[p] == lo[p - 1]:  # p is not prime
        continue
      cnt, q = lo[p - 1], p * p
      m = min(r, n // q)
      k = min(m, r // p)
      for i in range(1, k + 1):
        hi[i] -= hi[i * p] - cnt
      for i in range(k + 1, m + 1):
        hi[i] -= lo[n // (i * p)] - cnt
      for v in range(r, q - 1, -1):
        lo[v] -= lo[v // p] - cnt
    self.n, self.r, self.lo, self.hi = n, r, lo, hi

  def __call__(self, v):
    return self.lo[v] if v <= self.r else self.hi[self.n // v]
