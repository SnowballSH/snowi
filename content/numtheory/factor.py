"""
Title: Primality and factoring
Description: `is_prime(n)`: deterministic Miller-Rabin for all
  $n < 2^64$. `factor(n)`: the prime factors of `n` in increasing
  order, with multiplicity, by Pollard-Brent rho.
Precondition: $n >= 1$; `factor(1)` is empty. Above $2^64$ the
  fixed bases are not proven: `is_prime` may pass a composite, so add
  random bases if inputs may be adversarial.
Time: `is_prime`: 7 modular powers. `factor`: $O(n^(1/4))$
  expected. PyPy: 0.3 ms on average for a random $n <= 10^18$. The
  worst case, two primes near $sqrt(n)$, takes about 10 ms below
  $2^63$ and 35 ms up to $2^64$, where every value is a big int.
Usage:
  p = is_prime(10**9 + 7)
  f = factor(360)  # [2, 2, 2, 3, 3, 5]
Tweak: Group into (prime, exponent) pairs with
  `sorted(Counter(factor(n)).items())`.
Example:
  >>> is_prime(2**61 - 1), is_prime(3825123056546413051)
  (True, False)
  >>> factor(10**18 - 1)
  [3, 3, 3, 3, 7, 11, 13, 19, 37, 52579, 333667]
  >>> from collections import Counter
  >>> divs = [1]  # all divisors of n, unsorted
  >>> for p, e in sorted(Counter(factor(360)).items()):
  ...   divs = [d * p**k for d in divs for k in range(e + 1)]
  >>> len(divs), sum(divs)
  (24, 1170)
Notes: Factor once, then derive what you need: the divisors (at
  most 103680 for $n <= 10^18$), $phi(n) = n product_p (1 - 1/p)$,
  or the multiplicative order of `a` modulo a prime `p`: start from
  $o = p - 1$ and divide out each prime factor `q` while
  $a^(o/q) = 1$. Rho splits off the smallest prime factor `r` in
  about $sqrt(r)$ steps, so the worst case is the largest balanced
  semiprime the input allows: budget for it, not for a random `n`.
Source: Miller-Rabin bases by Jim Sinclair; Brent (1980); kactl
  MillerRabin and Factor (CC0).
Status: stress-tested against trial division and random-base
  Miller-Rabin
"""

from math import gcd


def is_prime(n):
  if n < 2 or n % 2 == 0:
    return n == 2
  d, s = n - 1, 0
  while d % 2 == 0:
    d, s = d // 2, s + 1
  for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
    x = pow(a, d, n)
    if a % n == 0 or x == 1 or x == n - 1:
      continue
    for _ in range(s - 1):
      x = x * x % n
      if x == n - 1:
        break
    else:
      return False
  return True


def rho(n):  # a proper divisor of an odd composite n
  c = 1
  while True:
    x = y = ys = 2
    g = q = r = 1
    while g == 1:
      x = y
      for _ in range(r):
        y = (y * y + c) % n
      k = 0
      while k < r and g == 1:
        ys = y
        for _ in range(min(128, r - k)):
          y = (y * y + c) % n
          q = q * abs(x - y) % n
        g = gcd(q, n)
        k += 128
      r *= 2
    if g == n:  # redo the last batch step by step
      g = 1
      while g == 1:
        ys = (ys * ys + c) % n
        g = gcd(abs(x - ys), n)
    if g != n:
      return g
    c += 1


def factor(n):
  res, todo = [], [n]
  while todo:
    m = todo.pop()
    if m == 1:
      continue
    if is_prime(m):
      res.append(m)
      continue
    d = rho(m) if m % 2 else 2
    todo += [d, m // d]
  return sorted(res)
