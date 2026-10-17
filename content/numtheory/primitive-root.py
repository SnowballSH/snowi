"""
Title: Primitive root
Description: The smallest primitive root `g` of a prime `p`.
Precondition: `p` prime.
Time: One `factor(p - 1)` plus $O(g log p)$ per prime factor; the
  smallest root is at most 113 for $p < 10^8$. PyPy: about 2 ms for
  a random $p < 2^63$.
Usage:
  g = primitive_root(998244353)  # 3
Uses: numtheory/factor
Example:
  >>> primitive_root(10**9 + 7), primitive_root(2**61 - 1)
  (5, 37)
  >>> p, a = 13, 3  # order of a mod p: smallest e >= 1, a^e = 1
  >>> e = p - 1
  >>> for q in factor(p - 1):
  ...   if pow(a, e // q, p) == 1:
  ...     e //= q
  >>> e
  3
  >>> p, k = 998244353, 1 << 20  # k-th root of unity, k | p - 1
  >>> w = pow(primitive_root(p), (p - 1) // k, p)
  >>> pow(w, k, p), pow(w, k // 2, p) == p - 1
  (1, True)
Notes: Primitive roots modulo `n` exist iff `n` is 1, 2, 4,
  $p^k$ or $2 p^k$ for an odd prime `p`. If `g` is one mod `p`,
  then `g` or $g + p$ is one mod every $p^k$ ($g + p$ exactly
  when $g^(p - 1) equiv 1 thick (mod p^2)$); mod $2 p^k$ take the
  odd one of `g` and $g + p^k$. There are $phi(p - 1)$ primitive
  roots, and $g^k$ is one iff $gcd(k, p - 1) = 1$. The order of
  any `a` divides $p - 1$, and $x^k equiv 1$ has $gcd(k, p - 1)$
  solutions. Cues: NTT roots of unity, discrete logs that turn
  products into sums (a product over a multiset becomes a
  convolution of exponents, UCup S4-21 J), counting `k`-th power
  residues.
Source: Classic (Gauss); cp-algorithms "Primitive root".
Status: stress-tested against brute-force orders
"""


def primitive_root(p):
  if p == 2:
    return 1
  qs = set(factor(p - 1))
  g = 2
  while any(pow(g, (p - 1) // q, p) == 1 for q in qs):
    g += 1
  return g
