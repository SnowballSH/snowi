"""
Title: Modular square root
Description: A square root of `a` modulo a prime `p` by
  Tonelli-Shanks: returns the smaller root $r <= p - r$ with
  $r^2 equiv a thick (mod p)$, or -1 if `a` is not a square.
Precondition: `p` prime; `a` any integer.
Time: $O(log^2 p)$ worst case, when $2^k | p - 1$ for large $k$.
  PyPy: about 4 µs for $p = 998244353$.
Usage:
  r = mod_sqrt(a, p)  # -1, or roots r and p - r
Tweak: For $p equiv 3 thick (mod 4)$ the root is simply
  $a^((p + 1)\/4)$ (check its square).
Example:
  >>> mod_sqrt(2, 7), mod_sqrt(3, 7), mod_sqrt(-1, 13)
  (3, -1, 5)
  >>> p, b, c = 10**9 + 7, 5, 6  # x^2 + b x + c = 0 mod p
  >>> r = mod_sqrt(b * b - 4 * c, p)  # -1 means no root
  >>> h = (p + 1) // 2  # inverse of 2
  >>> sorted((-b + s) * h % p for s in (r, -r))
  [1000000004, 1000000005]
Notes: The roots are `r` and $p - r$, distinct unless $r = 0$ or
  $p = 2$.
  Legendre symbol: $a^((p - 1)\/2) mod p$ is 1 for a nonzero
  square, $p - 1$ for a non-square, 0 for $a equiv 0$; half of the
  nonzero residues are squares. Solve $x^2 + b x + c equiv 0$ for
  odd `p` as $x = (-b plus.minus r)\/2$ with $r^2 equiv b^2 - 4c$.
  Modulo $p^k$ ($p$ odd, $p divides.not a$) lift a root with
  $r <- r - (r^2 - a)\/(2r)$, which doubles the correct digits;
  for a composite modulus combine prime powers with CRT. Cues: the
  constant term of a power series square root, "is $x^2 equiv a$
  solvable", $sqrt(5)$ for Binet's Fibonacci formula (for odd `p`
  it exists iff $p equiv plus.minus 1 thick (mod 5)$ or $p = 5$).
Source: Tonelli (1891), Shanks (1973).
Status: stress-tested against all squares modulo small primes
"""


def mod_sqrt(a, p):
  a %= p
  if a < 2:  # also p = 2
    return a
  if pow(a, (p - 1) // 2, p) != 1:
    return -1
  s, q, z = 0, p - 1, 2
  while q % 2 == 0:
    s, q = s + 1, q // 2
  while pow(z, (p - 1) // 2, p) == 1:
    z += 1
  c, x, t = pow(z, q, p), pow(a, (q + 1) // 2, p), pow(a, q, p)
  while t != 1:  # x^2 = a t, t has order 2^i < 2^s
    i, u = 0, t
    while u != 1:
      i, u = i + 1, u * u % p
    b = pow(c, 1 << (s - i - 1), p)
    x, c = x * b % p, b * b % p
    t, s = t * c % p, i
  return min(x, p - x)
