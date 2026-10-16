"""
Title: Extended Euclid and CRT
Description: `euclid(a, b)` returns `(g, x, y)` with $g = gcd(a, b)$
  and $a x + b y = g$. `inv_mod(a, m)` is the inverse of `a` modulo
  any `m`, or -1 if $gcd(a, m) != 1$. `crt(a, m, b, n)` solves
  $x equiv a thick (mod m)$, $x equiv b thick (mod n)$ for any
  moduli, coprime or not: it returns $(x, op("lcm")(m, n))$ with
  $0 <= x < op("lcm")$, or $(-1, -1)$ if there is no solution.
  `crt_all` folds lists of residues and moduli the same way.
Precondition: `euclid`: $a, b >= 0$. Moduli are positive; residues
  may be any integers.
Time: $O(log min(a, b))$ per call.
Usage:
  g, x, y = euclid(a, b)
  inv = inv_mod(3, 10)  # 7
  r, l = crt(2, 4, 4, 6)  # (10, 12)
  r, l = crt_all([1, 2, 3], [2, 3, 5])  # (23, 30)
Example:
  >>> euclid(240, 46)
  (2, -9, 47)
  >>> inv_mod(4, 10), inv_mod(-3, 10)
  (-1, 3)
  >>> crt(2, 4, 3, 6), crt(-1, 10**18, 1, 3)
  ((-1, -1), (1999999999999999999, 3000000000000000000))
  >>> a, b, c = 6, -10, 8  # a x + b y = c, a or b nonzero
  >>> g, x, y = euclid(abs(a), abs(b))
  >>> k = c // g  # needs c % g == 0, else no solution
  >>> x *= k * (-1 if a < 0 else 1)
  >>> y *= k * (-1 if b < 0 else 1)
  >>> (x, y), a * x + b * y == c
  ((8, 4), True)
Notes: A linear congruence $a x equiv c thick (mod m)$ is solvable
  iff $g = gcd(a, m)$ divides `c`, with `g` solutions modulo `m`:
  divide everything by `g` and multiply by `inv_mod(a/g, m/g)`. All
  solutions of $a x + b y = c$ are $(x + k b/g, y - k a/g)$. Use CRT
  when the modulus is composite (solve modulo each prime power, then
  combine) or when a quantity is known modulo several periods.
Source: kactl Euclid and CRT (CC0); AtCoder Library inv_mod/crt
  (CC0).
Status: stress-tested against brute-force search
"""


def euclid(a, b):
  x, y, u, v = 1, 0, 0, 1
  while b:
    q = a // b
    a, b, x, u, y, v = b, a - q * b, u, x - q * u, v, y - q * v
  return a, x, y


def inv_mod(a, m):
  g, x, _ = euclid(a % m, m)
  return x % m if g == 1 else -1


def crt(a, m, b, n):
  g, x, _ = euclid(m, n)
  if (b - a) % g:
    return -1, -1
  l = m // g * n
  return (a + (b - a) // g * x % (n // g) * m) % l, l


def crt_all(rs, ms):
  x, l = 0, 1
  for r, m in zip(rs, ms):
    x, l = crt(x, l, r, m)
    if l < 0:
      break
  return x, l
