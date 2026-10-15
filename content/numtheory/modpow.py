"""
Title: Modular arithmetic
Description: Python has modular power built in: `pow(b, e, m)`, and
  `pow(a, -1, m)` for the inverse when `gcd(a, m) = 1`. `MOD` is the
  modulus other snippets assume.
Precondition: `MOD` prime for the Fermat inverse `modinv`.
Time: $O(log e)$.
Usage:
  x = pow(3, 10**12, MOD); y = modinv(x)
Tweak: Change `MOD`. Keep values reduced: big ints are slow in PyPy.
Example:
  >>> modinv(2) * 2 % MOD
  1
  >>> pow(5, -1, 7)
  3
Status: stress-tested against repeated multiplication
"""

MOD = 998244353


def modinv(a):
  return pow(a, MOD - 2, MOD)
