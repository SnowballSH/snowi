"""
Title: Modular arithmetic
Description: `MOD` and the Fermat inverse that other snippets use;
  `pow(a, -1, m)` inverts modulo any `m` coprime to `a`.
Precondition: `MOD` prime for the Fermat inverse `modinv`.
Time: $O(log e)$.
Usage:
  x = pow(3, 10**12, MOD); y = modinv(x)
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
