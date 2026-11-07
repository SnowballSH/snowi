"""
Title: Lagrange interpolation
Description: `lagrange(y, x)` returns $f(x)$, where $f$ is the
  polynomial of degree below `len(y)` with $f(i) = y_i$ for every
  index $i$. `interpolate(xs, ys)` returns the coefficients of the
  polynomial of degree below $n$ through $n$ points $(x_i, y_i)$,
  lowest degree first.
Precondition: Arithmetic modulo the prime `MOD`; values in
  `[0, MOD)`. `lagrange`: `y` not empty, `len(y)` $<=$ `MOD`, any
  integer `x`. `interpolate`: the `xs` are distinct and in
  `[0, MOD)`.
Time: $O(n)$ for `lagrange`, 0.2 s for $n = 10^6$; $O(n^2)$ for
  `interpolate`, 2 s for $n = 5000$.
Usage:
  s = lagrange([0, 1, 5, 14], 10)  # sum_(i <= 10) i^2: 385
  c = interpolate(xs, ys)  # f(t) = sum_k c[k] t^k
Tweak: Nodes $a, a + h, a + 2h, dots$: evaluate at $(x - a) / h$
  (division mod `MOD`). An exact integer answer without a
  modulus: interpolate modulo two or three primes and combine by
  CRT; never in floating point.
Uses: numtheory/modpow
Example:
  >>> k, n = 3, 10**18  # sum_(i = 1)^n i^k: degree k + 1
  >>> y = [0] * (k + 2)
  >>> for i in range(1, k + 2):
  ...   y[i] = (y[i - 1] + pow(i, k, MOD)) % MOD
  >>> lagrange(y, n) == (n * (n + 1) // 2) ** 2 % MOD
  True
  >>> interpolate([1, 2, 3], [2, 5, 10])  # t^2 + 1
  [1, 0, 1]
Notes: The answer is a polynomial in a huge parameter $K$ when a
  DP over values $0, dots, K$ only adds, multiplies and takes
  prefix sums: a prefix sum raises the degree by one, a product
  adds degrees. Compute the answer for $K = 0, dots, d$ by brute
  force or the slow DP and call `lagrange`. Take two extra samples
  and check that they are reproduced: a wrong degree bound fails
  silently. With $min(K, a_i)$ or thresholds the answer is only
  piecewise polynomial; sample beyond the last breakpoint (or per
  piece), and above any small $K$ the statement treats specially.
  Cues: "$K$ up to $10^9$" next to $n <= 3000$, "values in
  $[1, K]$", "$sum i^k$", child $<=$ parent labellings of a tree
  (degree = size). A DP dimension that only adds under merges can
  be replaced by a variable $z$: run the DP at $n + 1$ numbers $z$
  and interpolate.
Source: Classical Lagrange formula with prefix and suffix
  products.
Status: stress-tested against Horner evaluation of random
  polynomials
"""


def lagrange(y, x):
  d, x = len(y) - 1, x % MOD
  suf, ifact = [1] * (d + 2), [1] * (d + 1)
  for i in range(d, -1, -1):
    suf[i] = suf[i + 1] * (x - i) % MOD
  for i in range(1, d + 1):
    ifact[d] = ifact[d] * i % MOD
  ifact[d] = pow(ifact[d], MOD - 2, MOD)
  for i in range(d, 0, -1):
    ifact[i - 1] = ifact[i] * i % MOD
  res, pre = 0, 1  # pre = prod_(j < i) (x - j)
  for i in range(d + 1):
    t = y[i] * pre % MOD * suf[i + 1] % MOD * ifact[i] % MOD
    t = t * ifact[d - i] % MOD
    res += -t if (d - i) & 1 else t
    pre = pre * (x - i) % MOD
  return res % MOD


def interpolate(xs, ys):
  n = len(xs)
  P, Q, res = [1] + [0] * n, [0] * n, [0] * n
  for j, xj in enumerate(xs):  # P = prod_j (z - x_j)
    for k in range(j + 1, 0, -1):
      P[k] = (P[k - 1] - xj * P[k]) % MOD
    P[0] = -xj * P[0] % MOD
  for xi, yi in zip(xs, ys):
    c = den = 0  # Q = P / (z - x_i), den = Q(x_i)
    for k in range(n, 0, -1):
      c = (P[k] + c * xi) % MOD
      Q[k - 1] = c
      den = (den * xi + c) % MOD
    c = yi * pow(den, MOD - 2, MOD) % MOD
    for k in range(n):
      res[k] = (res[k] + c * Q[k]) % MOD
  return res
