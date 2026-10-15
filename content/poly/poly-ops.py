"""
Title: Power series operations
Description: Formal power series modulo $x^n$ and `MOD`:
  `poly_inv` gives $1 / a$, `poly_log` gives $ln a$, `poly_exp`
  gives $exp a$ and `poly_pow` gives $a^k$, each as `n`
  coefficients. `poly_divmod(a, b)` divides polynomials:
  $a = q b + r$ with `len(q) = len(a) - len(b) + 1` (empty when
  that is not positive) and `len(r) = len(b) - 1`.
Precondition: Entries in `[0, MOD)`, `n` $>= 0$. `poly_inv` needs
  $a_0 != 0$, `poly_log` needs $a_0 = 1$, `poly_exp` needs
  $a_0 = 0$; `poly_pow` takes any `a` and $k >= 0$. `poly_divmod`
  needs `b[-1]` $!= 0$. Lengths up to $2^22$ (NTT size).
Time: $O(n log n)$. For $n = 5 dot 10^5$: inverse 0.55 s, log 0.95
  s, exp 1.2 s, pow 2.2 s (`conv` of that length: 0.35 s). Newton
  steps and NTT sizes round up to powers of two: $n = 2^18$ takes
  half as long, $n = 2^18 + 1$ as long as $5 dot 10^5$. Below
  $n approx 3000$ the $O(n^2)$ recurrences of the generating
  functions note are simpler.
Usage:
  b = poly_inv(a, n)  # a b = 1 mod x^n
  e = poly_exp(poly_log(a, n), n)  # e = a mod x^n
  q, r = poly_divmod(a, b)
Tweak: Square root when $a_0 = 1$:
  `poly_pow(a, (MOD + 1) // 2, n)` (an exponent in the field works
  because $a_0 = 1$). Otherwise strip $x^v$ and the leading
  coefficient $c$ as `poly_pow` does, take that root and multiply
  by $x^(v \/ 2)$ and `mod_sqrt(c, MOD)` (`numtheory/mod-sqrt`);
  there is no root if $v$ is odd or $c$ is not a square. Another
  NTT prime: change `MOD` (and the root 3 in `ntt_roots` if
  needed).
Uses: poly/ntt
Example:
  >>> poly_inv([1, MOD - 1], 5)  # 1 / (1 - x)
  [1, 1, 1, 1, 1]
  >>> e = poly_exp([0, 1], 5)  # 1 / k!
  >>> [x * f % MOD for x, f in zip(e, [1, 1, 2, 6, 24])]
  [1, 1, 1, 1, 1]
  >>> poly_pow([0, 0, 1, 1], 3, 8)  # x^6 (1 + x)^3
  [0, 0, 0, 0, 0, 0, 1, 3]
  >>> h = (MOD + 1) // 2
  >>> s = poly_pow([1, MOD - 4], h, 7)  # sqrt(1 - 4x)
  >>> [-x * modinv(2) % MOD for x in s[1:]]  # Catalan numbers
  [1, 1, 2, 5, 14, 42]
  >>> poly_divmod([5, 0, 3, 1], [1, 1])  # divide by x + 1
  ([998244351, 2, 1], [7])
  >>> n = 8  # partitions: exp of sum_k sum_j x^(jk) / j
  >>> inv = [1] * n  # inv[j] = 1 / j
  >>> for i in range(2, n):
  ...   inv[i] = (MOD - MOD // i) * inv[MOD % i] % MOD
  >>> L = [0] * n
  >>> for k in range(1, n):
  ...   for j in range(1, (n - 1) // k + 1):
  ...     L[j * k] += inv[j]
  >>> poly_exp([x % MOD for x in L], n)
  [1, 1, 2, 3, 5, 7, 11, 15]
Notes: Newton doubles the number of correct terms per step:
  $b <- b (2 - a b)$ for the inverse; exp keeps $g = 1 / f$
  alongside and sets $f <- f (1 + a - ln f)$. The products are
  cyclic of length $2m$ (`cyc_mul`): the wrapped terms only land
  where the result is already known. $ln a = integral a' / a$.
  `poly_pow` strips the leading $x^v$ and coefficient $c$, then
  takes $exp(k ln a)$; the result is 0 when $v k >= n$. Reduce $k$
  mod `MOD` only inside the log, never in $c^k$ or in $v k$. Cues:
  labelled sets of connected pieces (exp), connected from all
  (log), sequences $1 / (1 - A)$ (inverse), $A^k$ for huge $k$,
  "remainder of a polynomial division" (divmod). Build inputs in
  $O(n log n)$ too: a `modinv` per term of a harmonic double loop
  costs more than the `poly_exp`; use the linear inverse table of
  the Example.
Source: cp-algorithms, Operations on polynomials and series; exp
  with the inverse kept alongside after Bernstein, Removing
  redundancy in high-precision Newton iteration.
Status: stress-tested against $O(n^2)$ series arithmetic
"""


def cyc_mul(a, B):  # a b mod (x^len(a) - 1) for B = ntt(b)
  ntt(a)  # in place: a is overwritten
  a = [x * y % MOD for x, y in zip(a, B)]
  intt(a)
  return a


def poly_inv_step(a, b):  # 1 / a: mod x^k -> mod x^(2k)
  k = len(b)
  f, g = a[: 2 * k], b + [0] * k
  f += [0] * (2 * k - len(f))
  ntt(g)
  f = cyc_mul(f, g)
  f[:k] = [0] * k  # a b - 1
  f = cyc_mul(f, g)
  return b + [-x % MOD for x in f[k:]]


def poly_inv(a, n):
  b = [pow(a[0], MOD - 2, MOD)]
  while len(b) < n:
    b = poly_inv_step(a, b)
  return b[:n]


def poly_log(a, n):
  if n <= 1:
    return [0] * n
  d, inv = [0] * (n - 1), [1] * n
  for i in range(1, min(len(a), n)):
    d[i - 1] = a[i] * i % MOD
  for i in range(2, n):
    inv[i] = (MOD - MOD // i) * inv[MOD % i] % MOD
  d = conv(d, poly_inv(a, n - 1))
  return [0] + [d[i - 1] * inv[i] % MOD for i in range(1, n)]


def poly_exp(a, n):
  inv = [1] * (2 * n + 1)
  for i in range(2, 2 * n + 1):
    inv[i] = (MOD - MOD // i) * inv[MOD % i] % MOD
  da = [i * x % MOD for i, x in enumerate(a)][1:]  # a'
  f, g, m = [1], [1], 1  # f = exp a, g = 1 / f
  while m < n:
    if len(g) < m:
      g = poly_inv_step(f, g)
    F, G, P = f + [0] * m, g + [0] * m, da[: 2 * m]
    P += [0] * (2 * m - len(P))
    ntt(F)
    ntt(G)
    P = cyc_mul(P, F)
    P[: m - 1] = [0] * (m - 1)  # f a' - f'
    P = cyc_mul(P, G)  # a' - (ln f)'
    T = [P[i - 1] * inv[i] % MOD for i in range(m, 2 * m)]
    f += cyc_mul([0] * m + T, F)[m:]  # f + f (a - ln f)
    m *= 2
  return f[:n]


def poly_pow(a, k, n):
  if k == 0:
    return [1] * min(n, 1) + [0] * (n - 1)
  v = 0
  while v < len(a) and a[v] == 0:
    v += 1
  if v == len(a) or v * k >= n:
    return [0] * n
  m, c = n - v * k, pow(a[v], MOD - 2, MOD)
  b = poly_log([x * c % MOD for x in a[v : v + m]], m)
  b = poly_exp([x * (k % MOD) % MOD for x in b], m)
  ck = pow(a[v], k, MOD)
  return [0] * (n - m) + [x * ck % MOD for x in b]


def poly_divmod(a, b):
  n, m = len(a), len(b)
  q, k = [], n - m + 1
  if k > 0:
    q = conv(a[::-1][:k], poly_inv(b[::-1], k))[:k][::-1]
  r = a[: m - 1] + [0] * (m - 1 - n)
  for i, x in enumerate(conv(b, q)[: m - 1]):
    r[i] = (r[i] - x) % MOD
  return q, r
