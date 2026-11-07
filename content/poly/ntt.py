"""
Title: Number theoretic transform
Description: `conv(a, b)` multiplies polynomials modulo `MOD`.
  `ntt(a)` replaces `a` by its values at the powers of a primitive
  `len(a)`-th root of unity, in bit-reversed order, and `intt(a)`
  undoes it.
Precondition: `MOD` $= 998244353 = 119 dot 2^23 + 1$ with primitive
  root 3. Entries in `[0, MOD)`. `len(a)` is a power of two, at most
  $2^23$, so `conv` needs `len(a) + len(b) - 1` $<= 2^23$.
Time: $O(n log n)$. A $2^18 times 2^18$ product takes about 0.6 s.
Usage:
  c = conv(a, b)  # len(a) + len(b) - 1 entries
Tweak: Other primes: the table in Polynomial techniques.
Uses: numtheory/modpow
Example:
  >>> conv([1, 2], [3, 4])
  [3, 10, 8]
  >>> a, b = [1, 2, 0, 0], [5, 1, 1, 1]
  >>> ntt(a); ntt(b)
  >>> c = [x * y % MOD for x, y in zip(a, b)]
  >>> intt(c); c  # cyclic: 2 x^4 wraps to x^0
  [7, 11, 3, 3]
  >>> a, b = [1, 2, 3], [4, 5, 6, 7, 8]
  >>> conv(a[::-1], b)[2:5]  # sum_i a_i b_(i + d), d = 0, 1, 2
  [32, 38, 44]
Notes: Radix 4: forward by decimation in frequency, inverse in time,
  so neither needs a bit-reversal pass, and one flat loop per level,
  which PyPy compiles better than nested loops. `conv` falls back to
  the schoolbook product up to 40 entries, which matters for
  divide-and-conquer products of many small polynomials. Cues: "for
  every $k$" with a sum over $i + j = k$, a product of many small
  polynomials, the modulus 998244353 itself. For a correlation
  $sum_i a_i b_(i + d)$ reverse `a`; the answer for `d` is
  `c[len(a) - 1 + d]`. Reduce negative values first.
Source: Radix-4 layout after the AtCoder Library (CC0).
Status: stress-tested against the schoolbook product
"""

NTT_RT = [0, 1]  # NTT_RT[k + j] = w^j, w a primitive 2k-th root


def ntt_roots(n):
  rt = NTT_RT
  while len(rt) < max(n, 4):
    k = len(rt)
    z = pow(3, (MOD - 1) // (2 * k), MOD)
    for i in range(k, 2 * k):
      rt.append(rt[i >> 1] * z % MOD if i & 1 else rt[i >> 1])
  return rt


def ntt(a):
  n, rt = len(a), ntt_roots(len(a))
  m, im = n >> 2, rt[3]  # im^2 = -1
  while m:
    m2, m3 = 2 * m, 3 * m
    for t in range(n >> 2):
      j = t & (m - 1)
      p, w2, w1 = 4 * t - 3 * j, rt[m + j], rt[m2 + j]
      w3 = w1 * w2 % MOD
      x0, x1, x2, x3 = a[p], a[p + m], a[p + m2], a[p + m3]
      s, d, u = x0 + x2, x0 - x2, x1 + x3
      e = (x1 - x3) * im % MOD
      a[p], a[p + m] = (s + u) % MOD, (s - u) * w2 % MOD
      a[p + m2] = (d + e) * w1 % MOD
      a[p + m3] = (d - e) * w3 % MOD
    m >>= 2
  if n.bit_length() % 2 == 0:
    for i in range(0, n, 2):
      x, y = a[i], a[i + 1]
      a[i], a[i + 1] = (x + y) % MOD, (x - y) % MOD


def intt(a):
  n, rt = len(a), ntt_roots(len(a))
  m, im = 1, rt[3]
  if n.bit_length() % 2 == 0:
    for i in range(0, n, 2):
      a[i], a[i + 1] = a[i] + a[i + 1], a[i] - a[i + 1]
    m = 2
  while m < n:
    m2, m3 = 2 * m, 3 * m
    for t in range(n >> 2):
      j = t & (m - 1)
      p, w2, w1 = 4 * t - 3 * j, rt[m + j], rt[m2 + j]
      w3 = w1 * w2 % MOD
      x0, x1 = a[p] % MOD, a[p + m] * w2 % MOD
      x2, x3 = a[p + m2] * w1 % MOD, a[p + m3] * w3 % MOD
      s, d, u = x0 + x1, x0 - x1, x2 + x3
      e = (x2 - x3) * im % MOD
      a[p], a[p + m2] = s + u, s - u  # |a[p]| < 4 MOD
      a[p + m], a[p + m3] = d + e, d - e
    m *= 4
  a[1:] = a[:0:-1]  # forward roots: undo k -> -k
  inv = pow(n, MOD - 2, MOD)
  for i in range(n):
    a[i] = a[i] * inv % MOD


def conv(a, b):
  if not a or not b:
    return []
  s = len(a) + len(b) - 1
  if min(len(a), len(b)) <= 40:
    c = [0] * s
    for i, x in enumerate(a):
      for j, y in enumerate(b):
        c[i + j] = (c[i + j] + x * y) % MOD
    return c
  n = 1 << (s - 1).bit_length()
  a, b = a + [0] * (n - len(a)), b + [0] * (n - len(b))
  ntt(a)
  ntt(b)
  for i in range(n):
    a[i] = a[i] * b[i] % MOD
  intt(a)
  return a[:s]
