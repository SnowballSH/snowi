"""
Title: Convolution modulo any m
Description: `conv_mod(a, b, m)` multiplies polynomials modulo any
  `m`, prime or not: `c[k]` is $sum_(i + j = k) a_i b_j mod m$. It
  splits each value into two 15-bit halves and needs four FFTs.
Precondition: $2 <= m < 2^31$, entries in `[0, m)`. Exact for
  result lengths up to $2^20$ when $m < 2^30$ (e.g. $10^9 + 7$),
  up to $2^19$ when $m < 2^31$; random inputs allow more.
Time: $O(n log n)$. A $2^18 times 2^18$ product takes about 0.4 s.
Usage:
  c = conv_mod(a, b, 10**9 + 7)
Tweak: For $m = 998244353$ use the faster `conv`. Three NTT primes
  with CRT are an exact alternative while $n m^2 < 2^86$, but two to
  three times slower in PyPy.
Uses: poly/fft
Example:
  >>> conv_mod([10**9, 5], [10**9, 3], 10**9 + 7)
  [49, 999999951, 15]
  >>> c = conv_mod([10**9] * 50, [10**9] * 50, 10**9 + 7)
  >>> c[0], c[49], len(c)  # 10^9 = -7: each product is 49
  (49, 2450, 99)
Notes: Values are first shifted into $(-m/2, m/2]$ and split into
  balanced digits in $[-2^14, 2^14)$. Their products then cancel in
  the sums, which makes the floats hundreds of times more accurate on
  random inputs than unsigned halves; values near 0 or $m$ become
  small. The worst case is all values near $m/2$. The two packed
  transforms are in bit-reversed order, where frequency $-k$ of
  position `p` sits at `p ^ (t - 1)`, `t` the highest power of two
  $<= p$. Up to 40 entries it multiplies directly. Cues: an answer
  modulo $10^9 + 7$ that needs a convolution, or a small modulus
  like 1009 in a divide-and-conquer product.
Source: kactl FastFourierTransformMod (CC0); balanced digits are a
  standard precision trick.
Status: stress-tested against the schoolbook product
"""


def conv_mod(a, b, m):
  if not a or not b:
    return []
  s = len(a) + len(b) - 1
  if min(len(a), len(b)) <= 40:
    c = [0] * s
    for i, x in enumerate(a):
      for j, y in enumerate(b):
        c[i + j] = (c[i + j] + x * y) % m
    return c
  n = 1 << (s - 1).bit_length()
  lr, li, rr, ri = [0.0] * n, [0.0] * n, [0.0] * n, [0.0] * n
  for xs, high, low in ((a, lr, li), (b, rr, ri)):
    for i, x in enumerate(xs):
      x -= m if 2 * x > m else 0
      y = ((x + 16384) & 32767) - 16384
      high[i], low[i] = float((x - y) >> 15), float(y)
  fft(lr, li)
  fft(rr, ri)
  ur, ui, vr, vi = [0.0] * n, [0.0] * n, [0.0] * n, [0.0] * n
  t = 1
  for p in range(n):
    if p == 2 * t:
      t = p
    q = p ^ (t - 1)  # frequency -k
    ar, ai = lr[p] + lr[q], li[p] - li[q]  # 2 * FFT(high)
    br, bi = li[p] + li[q], lr[q] - lr[p]  # 2 * FFT(low)
    xr, xi = rr[p] * 0.5, ri[p] * 0.5
    ur[p], ui[p] = ar * xr - ai * xi, ar * xi + ai * xr
    vr[p], vi[p] = br * xr - bi * xi, br * xi + bi * xr
  ifft(ur, ui)  # high * high + i high * low
  ifft(vr, vi)  # low * high + i low * low
  c = [0] * s
  for i in range(s):
    x, y = round(ur[i]), round(ui[i]) + round(vr[i])
    c[i] = ((((x % m << 15) + y) % m << 15) + round(vi[i])) % m
  return c
