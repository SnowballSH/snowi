"""
Title: Fast Fourier transform
Description: `fft_conv(a, b)` multiplies polynomials with real
  coefficients in floats. `fft(re, im)` evaluates the complex
  vector `re + i im` in place at the powers of
  $omega = e^((2 pi i) / n)$, in bit-reversed order, and
  `ifft(re, im)` undoes it.
Precondition: `len(re)` is a power of two. Results round exactly
  while $norm(a) norm(b) < 4 dot 10^14$ (Euclidean norms), which
  holds if $n dot max abs(a_i) dot max abs(b_j) < 4 dot 10^14$ for
  the longer length $n$, e.g. $2^18$ values up to $3 dot 10^4$.
  Small results do not help: with mixed signs the error still
  follows the inputs. Beyond, use `conv_mod`.
Time: $O(n log n)$. A $2^18 times 2^18$ product takes about 0.2 s.
Usage:
  c = [round(x) for x in fft_conv(a, b)]
Example:
  >>> [round(x) for x in fft_conv([1, 2, 3], [4, 5])]
  [4, 13, 22, 15]
  >>> re, im = [1.0, 2.0, 0.0, 0.0], [0.0] * 4
  >>> fft(re, im)  # 1 + 2x at 1, -1, i, -i (bit-reversed order)
  >>> [round(x) for x in re], [round(x) for x in im]
  ([3, -1, 1, 1], [0, 0, 2, -2])
  >>> ifft(re, im); [round(x) for x in re]
  [1, 2, 0, 0]
Notes: Real and imaginary parts live in separate float lists, which
  PyPy stores unboxed; a list of `complex` is several times slower.
  Forward is decimation in frequency and inverse in time, so neither
  needs a bit-reversal pass, and each level is one flat loop, which
  PyPy runs faster than nested loops. The bound is empirical: the
  worst case found (unbalanced magnitudes, a result just past a
  power of two) rounds within 0.41 at $4 dot 10^14$ and fails near
  $4.7 dot 10^14$. Cues: integer convolutions whose exact values
  exceed the NTT modulus (counts of pair sums, wildcard string
  matching via $sum (x - y)^2 x y$), real-valued signals.
Source: Classic radix-2 Cooley-Tukey; layout as in kactl (CC0).
Status: stress-tested against the schoolbook product
"""

from math import cos, pi, sin

# FFT_RE[k + j] + i FFT_IM[k + j] = e^(i pi j / k)
FFT_RE, FFT_IM = [0.0, 1.0], [0.0, 0.0]


def fft_roots(n):
  rr, ri = FFT_RE, FFT_IM
  while len(rr) < n:
    k = len(rr)
    for j in range(k):
      rr.append(cos(pi * j / k))
      ri.append(sin(pi * j / k))
  return rr, ri


def fft(re, im):
  n = len(re)
  rr, ri = fft_roots(n)
  k = n >> 1
  while k:
    for t in range(n >> 1):
      j = t & (k - 1)
      p = 2 * t - j
      q, wr, wi = p + k, rr[k + j], ri[k + j]
      xr, xi = re[p] - re[q], im[p] - im[q]
      re[p] += re[q]
      im[p] += im[q]
      re[q], im[q] = xr * wr - xi * wi, xr * wi + xi * wr
    k >>= 1


def ifft(re, im):
  n = len(re)
  rr, ri = fft_roots(n)
  k = 1
  while k < n:
    for t in range(n >> 1):
      j = t & (k - 1)
      p = 2 * t - j
      q, wr, wi = p + k, rr[k + j], ri[k + j]
      xr = re[q] * wr + im[q] * wi  # times conj(w)
      xi = im[q] * wr - re[q] * wi
      re[q], im[q] = re[p] - xr, im[p] - xi
      re[p] += xr
      im[p] += xi
    k *= 2
  for i in range(n):
    re[i] /= n
    im[i] /= n


def fft_conv(a, b):
  if not a or not b:
    return []
  s = len(a) + len(b) - 1
  n = 1 << (s - 1).bit_length()
  ar = [float(x) for x in a] + [0.0] * (n - len(a))
  br = [float(x) for x in b] + [0.0] * (n - len(b))
  ai, bi = [0.0] * n, [0.0] * n
  fft(ar, ai)
  fft(br, bi)
  for i in range(n):
    x, y, u, v = ar[i], ai[i], br[i], bi[i]
    ar[i], ai[i] = x * u - y * v, x * v + y * u
  ifft(ar, ai)
  return ar[:s]
