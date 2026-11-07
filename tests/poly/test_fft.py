import cmath
import math

from testlib import Timer, check, check_eq, iters, load, rnd, rnd_list, rng

lib = load("poly/fft")


def brute_conv(a, b):
  if not a or not b:
    return []
  c = [0] * (len(a) + len(b) - 1)
  for i, x in enumerate(a):
    for j, y in enumerate(b):
      c[i + j] += x * y
  return c


def bitrev(k, bits):
  return int(format(k, "b").zfill(bits)[::-1], 2) if bits else 0


def check_constant(n, m, v, w):
  c = lib.fft_conv([v] * n, [w] * m)
  want = [min(k + 1, n, m, n + m - 1 - k) * v * w for k in range(n + m - 1)]
  check_eq([round(x) for x in c], want, f"constant {n}x{m}")


def check_alternating(n, m, v, w):
  """a_i = (-1)^i v, b = m copies of w: tiny results, large errors."""
  a = [-v if i % 2 else v for i in range(n)]
  c = lib.fft_conv(a, [w] * m)
  want = []
  for k in range(n + m - 1):
    lo, hi = max(0, k - m + 1), min(k, n - 1)
    want.append(0 if (hi - lo) % 2 else -v * w if lo % 2 else v * w)
  check_eq([round(x) for x in c], want, f"alternating {n}x{m}")


for it in range(iters(500)):
  n = rnd(0, 20 if it % 4 else 100)
  m = rnd(0, 20 if it % 4 else 100)
  # Magnitudes may differ a lot, e.g. weights times a 0/1 mask.
  ha = hb = 1000 if it % 3 else 10**6
  if it % 4 == 1:
    ha, hb = 10**9, 1
  if it % 4 == 3:
    ha, hb = 1, 10**9
  a, b = rnd_list(n, -ha, ha), rnd_list(m, -hb, hb)
  check_eq([round(x) for x in lib.fft_conv(a, b)], brute_conv(a, b))
  x = [rng.uniform(-1, 1) for _ in range(n)]
  y = [rng.uniform(-1, 1) for _ in range(m)]
  got, want = lib.fft_conv(x, y), brute_conv(x, y)
  check_eq(len(got), len(want))
  check(all(abs(u - v) < 1e-9 for u, v in zip(got, want)))

for lg in range(8):
  n = 1 << lg
  re = [rng.uniform(-1, 1) for _ in range(n)]
  im = [rng.uniform(-1, 1) for _ in range(n)]
  fr, fi = re[:], im[:]
  lib.fft(fr, fi)
  for k in range(n):
    w = cmath.exp(2j * math.pi * k / n)
    val = sum(complex(re[j], im[j]) * w**j for j in range(n))
    p = bitrev(k, lg)
    check(abs(complex(fr[p], fi[p]) - val) < 1e-9, f"dft n={n} k={k}")
  lib.ifft(fr, fi)
  check(all(abs(u - v) < 1e-12 for u, v in zip(fr + fi, re + im)))

# Precision bound from the header: |a| |b| < 4e14, e.g. 2^18 values up
# to 3e4; 2^19 values up to 2.7e4 give 3.8e14.
with Timer("fft_conv 2^19 x 2^19, worst case"):
  check_constant(1 << 19, 1 << 19, 27000, 27000)
check_alternating(1 << 19, 1 << 19, 27000, 27000)
check_constant(1 << 18, 1 << 18, 30000, 30000)
check_constant(1, 1 << 10, 30000, 30000)
check_constant(1, 1, 30000, 30000)
check_constant(1 << 10, 1, 10**9, 0)  # exact zeros
# Unbalanced inputs inside the bound |a| |b| < 4e14.
check_constant(1 << 17, 1 << 17, 10**9, 1)
check_constant(1 << 18, 1 << 10, 141000, 141000)
check_alternating(1 << 18, 1 << 10, 141000, 141000)
# Near-constant values at the bound, result length past 2^20: every
# entry against the exact product mod 998244353 (poly/ntt).
ntt = load("poly/ntt")
n, h = 600001, 25900
a, b = rnd_list(n, h - h // 100, h), rnd_list(n, h - h // 100, h)
got = [round(x) % ntt.MOD for x in lib.fft_conv(a, b)]
check_eq(got, ntt.conv(a, b), "near-constant, |a||b| < 4e14")
with Timer("fft_conv 2^18 x 2^18"):
  n = 1 << 18
  a, b = rnd_list(n, 0, 30000), rnd_list(n, 0, 30000)
  c = lib.fft_conv(a, b)
for _ in range(20):
  k = rnd(0, 2 * n - 2)
  lo, hi = max(0, k - n + 1), min(k, n - 1)
  check_eq(round(c[k]), sum(a[i] * b[k - i] for i in range(lo, hi + 1)))
print("OK poly/fft")
