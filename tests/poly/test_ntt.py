from testlib import Timer, check_eq, iters, load, rnd, rnd_list

lib = load("poly/ntt")
MOD = lib.MOD


def brute_conv(a, b):
  if not a or not b:
    return []
  c = [0] * (len(a) + len(b) - 1)
  for i, x in enumerate(a):
    for j, y in enumerate(b):
      c[i + j] = (c[i + j] + x * y) % MOD
  return c


def brute_dft_bitrev(a):
  """a(w^k) at position rev(k), w = 3^((MOD - 1) / n)."""
  n = len(a)
  bits = n.bit_length() - 1
  w = pow(3, (MOD - 1) // n, MOD)
  res = [0] * n
  for k in range(n):
    x = pow(w, k, MOD)
    val = sum(c * pow(x, j, MOD) for j, c in enumerate(a)) % MOD
    rev = int(format(k, "b").zfill(bits)[::-1], 2) if bits else 0
    res[rev] = val
  return res


for it in range(iters(600)):
  n = rnd(0, 50 if it % 3 else 120)
  m = rnd(0, 50 if it % 3 else 120)
  hi = MOD - 1 if it % 5 else 3
  a, b = rnd_list(n, 0, hi), rnd_list(m, 0, hi)
  if it % 7 == 0:
    a, b = [MOD - 1] * n, [MOD - 1] * m
  check_eq(lib.conv(a, b), brute_conv(a, b))
  # Correlation from the header: sum_i a_i b_(i + d) at
  # c[len(a) - 1 + d] of conv(reversed a, b).
  c = lib.conv(a[::-1], b)
  for d in range(m - n + 1 if n else 0):
    want = sum(a[i] * b[i + d] for i in range(n)) % MOD
    check_eq(c[n - 1 + d], want)

for lg in range(8):
  n = 1 << lg
  a = rnd_list(n, 0, MOD - 1)
  f = a[:]
  lib.ntt(f)
  check_eq(f, brute_dft_bitrev(a), f"n={n}")
  lib.intt(f)
  check_eq(f, a, f"round trip n={n}")
  b = rnd_list(n, 0, MOD - 1)
  g = b[:]
  lib.ntt(f)
  lib.ntt(g)
  h = [x * y % MOD for x, y in zip(f, g)]
  lib.intt(h)
  full = brute_conv(a, b)
  cyc = [sum(full[k::n]) % MOD for k in range(n)]
  check_eq(h, cyc, f"cyclic n={n}")

with Timer("conv 2^18 x 2^18"):
  n = 1 << 18
  a, b = rnd_list(n, 0, MOD - 1), rnd_list(n, 0, MOD - 1)
  c = lib.conv(a, b)
for _ in range(20):
  k = rnd(0, 2 * n - 2)
  lo, hi = max(0, k - n + 1), min(k, n - 1)
  want = sum(a[i] * b[k - i] for i in range(lo, hi + 1)) % MOD
  check_eq(c[k], want)
n, m = 1 << 16, (1 << 16) + 3
c = lib.conv([1] * n, [1] * m)
check_eq(c, [min(k + 1, n, m, n + m - 1 - k) for k in range(n + m - 1)])
print("OK poly/ntt")
