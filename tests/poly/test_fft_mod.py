from testlib import Timer, check_eq, iters, load, rnd, rnd_list

lib = load("poly/fft-mod")


def brute_conv(a, b, m):
  if not a or not b:
    return []
  c = [0] * (len(a) + len(b) - 1)
  for i, x in enumerate(a):
    for j, y in enumerate(b):
      c[i + j] = (c[i + j] + x * y) % m
  return c


def check_constant(n, m, v):
  """All entries v. Worst for rounding: v near m / 2; v = m - 1 needs
  the balanced split."""
  c = lib.conv_mod([v] * n, [v] * n, m)
  sq = v * v % m
  check_eq(c, [min(k + 1, 2 * n - 1 - k) * sq % m for k in range(2 * n - 1)])


MODS = [10**9 + 7, 998244353, 2**31 - 1, 2**30 - 1, 2**30, 2, 3, 1009]
for it in range(iters(500)):
  m = MODS[rnd(0, len(MODS) - 1)] if it % 3 else rnd(2, 2**31 - 1)
  n = rnd(0, 50 if it % 3 else 120)
  k = rnd(0, 50 if it % 3 else 120)
  a, b = rnd_list(n, 0, m - 1), rnd_list(k, 0, m - 1)
  if it % 5 == 0:
    v = [0, m - 1, m // 2, m // 2 + 1]
    a, b = [v[it % 4]] * n, [v[it // 5 % 4]] * k
  check_eq(lib.conv_mod(a, b, m), brute_conv(a, b, m), f"m={m}")

# Limits from the header: length 2^20 below 2^30, 2^19 below 2^31.
with Timer("conv_mod 2^19 x 2^19, worst case"):
  check_constant(1 << 19, 2**30 - 1, 2**29 - 1)
check_constant(1 << 18, 2**31 - 1, 2**30)
check_constant(1 << 18, 2**31 - 1, 2**31 - 2)
with Timer("conv_mod 2^18 x 2^18"):
  n, m = 1 << 18, 10**9 + 7
  a, b = rnd_list(n, 0, m - 1), rnd_list(n, 0, m - 1)
  c = lib.conv_mod(a, b, m)
for _ in range(20):
  k = rnd(0, 2 * n - 2)
  lo, hi = max(0, k - n + 1), min(k, n - 1)
  want = sum(a[i] * b[k - i] for i in range(lo, hi + 1)) % m
  check_eq(c[k], want)
print("OK poly/fft-mod")
