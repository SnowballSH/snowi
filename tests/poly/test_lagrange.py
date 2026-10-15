from testlib import Timer, check_eq, iters, load, rnd, rnd_list, rng

lib = load("poly/lagrange")
MOD = lib.MOD


def horner(c, x):
  r = 0
  for v in reversed(c):
    r = (r * x + v) % MOD
  return r


for it in range(iters(1500)):
  n = rnd(1, 12 if it % 4 else 40)
  c = [MOD - 1] * n if it % 9 == 0 else rnd_list(n, 0, MOD - 1)
  y = [horner(c, i) for i in range(n)]
  for x in (rnd(0, n - 1), rnd(0, MOD - 1), rnd(MOD, 10**30), -rnd(0, 10**30)):
    check_eq(lib.lagrange(y, x), horner(c, x), f"x={x}")
  xs = [0] if it % 3 == 0 else []
  if it % 5 == 0:
    xs.append(MOD - 1)
  xs = xs[:n]
  while len(xs) < n:
    x = rnd(0, MOD - 1) if it % 2 else rnd(0, 3 * n)
    if x not in xs:
      xs.append(x)
  rng.shuffle(xs)
  check_eq(lib.interpolate(xs, [horner(c, x) for x in xs]), c)

check_eq(lib.interpolate([], []), [])
check_eq(lib.lagrange([7], 123456789012), 7)
check_eq(lib.lagrange([0, 1, 5, 14], 10), 385)
for k in range(8):
  y = [0] * (k + 2)
  for i in range(1, k + 2):
    y[i] = (y[i - 1] + pow(i, k, MOD)) % MOD
  brute = 0
  for m in range(60):
    brute = (brute + pow(m, k, MOD)) % MOD if m else 0
    check_eq(lib.lagrange(y, m), brute)

d = 10**6
y = [pow(i + 5, d, MOD) for i in range(d + 1)]
with Timer("lagrange d = 1e6"):
  x = rnd(0, 10**18)
  got = lib.lagrange(y, x)
check_eq(got, pow(x + 5, d, MOD))
n = 2000
c, xs = rnd_list(n, 0, MOD - 1), [i * i for i in range(n)]
ys = [horner(c, x) for x in xs]
with Timer("interpolate n = 2000"):
  got = lib.interpolate(xs, ys)
check_eq(got, c)
print("OK poly/lagrange")
