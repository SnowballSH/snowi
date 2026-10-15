from math import gcd

from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("numtheory/euclid")


def brute_inv(a, m):
  for x in range(m):
    if a * x % m == 1 % m:
      return x
  return -1


def brute_crt(rs, ms):
  l = 1
  for m in ms:
    l = l * m // gcd(l, m)
  for x in range(l):
    if all(x % m == r % m for r, m in zip(rs, ms)):
      return x, l
  return -1, -1


def check_euclid(a, b):
  g, x, y = lib.euclid(a, b)
  check_eq(g, gcd(a, b))
  check_eq(a * x + b * y, g)
  if g:
    check(abs(x) <= max(b // g, 1) and abs(y) <= max(a // g, 1))


for _ in range(iters(5000)):
  check_euclid(rnd(0, 30), rnd(0, 30))
  check_euclid(rnd(0, 10**30), rnd(0, 10**30))
  g = rnd(1, 1000)
  check_euclid(g * rnd(0, 10**18), g * rnd(0, 10**18))
check_euclid(0, 0)

for m in range(1, 50):
  for a in range(-m, 2 * m):
    check_eq(lib.inv_mod(a, m), brute_inv(a % m, m))
for _ in range(iters(5000)):
  m, a = rnd(1, 10**30), rnd(-(10**30), 10**30)
  v = lib.inv_mod(a, m)
  if gcd(a, m) != 1:
    check_eq(v, -1)
  else:
    check(0 <= v < m and a * v % m == 1 % m)

for m in range(1, 13):
  for n in range(1, 13):
    for a in range(-m, m):
      for b in range(n):
        check_eq(lib.crt(a, m, b, n), brute_crt([a, b], [m, n]))
for _ in range(iters(2000)):
  k = rnd(0, 4)
  ms = [rnd(1, 10) for _ in range(k)]
  rs = [rnd(-20, 20) for _ in range(k)]
  check_eq(lib.crt_all(rs, ms), brute_crt(rs, ms))
for _ in range(iters(5000)):
  g = rnd(1, 10**9)
  m, n = g * rnd(1, 10**12), g * rnd(1, 10**12)
  g = gcd(m, n)
  a, b = rnd(0, m - 1), rnd(0, n - 1)
  if rnd(0, 1):
    b = (a + g * rnd(0, n)) % n
  x, l = lib.crt(a, m, b, n)
  if (b - a) % g:
    check_eq((x, l), (-1, -1))
  else:
    check_eq(l, m // g * n)
    check(0 <= x < l and x % m == a and x % n == b)
check_eq(lib.crt_all([], []), (0, 1))
check_eq((lib.inv_mod(3, 10), lib.crt(2, 4, 4, 6)), (7, (10, 12)))
check_eq(lib.crt_all([1, 2, 3], [2, 3, 5]), (23, 30))
check_eq(lib.crt_all([2, 3], [4, 6]), (-1, -1))

with Timer("2e5 crt + inv_mod, 1e9 to 1e18 moduli"):
  acc = 0
  for _ in range(100_000):
    m, n = rnd(1, 10**9), rnd(1, 10**9)
    acc ^= lib.crt(rnd(0, m - 1), m, rnd(0, n - 1), n)[1]
    acc ^= lib.inv_mod(rnd(0, 10**18), 10**18 + 9)
  check(acc != 42)
print("OK numtheory/euclid")
