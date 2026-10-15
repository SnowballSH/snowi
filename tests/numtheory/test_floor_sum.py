from math import gcd

from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("numtheory/floor-sum")
floor_sum = lib.floor_sum


def brute(n, m, a, b):
  return sum((a * i + b) // m for i in range(n))


for _ in range(iters(20000)):
  n, m = rnd(0, 30), rnd(1, 30)
  a, b = rnd(-60, 60), rnd(-60, 60)
  check_eq(floor_sum(n, m, a, b), brute(n, m, a, b))
for _ in range(iters(2000)):
  n, m = rnd(0, 500), rnd(1, 1 << rnd(0, 62))
  a, b = rnd(-(1 << 62), 1 << 62), rnd(-(1 << 70), 1 << 70)
  check_eq(floor_sum(n, m, a, b), brute(n, m, a, b))
# If gcd(a, m) = 1: floor_sum(m, m, a, b) = (a - 1)(m - 1) / 2 + b.
for _ in range(iters(2000)):
  m, a = rnd(1, 1 << rnd(1, 80)), rnd(-(1 << 80), 1 << 80)
  if gcd(a, m) == 1:
    b = rnd(-(1 << 90), 1 << 90)
    check_eq(floor_sum(m, m, a, b), (a - 1) * (m - 1) // 2 + b)
check_eq(floor_sum(0, 5, 3, 2), 0)
check_eq(floor_sum(1, 1, -7, -9), -9)
check_eq(floor_sum(10**18, 1, 0, 10**18), 10**36)

for _ in range(iters(2000)):  # Tweak: sum over i in [l, r)
  m, a, b = rnd(1, 20), rnd(-20, 20), rnd(-20, 20)
  l = rnd(-10, 10)
  r = rnd(l, 20)
  want = sum((a * i + b) // m for i in range(l, r))
  check_eq(floor_sum(r - l, m, a, a * l + b), want)

with Timer("2e5 calls, n, m < 2^32"):
  acc = 0
  for _ in range(200_000):
    m, n = rnd(1, (1 << 32) - 1), rnd(0, (1 << 32) - 1)
    acc ^= floor_sum(n, m, rnd(-m, m), rnd(-m, m))
  check(acc != 42)
print("OK numtheory/floor-sum")
