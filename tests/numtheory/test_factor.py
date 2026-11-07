from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("numtheory/factor")


def oracle_prime(n):
  """Independent Miller-Rabin with 24 random bases."""
  if n < 4:
    return n >= 2
  if n % 2 == 0:
    return False
  d, s = n - 1, 0
  while d % 2 == 0:
    d //= 2
    s += 1
  for _ in range(24):
    x = pow(rnd(2, n - 2), d, n)
    ok = x in (1, n - 1)
    for _ in range(s - 1):
      if ok:
        break
      x = x * x % n
      ok = x == n - 1
    if not ok:
      return False
  return True


def trial_factor(x):
  f, d = [], 2
  while d * d <= x:
    while x % d == 0:
      f.append(d)
      x //= d
    d += 1
  if x > 1:
    f.append(x)
  return f


def random_prime(lo, hi):
  while True:
    p = rnd(lo, hi)
    if oracle_prime(p):
      return p


def check_factor(n):
  f = lib.factor(n)
  check(f == sorted(f))
  prod = 1
  for p in f:
    check(oracle_prime(p) and lib.is_prime(p))
    prod *= p
  check_eq(prod, n)


n = 100_000
composite = bytearray(n + 1)
for i in range(2, n + 1):
  check_eq(lib.is_prime(i), not composite[i])
  if not composite[i]:
    composite[i * i :: i] = b"\1" * len(range(i * i, n + 1, i))
check(lib.is_prime(10**9 + 7) and lib.factor(360) == [2, 2, 2, 3, 3, 5])
check(not lib.is_prime(0) and not lib.is_prime(1) and not lib.is_prime(-7))
pseudoprimes = [2047, 3215031751, 4759123141, 1122004669633, 3474749660383]
pseudoprimes += [341550071728321, 3825123056546413051, 561, 41041]
for n in pseudoprimes + [825265, 321197185, 5394826801]:
  check(not lib.is_prime(n))
  check_factor(n)
# Strong pseudoprimes to all bases but one, in order of the missing
# base: dropping any base must fail here.
almost = [1921077350011, 1411807385341, 443538368977861, 4341937413061]
almost += [5517315475561, 3933464309633, 107528788110061]
for n in almost:
  check(not lib.is_prime(n))
for _ in range(iters(5000)):
  n = rnd(1, (1 << 64) >> rnd(0, 63))
  check_eq(lib.is_prime(n), oracle_prime(n))
for _ in range(iters(200)):
  p = random_prime(2, (1 << 64) >> rnd(0, 62))
  check(lib.is_prime(p))
  check(not lib.is_prime(p * p))

for n in range(1, 3000):
  check_eq(lib.factor(n), trial_factor(n))
for _ in range(iters(300)):
  n = rnd(1, 10**9)
  check_eq(lib.factor(n), trial_factor(n))
for _ in range(iters(500)):
  check_factor(rnd(1, (1 << 63) >> rnd(0, 62)))
for n in [(1 << 63) - 1, (1 << 64) - 1, 1 << 62, 3**39, 10**18, 10**18 + 1]:
  check_factor(n)
check_factor(897612484786617600)
for _ in range(iters(10)):  # prime squares and cubes
  p, q = random_prime(2 * 10**9, 3 * 10**9), random_prime(10**6, 2 * 10**6)
  check_eq(lib.factor(p * p), [p, p])
  check_eq(lib.factor(q * q * q), [q, q, q])

with Timer("1e4 is_prime on 64-bit n"):
  cnt = sum(lib.is_prime(rnd(1, 1 << 64)) for _ in range(10_000))
  check(100 < cnt < 500)
with Timer("1e3 factor of random n <= 1e18"):
  for _ in range(1000):
    check_factor(rnd(1, 10**18))
ps = [random_prime(9 * 10**8, 10**9) for _ in range(50)]
qs = [random_prime(10**9, 11 * 10**8) for _ in range(50)]
with Timer("50 semiprimes p q with p, q near 1e9"):
  for p, q in zip(ps, qs):
    check_eq(lib.factor(p * q), [p, q])
# Balanced semiprimes at the top of the range, the worst case: below
# 2^63, and up to 2^64 where every value is a big int.
for lo, hi, label in [
  (3 * 10**9, 3037000499, "near 2^63"),
  (4 * 10**9, (1 << 32) - 1, "in [2^63, 2^64)"),
]:
  ps = [random_prime(lo, hi) for _ in range(40)]
  with Timer("20 semiprimes " + label):
    for p, q in zip(ps[::2], ps[1::2]):
      check_eq(lib.factor(p * q), sorted([p, q]))
print("OK numtheory/factor")
