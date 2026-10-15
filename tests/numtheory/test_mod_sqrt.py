from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("numtheory/mod-sqrt")


def is_prime_slow(n):
  return n >= 2 and all(n % d for d in range(2, int(n**0.5) + 1))


def check_prime(p):
  """Squares give back a root; a non-square times a square gives -1."""
  y = rnd(0, p - 1)
  r = lib.mod_sqrt(y * y, p)
  check_eq(r, min(y, p - y), p)
  non = 2
  while pow(non, (p - 1) // 2, p) != p - 1:
    non += 1
  if y:
    check_eq(lib.mod_sqrt(non * y * y, p), -1, p)


for p in range(2, 700):
  if not is_prime_slow(p):
    continue
  root = [-1] * p
  for x in range(p - 1, -1, -1):
    root[x * x % p] = x
  for a in range(p):
    check_eq(lib.mod_sqrt(a, p), root[a], (a, p))
  check_eq(lib.mod_sqrt(p + 1, p), 1)
  check_eq(lib.mod_sqrt(-p, p), 0)
  check_eq(lib.mod_sqrt(-1, p), root[p - 1])
# Primes k 2^s + 1 with large s are the slowest case.
for p in [998244353, 469762049, 2013265921, 167772161, 2**61 - 1]:
  for _ in range(100):
    check_prime(p)
for p in [4179340454199820289, 2**63 - 25, 2**64 - 59, 2**127 - 1]:
  for _ in range(20):
    check_prime(p)
for _ in range(iters(1000)):
  p = rnd(2, 10**6 if rnd(0, 1) else 3 * 10**9)
  while not is_prime_slow(p):
    p += 1
  check_prime(p)
check_eq(lib.mod_sqrt(2, 998244353), 116195171)

with Timer("2e5 roots mod 998244353"):
  p, found = 998244353, 0
  for _ in range(200_000):
    a = rnd(0, p - 1)
    r = lib.mod_sqrt(a, p)
    if r >= 0:
      found += 1
      check_eq(r * r % p, a)
  check(98_000 < found < 102_000)
print("OK numtheory/mod-sqrt")
