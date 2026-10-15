from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("numtheory/primitive-root")


def is_prime_slow(n):
  return n >= 2 and all(n % d for d in range(2, int(n**0.5) + 1))


def brute_order(a, p):
  cur, e = a % p, 1
  while cur != 1:
    cur, e = cur * a % p, e + 1
  return e


def distinct_primes(n):
  f, d = [], 2
  while d * d <= n:
    if n % d == 0:
      f.append(d)
      while n % d == 0:
        n //= d
    d += 1
  return f + [n] if n > 1 else f


def generates(g, p, qs):
  return all(pow(g, (p - 1) // q, p) != 1 for q in qs)


def check_large(p, qs):
  """The answer generates and no smaller g does, given primes of p - 1."""
  g = lib.primitive_root(p)
  check(generates(g, p, qs), p)
  for h in range(2, g):
    check(not generates(h, p, qs), (p, h))


def order_example(a, p):
  e = p - 1
  for q in lib.factor(p - 1):
    if pow(a, e // q, p) == 1:
      e //= q
  return e


for p in range(2, 3000):
  if is_prime_slow(p):
    g = 1
    while brute_order(g, p) != p - 1:
      g += 1
    check_eq(lib.primitive_root(p), g, p)
    if p < 300:
      for a in range(1, p):
        check_eq(order_example(a, p), brute_order(a, p))
check_eq(lib.primitive_root(998244353), 3)
check_eq(lib.primitive_root(10**9 + 7), 5)
check_eq(lib.primitive_root(2**61 - 1), 37)
check_eq(lib.primitive_root(90441961), 113)
for _ in range(iters(20)):
  p = rnd(2, 10**11)
  while not lib.is_prime(p):
    p += 1
  check_large(p, distinct_primes(p - 1))
for p in [4179340454199820289, 2**63 - 25, 2**64 - 59, 10**18 + 3]:
  check_large(p, set(lib.factor(p - 1)))
check_large(2**127 - 1, set(lib.factor(2**127 - 2)))  # big ints

ps = []
while len(ps) < 300:
  p = rnd(1, 1 << 63)
  if lib.is_prime(p):
    ps.append(p)
with Timer("300 random 63-bit primes"):
  for p in ps:
    check_large(p, set(lib.factor(p - 1)))
print("OK numtheory/primitive-root")
