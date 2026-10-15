from math import isqrt

from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("numtheory/prime-count")

LIM = 1_000_000
is_p = bytearray([1]) * (LIM + 1)
is_p[0] = is_p[1] = 0
for i in range(2, isqrt(LIM) + 1):
  if is_p[i]:
    is_p[i * i :: i] = bytes(len(range(i * i, LIM + 1, i)))
pi_table = [0] * (LIM + 1)
for i in range(1, LIM + 1):
  pi_table[i] = pi_table[i - 1] + is_p[i]


def check_all(n):
  """Every value the class promises: v <= sqrt(n) and v = n // i."""
  pc = lib.PrimeCount(n)
  check_eq(pc.r, isqrt(n))
  for v in range(pc.r + 1):
    check_eq(pc(v), pi_table[v], (n, v))
  i = 1
  while i <= n:
    check_eq(pc(n // i), pi_table[n // i], (n, i))
    i = n // (n // i) + 1


for n in range(2000):
  check_all(n)
for _ in range(iters(200)):
  check_all(rnd(0, LIM))
for n in [LIM, LIM - 1, 1 << 19, 999_999]:
  check_all(n)
for r in range(2, 1001, 7):  # squares and their neighbours
  for n in (r * r - 1, r * r, r * r + 1):
    check_all(n)

# Example: numbers with exactly 4 divisors (p^3 or p q, p < q).
divs = [0] * 20001
for d in range(1, 20001):
  for k in range(d, 20001, d):
    divs[k] += 1
brute = 0
for n in range(1, 20001):
  brute += divs[n] == 4
  if n % 97 and n > 200:
    continue
  pc, ans = lib.PrimeCount(n), 0
  for p in range(2, pc.r + 1):
    if pc(p) > pc(p - 1):
      ans += (p**3 <= n) + pc(n // p) - pc(p)
  check_eq(ans, brute, n)


# Tweak: sum of primes, built as the header describes.
def prime_sum(n):
  r = isqrt(n)
  lo = [max(i * (i + 1) // 2 - 1, 0) for i in range(r + 1)]
  hi = [0] + [(n // i) * (n // i + 1) // 2 - 1 for i in range(1, r + 1)]
  for p in range(2, r + 1):
    if lo[p] == lo[p - 1]:
      continue
    s, q = lo[p - 1], p * p
    for i in range(1, min(r, n // q) + 1):
      d = i * p
      hi[i] -= p * ((hi[d] if d <= r else lo[n // d]) - s)
    for v in range(r, q - 1, -1):
      lo[v] -= p * (lo[v // p] - s)
  return hi[1] if n > r else lo[n]


prefix = 0
for n in range(1, 3000):
  prefix += n * is_p[n]
  check_eq(prime_sum(n), prefix, n)
check_eq(prime_sum(10**9), 24739512092254535)

with Timer("pi(1e10)"):
  check_eq(lib.PrimeCount(10**10)(10**10), 455052511)
with Timer("pi(1e11)"):
  pc = lib.PrimeCount(10**11)
  check_eq(pc(10**11), 4118054813)
  check_eq(pc(10**11 // 10), 455052511)
  check_eq(pc(10**11 // 100), 50847534)
check(pc(1) == 0 and pc(2) == 1)
print("OK numtheory/prime-count")
