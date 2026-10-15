from math import gcd

from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("numtheory/discrete-log")


def brute_log(a, b, m):
  """a^x mod m takes at most m distinct values before it cycles."""
  cur, b = 1 % m, b % m
  for x in range(m + 1):
    if cur == b:
      return x
    cur = cur * a % m
  return -1


def is_prime_slow(n):
  return n >= 2 and all(n % d for d in range(2, int(n**0.5) + 1))


def tweak_logs(a, m, t):
  """The Tweak's many-query variant, for every b: one table of a^j,
  then giant steps by a^(-t)."""
  table, cur = {}, 1
  for j in range(t):
    table.setdefault(cur, j)  # the smallest j for each value
    cur = cur * a % m
  step, out = pow(a, -t, m), []
  for b in range(m):
    c, x = b, -1
    for i in range(m // t + 1):  # at most m / t giant steps
      if c in table:
        x = i * t + table[c]
        break
      c = c * step % m
    out.append(x)
  return out


def check_large(a, x, m):
  """The answer solves the equation and is at most a known x."""
  b = pow(a, x, m)
  got = lib.discrete_log(a, b, m)
  check(0 <= got <= x, (a, x, m, got))
  check_eq(pow(a, got, m), b)


for m in range(1, 45):
  for a in range(m):
    for b in range(m):
      check_eq(lib.discrete_log(a, b, m), brute_log(a, b, m), (a, b, m))
for _ in range(iters(2000)):
  m = rnd(1, 3000 if rnd(0, 1) else 300)
  if rnd(0, 2) == 0:
    m = 1 << rnd(0, 11)
  elif rnd(0, 2) == 0:
    m = rnd(1, 30) * rnd(1, 30) * 4
  a, b = rnd(0, m - 1), rnd(0, m - 1)
  if rnd(0, 1):
    b = pow(a, rnd(0, 2 * m), m)
  check_eq(lib.discrete_log(a, b, m), brute_log(a, b, m), (a, b, m))
# Unreduced and negative inputs.
check_eq(lib.discrete_log(20, 98, 17), 4)
check_eq(lib.discrete_log(-14, -4, 17), 4)
check_eq(lib.discrete_log(0, 0, 1), 0)
check_eq(lib.discrete_log(0, 1, 7), 0)
check_eq(lib.discrete_log(0, 0, 7), 1)
check_eq(lib.discrete_log(5, 0, 7), -1)
check_eq(lib.discrete_log(-1, 2, 3), 1)
check_eq(lib.discrete_log(-2, 2, 4), 1)
for m in [13, 16, 21, 27, 35]:  # the Tweak, coprime a only
  for a in range(1, m):
    if gcd(a, m) == 1:
      want = [lib.discrete_log(a, b, m) for b in range(m)]
      for t in range(1, m + 1):
        check_eq(tweak_logs(a, m, t), want, (a, m, t))
for m in range(1, 200):
  for a in range(1, m):
    if gcd(a, m) == 1:
      want = 1
      while pow(a, want, m) != 1 % m:
        want += 1
      check_eq(1 + lib.discrete_log(a, pow(a, -1, m), m), want)
p = 10**9 + 7  # 5 is a primitive root
check_eq(lib.discrete_log(5, p - 1, p), (p - 1) // 2)
check_eq(lib.discrete_log(25, 5, p), -1)
for _ in range(3):
  x = rnd(0, p - 2)
  check_eq(lib.discrete_log(5, pow(5, x, p), p), x)
for m in [2999999929, 3 * 10**9, 1 << 31, 998244353, 735134400]:
  for _ in range(3):
    check_large(rnd(0, m - 1), rnd(0, m), m)
for m in [10**12 + 39, 2**40]:  # above 3e9 as well
  check_large(rnd(2, m - 1), rnd(0, m), m)

primes = []
while len(primes) < 30:
  m = rnd(5 * 10**8, 10**9)
  if is_prime_slow(m):
    primes.append(m)
with Timer("60 logs mod ~1e9, half prime"):
  for m in primes + [rnd(5 * 10**8, 10**9) for _ in range(30)]:
    check_large(rnd(1, m - 1), rnd(0, m - 2), m)
print("OK numtheory/discrete-log")
