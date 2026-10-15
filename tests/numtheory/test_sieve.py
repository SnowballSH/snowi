from math import gcd

from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("numtheory/sieve")


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


def brute_mu(x):
  f = trial_factor(x)
  if len(set(f)) < len(f):
    return 0
  return -1 if len(f) % 2 else 1


def check_sieve(n):
  s = lib.Sieve(n)
  check_eq(len(s.spf), n + 1)
  check_eq(s.primes, [x for x in range(2, n + 1) if len(trial_factor(x)) == 1])
  check_eq((s.spf[0], s.spf[1], s.phi[0], s.mu[0]), (0, 0, 0, 0))
  for x in range(1, n + 1):
    f = trial_factor(x)
    check_eq(s.factorize(x), f)
    if x > 1:
      check_eq(s.spf[x], f[0])
    check_eq(s.phi[x], sum(1 for k in range(1, x + 1) if gcd(k, x) == 1))
    check_eq(s.mu[x], brute_mu(x))


for n in range(1, 40):
  check_sieve(n)
for n in (997, 1000, 1024):
  check_sieve(n)

check_eq(lib.Sieve(1000).factorize(360), [2, 2, 2, 3, 3, 5])
s = lib.Sieve(300)
for _ in range(iters(100)):
  N, M = rnd(1, 300), rnd(1, 60)
  cnt = sum(s.mu[d] * (N // d) * (M // d) for d in range(1, min(N, M) + 1))
  want = sum(
    1 for i in range(1, N + 1) for j in range(1, M + 1) if gcd(i, j) == 1
  )
  check_eq(cnt, want)

for n in (10**6, 10**7):
  with Timer(f"Sieve({n:.0e})"):
    big = lib.Sieve(n)
  # sum_d mu(d) (n/d) = 1 and sum_d phi(d) (n/d) = n (n+1) / 2
  mu, phi = big.mu, big.phi
  check_eq(sum(mu[d] * (n // d) for d in range(1, n + 1)), 1)
  check_eq(sum(phi[d] * (n // d) for d in range(1, n + 1)), n * (n + 1) // 2)
  check_eq(len(big.primes), {10**6: 78498, 10**7: 664579}[n])
  check(big.factorize(n - 1) == trial_factor(n - 1))
  del big, mu, phi
print("OK numtheory/sieve")
