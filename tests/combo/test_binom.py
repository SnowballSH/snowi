from math import comb

from testlib import Timer, check_eq, iters, load, rnd

lib = load("combo/binom")
MOD = lib.MOD
P = 80


def pascal(m):
  c = [[0] * P for _ in range(P)]
  for n in range(P):
    c[n][0] = 1 % m
    for k in range(1, n + 1):
      c[n][k] = (c[n - 1][k - 1] + c[n - 1][k]) % m
  return c


def brute_catalan(count):
  cat = [1] + [0] * (count - 1)
  for m in range(1, count):
    cat[m] = sum(cat[i] * cat[m - 1 - i] for i in range(m)) % MOD
  return cat


def test_small():
  pas, cat = pascal(MOD), brute_catalan(P // 2)
  for size in (0, 1, 2, 5, P - 1):
    bc = lib.Binom(size)
    for n in range(-3, size + 1):
      for k in range(-3, size + 3):
        inside = 0 <= k <= n
        check_eq(bc.C(n, k), pas[n][k] if inside else 0)
        want = 0
        if inside:
          want = 1
          for i in range(k):
            want = want * (n - i) % MOD
        check_eq(bc.perm(n, k), want)
    for k in range(size // 2 + 1):
      check_eq(bc.catalan(k), cat[k])
    for i in range(1, size + 1):
      check_eq(bc.inv(i) * i % MOD, 1)
    for i in range(size + 1):
      check_eq(bc.fact[i] * bc.ifact[i] % MOD, 1)
    check_eq(bc.multinom([]), 1)
    for _ in range(iters(300)):
      ks, s = [], 0
      while True:
        k = rnd(0, 6)
        if s + k > size or rnd(0, 4) == 0:
          break
        ks.append(k)
        s += k
      want, s = 1, 0
      for k in ks:  # choose the slots of each part in turn
        s += k
        want = want * pas[s][k] % MOD
      check_eq(bc.multinom(ks), want)


def test_exact():
  n = 3000
  bc = lib.Binom(n)
  for _ in range(iters(3000)):
    a = rnd(0, n)
    k = rnd(-2, a + 2)
    want = comb(a, k) % MOD if 0 <= k <= a else 0
    check_eq(bc.C(a, k), want)


def test_lucas():
  for p in (2, 3, 5, 7, 13, 97, 1009):
    lu = lib.Lucas(p)
    pas = pascal(p)
    for i in range(p):
      check_eq(lu.fact[i] * lu.ifact[i] % p, 1)
    for n in range(-2, P):
      for k in range(-2, P + 2):
        check_eq(lu.C(n, k), pas[n][k] if 0 <= k <= n else 0)
    for _ in range(iters(300)):
      n, k = rnd(0, 10**18), rnd(0, 25)
      if rnd(0, 1):
        k = n - k
      if rnd(0, 3) == 0:
        n = rnd(0, 2000)
        k = rnd(-5, n + 5)
      want = comb(n, k) % p if 0 <= k <= n else 0
      check_eq(lu.C(n, k), want)


def test_speed():
  with Timer("Binom(2e6) + 2e5 C queries"):
    n = 2_000_000
    bc = lib.Binom(n)
    total = 0
    for i in range(200_000):
      total += bc.C(n - i, i)
  check_eq(bc.C(n, 7), comb(n, 7) % MOD)
  check_eq(bc.C(n, n // 2), bc.C(n, n - n // 2))
  with Timer("Lucas(1000003) + 1e5 queries"):
    p = 1_000_003
    lu = lib.Lucas(p)
    for _ in range(100_000):
      lu.C(rnd(0, 10**18), rnd(0, 10**18))
  for _ in range(iters(100)):
    n, k = rnd(0, 10**18), rnd(0, 20)
    check_eq(lu.C(n, k), comb(n, k) % p)


test_small()
test_exact()
test_lucas()
test_speed()
print("OK combo/binom")
