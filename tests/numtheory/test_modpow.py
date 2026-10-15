from testlib import check_eq, iters, load, rnd

lib = load("numtheory/modpow")
MOD = lib.MOD

for _ in range(iters(2000)):
  b = rnd(1, MOD - 1)
  check_eq(b * lib.modinv(b) % MOD, 1)
  check_eq(lib.modinv(b), pow(b, -1, MOD))
print("OK numtheory/modpow")
