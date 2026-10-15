from testlib import CONTENT, Timer, check, check_eq, iters, load, rnd, rnd_list

lib = load("poly/subset-transforms")
MOD = lib.MOD


def brute_conv(a, b, op):  # the definition, O(4^n)
  c = [0] * len(a)
  for i, x in enumerate(a):
    for j, y in enumerate(b):
      k = i | j if op == "|" else i & j if op == "&" else i ^ j
      c[k] = (c[k] + x * y) % MOD
  return c


def brute_subset(a, b):
  c = [0] * len(a)
  for i, x in enumerate(a):
    for j, y in enumerate(b):
      if not i & j:
        c[i | j] = (c[i | j] + x * y) % MOD
  return c


def brute_transform(a, op):
  n, r = len(a), [0] * len(a)
  for S in range(n):
    for T in range(n):
      if op == "|" and T & ~S or op == "&" and S & ~T:
        continue
      odd = op == "^" and bin(S & T).count("1") % 2
      r[S] = (r[S] + (-a[T] if odd else a[T])) % MOD
  return r


for it in range(iters(500)):
  n = 1 << rnd(0, 4 if it % 5 else 6)
  hi = MOD - 1 if it % 3 else 2
  a, b = rnd_list(n, 0, hi), rnd_list(n, 0, hi)
  if it % 11 == 0:
    a, b = [MOD - 1] * n, [MOD - 1] * n
  for op in "|&^":
    t = a[:]
    lib.bit_transform(t, op)
    check_eq(t, brute_transform(a, op), op)
    lib.bit_transform(t, op, True)
    check_eq(t, a, op)
    check_eq(lib.bit_conv(a, b, op), brute_conv(a, b, op), op)
  check_eq(lib.subset_conv(a, b), brute_subset(a, b))


# The Tweak, applied to the snippet's own text: exact integers and max.
SOURCE = (CONTENT / "poly/subset-transforms.py").read_text()
TAIL = """  k = pow(n, MOD - 2, MOD) if inv and op == "^" else 1
  for i in range(n):
    a[i] = a[i] % MOD * k % MOD
"""


def tweaked(*edits):
  source = SOURCE
  for old, new in edits:
    check(old in source, f"Tweak target {old!r} is gone")
    source = source.replace(old, new)
  namespace = {"MOD": MOD}
  exec(compile(source, "tweaked", "exec"), namespace)
  return namespace


exact = tweaked(
  (TAIL, '  if inv and op == "^": a[:] = [x // n for x in a]\n'),
  ("x * y % MOD", "x * y"),
)
maxed = tweaked(
  ("a[j + s] += sg * a[j]", "a[j + s] = max(a[j + s], a[j])"),
  ("a[j] += sg * a[j + s]", "a[j] = max(a[j], a[j + s])"),
  (TAIL, ""),
)
for it in range(iters(200)):
  n, hi = 1 << rnd(0, 5), 4 * 10**9 if it % 2 else 3
  a, b = rnd_list(n, -hi, hi), rnd_list(n, -hi, hi)
  for op in "|&^":
    want = [0] * n
    for S in range(n):
      for T in range(n):
        if op == "^":
          want[S] += -a[T] if bin(S & T).count("1") % 2 else a[T]
        elif op == "|" and not T & ~S or op == "&" and not S & ~T:
          want[S] += a[T]
    t = a[:]
    exact["bit_transform"](t, op)
    check_eq(t, want, op)
    exact["bit_transform"](t, op, True)
    check_eq(t, a, op)
    c = [0] * n
    for i in range(n):
      for j in range(n):
        c[i | j if op == "|" else i & j if op == "&" else i ^ j] += a[i] * b[j]
    check_eq(exact["bit_conv"](a, b, op), c, op)
  for op in "|&":
    t = a[:]
    maxed["bit_transform"](t, op)
    for S in range(n):
      inside = [
        T for T in range(n) if (T | S == S if op == "|" else T & S == S)
      ]
      check_eq(t[S], max(a[T] for T in inside), op)

f = [0] * 7 + [1]
check_eq(lib.bit_conv(f, f, "&"), f)

n = 1 << 20
a = rnd_list(n, 0, MOD - 1)
with Timer("XOR transform and inverse, n = 2^20"):
  t = a[:]
  lib.bit_transform(t, "^")
  lib.bit_transform(t, "^", True)
check_eq(t, a)
n = 1 << 16
a, b = a[:n], rnd_list(n, 0, MOD - 1)
c = lib.bit_conv(a, b, "^")
for _ in range(20):
  S = rnd(0, n - 1)
  check_eq(c[S], sum(a[i] * b[i ^ S] for i in range(n)) % MOD)
with Timer("subset_conv n = 2^16"):
  c = lib.subset_conv(a, b)
for _ in range(20):
  S = T = rnd(0, n - 1)
  want = 0
  while True:
    want = (want + a[T] * b[S ^ T]) % MOD
    if T == 0:
      break
    T = (T - 1) & S
  check_eq(c[S], want)
print("OK poly/subset-transforms")
