import math

from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("misc/ternary-search")


def valley(n, flat):
  """Strictly down to a random minimum, then up (flat steps if flat)."""
  m = rnd(0, n - 1)
  f = [0] * n
  f[m] = rnd(-(10**12), 10**12)
  for i in range(m - 1, -1, -1):
    f[i] = f[i + 1] + rnd(1, 5)
  for i in range(m + 1, n):
    f[i] = f[i - 1] + rnd(0 if flat else 1, 5)
  return f


def convex(n):
  d = sorted(rnd(-6, 6) for _ in range(n - 1))
  f = [rnd(-99, 99)]
  for x in d:
    f.append(f[-1] + x)
  return f


def check_int(f):
  n, lo, calls = len(f), rnd(-(10**9), 10**9), [0]

  def g(x):
    check(lo <= x < lo + n, "called outside [lo, hi]")
    calls[0] += 1
    check(calls[0] <= 2 * (n - 1).bit_length(), "too many calls")
    return f[x - lo]

  check_eq(lib.ternary(lo, lo + n - 1, g) - lo, f.index(min(f)))


def check_real(a, b, f, cl, cr):
  """A minimiser up to what floats resolve: x within a few ulps of
  the minimisers [cl, cr], or f(x) within 1e-9 of the minimum."""
  calls = [0]

  def g(t):
    calls[0] += 1
    return f(t)

  x = lib.golden(a, b, g)
  check_eq(calls[0], 102)
  d = 1e-13 * max(1.0, abs(a), abs(b))
  fmin = f(cl)
  check(a - d <= x <= b + d, "outside [a, b]")
  check(
    cl - d <= x <= cr + d or f(x) <= fmin + 1e-9 * max(1.0, abs(fmin)),
    f"golden({a}, {b}) = {x}, minimisers [{cl}, {cr}]",
  )
  return x


def random_real():
  a = rng.uniform(-1e3, 1e3)
  b = a + rng.uniform(0, 1e3)
  if rnd(0, 3) == 0:
    a, b = a * 1e6, b * 1e6
  c, s = rng.uniform(a, b), rng.uniform(1e-3, 1e3)
  k, tol = rng.uniform(-1e6, 1e6), 1e-9 * max(1.0, b - a)
  case = rnd(0, 5)
  if case == 0:
    check_real(a, b, lambda x: s * (x - c) ** 2 + k, c, c)
  elif case == 1:
    x = check_real(a, b, lambda x: s * abs(x - c), c, c)
    check(abs(x - c) <= tol, "V shape")
  elif case == 2:  # linear on the left, cubic on the right
    f = lambda x: s * (c - x) + k if x < c else (x - c) ** 3 + k
    check_real(a, b, f, c, c)
  elif case == 3:  # flat bottom of half-width w
    w = rng.uniform(0, (b - a) / 2)
    f = lambda x: max(0.0, abs(x - c) - w)
    x = check_real(a, b, f, max(a, c - w), min(b, c + w))
    check(abs(x - c) <= w + tol, "flat bottom")
  elif case == 4:  # monotone: minimum at an endpoint
    check(check_real(a, b, lambda x: s * x, a, a) - a <= tol, "at a")
    check(b - check_real(a, b, lambda x: -s * x, b, b) <= tol, "at b")
  else:  # smooth, not a polynomial
    b = a + rng.uniform(0, 40)
    c = rng.uniform(a, b)
    f = lambda x: math.exp(x - c) + math.exp(c - x)
    check_real(a, b, f, c, c)


for _ in range(iters(2000)):
  n = rnd(1, 40)
  check_int(valley(n, False))
  check_int(valley(n, True))
  check_int(convex(n))
  random_real()

# Edges: one point, a constant function, a huge range, tuples.
check_eq(lib.ternary(5, 5, lambda x: 0), 5)
check_eq(lib.ternary(-3, 9, lambda x: 7), -3)
for _ in range(100):
  c = rnd(-(10**18), 10**18)
  big = 10**18
  check_eq(lib.ternary(-big, big, lambda x, c=c: abs(x - c)), c)
  check_eq(lib.ternary(-big, big, lambda x, c=c: (abs(x - c), 0)), c)
check_eq(lib.golden(2.5, 2.5, lambda x: x * x), 2.5)
for _ in range(100):  # golden with tuples, as in C++
  c = rng.uniform(-1e3, 1e3)
  x = lib.golden(-1e3, 1e3, lambda t, c=c: (abs(t - c), 0))
  check(abs(x - c) < 1e-9, "golden with tuples")

with Timer("1e5 integer searches over 1e9 values (6e6 calls)"):
  total = 0
  for _ in range(100_000):
    c = rnd(0, 10**9)
    f = lambda x, c=c: (x - c) * (x - c)
    total += lib.ternary(0, 10**9, f) - c
  check_eq(total, 0)
with Timer("5e4 golden searches (5e6 calls)"):
  for _ in range(50_000):
    c = rng.uniform(-1e3, 1e3)
    x = lib.golden(-1e3, 1e3, lambda t, c=c: abs(t - c))
    check(abs(x - c) < 1e-9, "golden speed case")
print("OK misc/ternary-search")
