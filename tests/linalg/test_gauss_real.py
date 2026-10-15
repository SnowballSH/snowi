from fractions import Fraction
from itertools import permutations, takewhile

from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("linalg/gauss-real")
gauss_real, det_real, solve_real = lib.gauss_real, lib.det_real, lib.solve_real


def perm_det(a):
  n, res = len(a), 0
  for perm in permutations(range(n)):
    term = 1
    for i in range(n):
      term *= a[i][perm[i]]
    inv = sum(perm[i] > perm[j] for i in range(n) for j in range(i + 1, n))
    res += -term if inv % 2 else term
  return res


def sub(a, rows, cols):
  return [
    [row[j] for j in range(len(row)) if cols >> j & 1]
    for i, row in enumerate(a)
    if rows >> i & 1
  ]


def minor_rank(a, m):
  n, best = len(a), 0
  for rows in range(1 << n):
    for cols in range(1 << m):
      k = bin(rows).count("1")
      if (
        k == bin(cols).count("1") and k > best and perm_det(sub(a, rows, cols))
      ):
        best = k
  return best


def cramer(a, b, m):
  """Exact unique solution from m independent rows."""
  for rows in range(1 << len(a)):
    if bin(rows).count("1") != m:
      continue
    s = sub(a, rows, (1 << m) - 1)
    d = perm_det(s)
    if d:
      bs = [b[i] for i in range(len(a)) if rows >> i & 1]
      x = []
      for j in range(m):
        t = [row[:j] + [bs[i]] + row[j + 1 :] for i, row in enumerate(s)]
        x.append(Fraction(perm_det(t), d))
      return x
  return None


def random_system():
  n, m, hi = rnd(0, 4), rnd(0, 4), rnd(1, 4)
  a = [[rnd(-hi, hi) for _ in range(m)] for _ in range(n)]
  if n > 1 and rnd(0, 2) == 0:  # dependent row
    i, j, f = rnd(0, n - 1), rnd(0, n - 1), rnd(-2, 2)
    a[i] = [f * v for v in a[j]]
  x0 = [rnd(-3, 3) for _ in range(m)]
  b = [sum(u * v for u, v in zip(row, x0)) for row in a]
  b = [v + (rnd(-2, 2) if rnd(0, 2) == 0 else 0) for v in b]
  return n, m, a, b


for _ in range(iters(2500)):
  n, m, a, b = random_system()
  ra = minor_rank(a, m)
  rab = minor_rank([row + [v] for row, v in zip(a, b)], m + 1)
  want_kind = 0 if ra < rab else 2 if ra < m else 1
  for exact in (False, True):  # Tweak: Fraction entries, eps = 0
    conv = Fraction if exact else float
    eps = 0 if exact else 1e-9
    af = [[conv(v) for v in row] for row in a]
    bf = [conv(v) for v in b]
    piv, d = gauss_real([row[:] for row in af], m, eps)
    check_eq(len(piv), ra)
    if n == m:
      want = perm_det(a)
      if exact:
        check_eq((d, det_real(af, eps)), (want, want))
      else:
        check(abs(d - want) < 1e-9 * max(1, abs(want)))
        check(abs(det_real(af, eps) - want) < 1e-9 * max(1, abs(want)))
    kind, x = solve_real(af, bf, m, eps)
    check_eq(kind, want_kind)
    if kind == 0:
      check_eq(x, None)
      continue
    check_eq(len(x), m)
    for row, v in zip(a, b):
      s = sum(u * w for u, w in zip(row, x))
      check(s == v if exact else abs(s - v) < 1e-9)
    if kind == 1:
      want = cramer(a, b, m)
      check(
        all(xv == w if exact else abs(xv - w) < 1e-9 for xv, w in zip(x, want))
      )

check_eq(solve_real([], [], 2), (2, [0.0, 0.0]))
check_eq(solve_real([[], []], [0.0, 1.0], 0), (0, None))
check_eq((det_real([]), det_real([[0.0, 1.0], [1.0, 0.0]])), (1, -1.0))
check_eq(det_real([[1.0, 2.0], [2.0, 4.0]]), 0)
tiny = Fraction(1, 10**10)  # exact mode reaches det_real too
check_eq(det_real([[tiny, 0], [0, 1]], 0), tiny)
check_eq(det_real([[1e-10, 0.0], [0.0, 1.0]], 1e-12), 1e-10)


def run_usage(**names):
  """Run the header's Usage lines in order, exactly as printed."""
  body = lib.__doc__.split("\nUsage:\n", 1)[1].split("\n")
  lines = takewhile(lambda line: line.startswith("  "), body)
  ns = dict(vars(lib), **names)
  exec("\n".join(line[2:] for line in lines), ns)
  return ns


# Usage top to bottom on one matrix: every call before the last
# one (gauss_real, in place) must still see the original matrix.
cases = [([[2, 1], [1, 3]], [3, 4])]
for _ in range(iters(300)):
  n, m, a, b = random_system()
  if n == m:
    cases.append((a, b))
for a, b in cases:
  n = len(a)
  u = run_usage(a=[[float(v) for v in row] for row in a], b=b, m=n)
  ra = minor_rank(a, n)
  rab = minor_rank([row + [v] for row, v in zip(a, b)], n + 1)
  check_eq(u["k"], 0 if ra < rab else 2 if ra < n else 1)
  if u["k"]:
    for row, v in zip(a, b):
      check(abs(sum(p * q for p, q in zip(row, u["x"])) - v) < 1e-9)
  check(abs(u["d"] - perm_det(a)) < 1e-9 * max(1, abs(perm_det(a))))
  check_eq(len(u["piv"]), ra)

for _ in range(iters(1000)):  # real entries: residuals
  n = rnd(1, 6)
  a = [[rng.uniform(-1e3, 1e3) for _ in range(n)] for _ in range(n)]
  b = [rng.uniform(-1e3, 1e3) for _ in range(n)]
  kind, x = solve_real(a, b, n)
  check_eq(kind, 1)
  for row, v in zip(a, b):
    check(abs(sum(u * w for u, w in zip(row, x)) - v) < 1e-6)

for n in (300, 500):
  a = [[rng.uniform(-1, 1) for _ in range(n)] for _ in range(n)]
  x0 = [rng.uniform(-10, 10) for _ in range(n)]
  b = [sum(u * v for u, v in zip(row, x0)) for row in a]
  with Timer(f"{n} x {n} solve"):
    kind, x = solve_real(a, b, n)
  check_eq(kind, 1)
  check(max(abs(u - v) for u, v in zip(x, x0)) < 1e-6)
print("OK linalg/gauss-real")
