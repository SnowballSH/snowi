from fractions import Fraction
from itertools import combinations

from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("linalg/simplex")
LPSolver = lib.LPSolver
INF = float("inf")


class Capped(LPSolver):
  """Fails the test on cycling instead of hanging."""

  pivots = 0

  def pivot(self, r, s):
    self.pivots += 1
    check(self.pivots < 1000, "too many pivots: cycling")
    LPSolver.pivot(self, r, s)


def solve_exact(rows, rhs):
  """Unique solution of the square system, or None (Fractions)."""
  n = len(rows)
  a = [[Fraction(v) for v in row] + [Fraction(r)] for row, r in zip(rows, rhs)]
  for c in range(n):
    p = next((i for i in range(c, n) if a[i][c]), None)
    if p is None:
      return None
    a[c], a[p] = a[p], a[c]
    for i in range(n):
      if i != c and a[i][c]:
        f = a[i][c] / a[c][c]
        a[i] = [u - f * v for u, v in zip(a[i], a[c])]
  return [a[i][n] / a[i][i] for i in range(n)]


def vertex_opt(A, b, c, box):
  """Exact optimum over the vertices of {Ax <= b, 0 <= x <= box}."""
  n = len(c)
  rows, rhs = [list(r) for r in A], list(b)
  for j in range(n):
    rows.append([-int(i == j) for i in range(n)])
    rhs.append(0)
    rows.append([int(i == j) for i in range(n)])
    rhs.append(box)
  best = None
  for tight in combinations(range(len(rows)), n):
    x = solve_exact([rows[i] for i in tight], [rhs[i] for i in tight])
    if x is None:
      continue
    if all(
      sum(u * v for u, v in zip(row, x)) <= r for row, r in zip(rows, rhs)
    ):
      val = sum(u * v for u, v in zip(c, x))
      if best is None or val > best:
        best = val
  return best


def check_solution(lp, A, b, c, val, x, exact, tol=1e-6):
  """Feasible, optimal and dual-certified (strong duality)."""
  m, n = len(b), len(c)

  def close(u, v):
    return u == v if exact else abs(u - v) <= tol * (1 + abs(u) + abs(v))

  def le(u, v):
    return u <= v or close(u, v)

  check_eq(len(x), n)
  check(all(le(0, v) for v in x))
  check(close(sum(u * v for u, v in zip(c, x)), val))
  for row, r in zip(A, b):
    check(le(sum(u * v for u, v in zip(row, x)), r))
  y = [0] * m
  for j, var in enumerate(lp.N):
    if var >= n:
      y[var - n] = lp.D[m][j]
  check(all(le(0, v) for v in y))
  check(close(sum(u * v for u, v in zip(b, y)), val))
  for j in range(n):
    check(le(c[j], sum(A[i][j] * y[i] for i in range(m))))


for _ in range(iters(1500)):
  n, m = rnd(1, 3), rnd(0, 4)
  A = [[rnd(-4, 4) if rnd(0, 2) else 0 for _ in range(n)] for _ in range(m)]
  b = [rnd(-6, 6) if rnd(0, 2) else 0 for _ in range(m)]
  c = [rnd(-4, 4) for _ in range(n)]
  if m > 1 and rnd(0, 3) == 0:
    A[1], b[1] = A[0][:], b[0]
  small, big = vertex_opt(A, b, c, 1000), vertex_opt(A, b, c, 2000)
  want = -INF if small is None else INF if big > small else small
  for exact in (False, True):  # Tweak: Fractions with eps = 0
    conv = Fraction if exact else float
    Ac = [[conv(v) for v in row] for row in A]
    bc, cc = [conv(v) for v in b], [conv(v) for v in c]
    lp = Capped(Ac, bc, cc, 0 if exact else 1e-9)
    val, x = lp.solve()
    if abs(want) == INF:
      check_eq((val, x), (want, None))
    else:
      check(val == want if exact else abs(val - want) < 1e-6)
      check_solution(lp, A, b, c, val, x, exact)

# Beale's example: Dantzig's rule alone cycles on it. In exact
# mode (eps = 0) a degenerate step has ratio exactly 0 = eps, and
# that must switch to Bland's rule too.
F = Fraction
A = [[F(1, 4), -60, F(-1, 25), 9], [F(1, 2), -90, F(-1, 50), 3], [0, 0, 1, 0]]
b, c = [0, 0, 1], [F(3, 4), -150, F(1, 50), -6]
lp = Capped([[float(v) for v in row] for row in A], b, [float(v) for v in c])
val, x = lp.solve()
check(abs(val - 0.05) < 1e-9)
check_solution(lp, A, b, c, val, x, False)
lp = Capped(A, [F(v) for v in b], c, 0)
val, x = lp.solve()
check_eq((val, x), (F(1, 20), [F(1, 25), 0, 1, 0]))
check_solution(lp, A, b, c, val, x, True)
check_eq(LPSolver([], [], [1]).solve(), (INF, None))
check_eq(LPSolver([], [], [-1, 0]).solve(), (0, [0, 0]))
check_eq(LPSolver([[0]], [-1], [1]).solve(), (-INF, None))

# Klee-Minty cube: Dantzig's rule visits 2^n - 1 vertices.
km = 10
K = [
  [2 ** (i - j + 1) if j < i else int(i == j) for j in range(km)]
  for i in range(km)
]
kb = [5 ** (i + 1) for i in range(km)]
kc = [2 ** (km - 1 - i) for i in range(km)]
val, x = LPSolver(K, kb, kc).solve()
check(abs(val / 5**km - 1) < 1e-9)


def speed(n, m):
  R = [[rng.uniform(-1, 1) for _ in range(n)] for _ in range(m)]
  x0 = [rng.uniform(0, 1) for _ in range(n)]
  b = [sum(u * v for u, v in zip(row, x0)) + rng.uniform(0, 1) for row in R]
  c = [rng.uniform(-1, 1) for _ in range(n)]
  for j in range(n):  # keep it bounded: x_j <= 2
    R.append([float(i == j) for i in range(n)])
    b.append(2.0)
  with Timer(f"random LP, {n} variables, {len(b)} constraints"):
    lp = LPSolver(R, b, c)
    val, x = lp.solve()
  check(abs(val) != INF)
  check_solution(lp, R, b, c, val, x, False)


speed(100, 100)
print("OK linalg/simplex")
