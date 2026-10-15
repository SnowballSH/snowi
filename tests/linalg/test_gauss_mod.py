from itertools import permutations, product, takewhile

from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("linalg/gauss-mod")
gauss, det, solve_linear, mat_inv = (
  lib.gauss,
  lib.det,
  lib.solve_linear,
  lib.mat_inv,
)
snippet_globals = lib.gauss.__globals__


def set_mod(p):
  snippet_globals["MOD"] = p


def perm_det(a, p):
  n, res = len(a), 0
  for perm in permutations(range(n)):
    term = 1
    for i in range(n):
      term = term * a[i][perm[i]] % p
    inv = sum(perm[i] > perm[j] for i in range(n) for j in range(i + 1, n))
    res += -term if inv % 2 else term
  return res % p


def mul(a, x, p):
  return [sum(u * v for u, v in zip(row, x)) % p for row in a]


def random_matrix(n, m, p):
  kind = rnd(0, 3)
  if kind == 3:  # rank <= k
    k = rnd(0, min(n, m))
    u = [[rnd(0, p - 1) for _ in range(k)] for _ in range(n)]
    v = [[rnd(0, p - 1) for _ in range(m)] for _ in range(k)]
    return [
      [sum(u[i][t] * v[t][j] for t in range(k)) % p for j in range(m)]
      for i in range(n)
    ]
  hi = [min(2, p - 1), p - 1, 1][kind]
  a = [[rnd(0, hi) for _ in range(m)] for _ in range(n)]
  if kind == 2:
    a = [[p - 1 if x else rnd(0, 1) for x in row] for row in a]
  if n > 1 and rnd(0, 3) == 0:
    a[rnd(0, n - 1)] = a[rnd(0, n - 1)][:]
  return a


def check_rref(a, m, piv):
  r = len(piv)
  check(all(piv[i] < piv[i + 1] for i in range(r - 1)))
  for i, c in enumerate(piv):
    check(not any(a[i][:c]))
    check_eq([row[c] for row in a], [int(k == i) for k in range(len(a))])
  check(not any(any(row[:m]) for row in a[r:]))


def check_inverse(a, p):
  n = len(a)
  inv = mat_inv(a)
  check_eq(inv is not None, perm_det(a, p) != 0)
  if inv is not None:
    for j in range(n):
      col = [inv[i][j] for i in range(n)]
      check_eq(mul(a, col, p), [int(i == j) for i in range(n)])


# Small primes: count every solution by enumeration (Tweak: MOD).
for p in (2, 3, 5, 7):
  set_mod(p)
  for _ in range(iters(300)):
    n, m = rnd(0, 4), rnd(0, 4)
    a = random_matrix(n, m, p)
    b = [rnd(0, p - 1) for _ in range(n)]
    if rnd(0, 1):
      b = mul(a, [rnd(0, p - 1) for _ in range(m)], p)
    sols = [x for x in product(range(p), repeat=m) if mul(a, x, p) == b]
    kernel = [x for x in product(range(p), repeat=m) if not any(mul(a, x, p))]
    g = [row[:] for row in a]
    piv, d = gauss(g, m)
    check_rref(g, m, piv)
    check_eq(p ** (m - len(piv)), len(kernel))
    r, x = solve_linear(a, b, m)
    if sols:
      check_eq(p ** (m - r), len(sols))
      check(tuple(x) in sols)
    else:
      check_eq((r, x), (-1, None))
    if n == m:
      check_eq(d, perm_det(a, p))
      check_eq(det(a), d)
      check_inverse(a, p)
set_mod(998244353)
MOD = 998244353


# The real modulus: ranks from minors, determinants by permutations.
def minor_rank(a, m):
  n, best = len(a), 0
  for rows in range(1 << n):
    for cols in range(1 << m):
      k = bin(rows).count("1")
      if k != bin(cols).count("1") or k <= best:
        continue
      sub = [
        [a[i][j] for j in range(m) if cols >> j & 1]
        for i in range(n)
        if rows >> i & 1
      ]
      if perm_det(sub, MOD):
        best = k
  return best


for _ in range(iters(1500)):
  n, m = rnd(0, 4), rnd(0, 4)
  a = random_matrix(n, m, MOD)
  g = [row[:] for row in a]
  piv, d = gauss(g, m)
  check_eq(len(piv), minor_rank(a, m))
  check_rref(g, m, piv)
  if n == m:
    check_eq(d, perm_det(a, MOD))
    check_inverse(a, MOD)
  b = [rnd(0, 2) for _ in range(n)]
  if rnd(0, 1):
    b = mul(a, [rnd(0, MOD - 1) for _ in range(m)], MOD)
  ra = minor_rank(a, m)
  rab = minor_rank([row + [v] for row, v in zip(a, b)], m + 1)
  r, x = solve_linear(a, b, m)
  check_eq(r, ra if ra == rab else -1)
  if r >= 0:
    check_eq(mul(a, x, MOD), b)


def run_usage(**names):
  """Run the header's Usage lines in order, exactly as printed."""
  body = lib.__doc__.split("\nUsage:\n", 1)[1].split("\n")
  lines = takewhile(lambda line: line.startswith("  "), body)
  ns = dict(vars(lib), **names)
  exec("\n".join(line[2:] for line in lines), ns)
  return ns


# Usage top to bottom on one matrix: every call before the last
# one (gauss, in place) must still see the original matrix.
cases = [([[2, 1], [1, 3]], [3, 4])]
for _ in range(iters(300)):
  n = rnd(1, 4)
  a = random_matrix(n, n, MOD)
  b = [rnd(0, 2) for _ in range(n)]
  if rnd(0, 1):
    b = mul(a, [rnd(0, MOD - 1) for _ in range(n)], MOD)
  cases.append((a, b))
for a, b in cases:
  n = len(a)
  u = run_usage(a=[row[:] for row in a], b=b, m=n)
  ra = minor_rank(a, n)
  rab = minor_rank([row + [v] for row, v in zip(a, b)], n + 1)
  check_eq(u["r"], ra if ra == rab else -1)
  if u["r"] >= 0:
    check_eq(mul(a, u["x"], MOD), b)
  check_eq(u["d"], perm_det(a, MOD))
  check_eq(u["inv"] is None, u["d"] == 0)
  if u["inv"] is not None:
    for j in range(n):
      col = [u["inv"][i][j] for i in range(n)]
      check_eq(mul(a, col, MOD), [int(i == j) for i in range(n)])
  check_eq(len(u["piv"]), ra)
  check_rref(u["a"], n, u["piv"])

check_eq(solve_linear([], [], 3), (0, [0, 0, 0]))
check_eq((mat_inv([]), det([])), ([], 1))

n = 300
a = [[rnd(0, MOD - 1) for _ in range(n)] for _ in range(n)]
x0 = [rnd(0, MOD - 1) for _ in range(n)]
b = mul(a, x0, MOD)
with Timer("300 x 300 solve"):
  r, x = solve_linear(a, b, n)
check_eq((r, x), (n, x0))
a = [row[:200] for row in a[:200]]
with Timer("200 x 200 inverse"):
  inv = mat_inv(a)
for j in range(3):
  col = [inv[i][j] for i in range(200)]
  check_eq(mul(a, col, MOD), [int(i == j) for i in range(200)])
print("OK linalg/gauss-mod")
