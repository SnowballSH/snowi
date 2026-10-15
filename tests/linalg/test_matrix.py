import random
from math import gcd

from testlib import CONTENT, Timer, check_eq, iters, load, rnd

lib = load("linalg/matrix")
mat_mul, mat_pow, MOD = lib.mat_mul, lib.mat_pow, lib.MOD


def naive_mul(a, b):
  rows, inner = len(a), len(b)
  cols = len(b[0]) if b else 0
  return [
    [sum(a[i][k] * b[k][j] for k in range(inner)) % MOD for j in range(cols)]
    for i in range(rows)
  ]


def naive_pow(a, e):
  n = len(a)
  r = [[int(i == j) for j in range(n)] for i in range(n)]
  for _ in range(e):
    r = naive_mul(r, a)
  return r


def count_walks(adj, s, t, length):
  if length == 0:
    return int(s == t)
  return sum(
    adj[s][v] * count_walks(adj, v, t, length - 1)
    for v in range(len(adj))
    if adj[s][v]
  )


def rand_matrix(n, m, hi):
  return [[rnd(0, hi) if rnd(0, 2) else 0 for _ in range(m)] for _ in range(n)]


for _ in range(iters(1000)):
  n = rnd(0, 6)
  hi = 3 if rnd(0, 1) else MOD - 1
  a = rand_matrix(n, n, hi)
  e = rnd(0, 25)
  check_eq(mat_pow(a, e), naive_pow(a, e))
  m, p = rnd(0, 5), rnd(0, 5)  # rectangular products
  b, c = rand_matrix(n, m, hi), rand_matrix(m, p, MOD - 1)
  check_eq(mat_mul(b, c), naive_mul(b, c))
  e1, e2 = rnd(0, 10**18), rnd(0, 10**18)
  check_eq(mat_pow(a, e1 + e2), naive_mul(mat_pow(a, e1), mat_pow(a, e2)))

for _ in range(iters(300)):  # walk counting
  n, length = rnd(1, 5), rnd(0, 7)
  adj = rand_matrix(n, n, 2)
  p = mat_pow(adj, length)
  for s in range(n):
    for t in range(n):
      check_eq(p[s][t], count_walks(adj, s, t, length))

for _ in range(iters(300)):  # huge e: permutation, shear, diagonal
  n = rnd(1, 8)
  perm = list(range(n))
  random.Random(rnd(0, 10**9)).shuffle(perm)
  e = rnd(0, 10**18)
  order = 1
  for i in range(n):
    length, j = 1, perm[i]
    while j != i:
      length, j = length + 1, perm[j]
    order = order * length // gcd(order, length)
  q = list(range(n))
  for _ in range(e % order):
    q = [perm[v] for v in q]
  pm = [[int(perm[i] == j) for j in range(n)] for i in range(n)]
  check_eq(
    mat_pow(pm, e), [[int(q[i] == j) for j in range(n)] for i in range(n)]
  )
  check_eq(mat_pow([[1, 1], [0, 1]], e), [[1, e % MOD], [0, 1]])
  d = [rnd(0, MOD - 1) for _ in range(n)]
  dg = [[d[i] if i == j else 0 for j in range(n)] for i in range(n)]
  want = [
    [pow(d[i], e, MOD) if i == j else 0 for j in range(n)] for i in range(n)
  ]
  check_eq(mat_pow(dg, e), want)

TINF = 1 << 62


def tropical_pow(best):
  """mat_pow with the Tweak applied literally: (min, +) for
  best = "min", (max, +) with -INF in place of INF for "max"."""
  code = (CONTENT / "linalg/matrix.py").read_text().split('"""', 2)[2]
  inf = "INF" if best == "min" else "-INF"
  edits = [
    ("acc = [0] * m", f"acc = [{inf}] * m"),
    ("if x:", f"if x != {inf}:"),
    (
      "acc[j] = (acc[j] + x * brow[j]) % mod",
      f"acc[j] = {best}(acc[j], x + brow[j])",
    ),
    (
      "[int(i == j) for j in range(n)]",
      f"[0 if i == j else {inf} for j in range(n)]",
    ),
  ]
  for old, new in edits:
    check_eq(code.count(old), 1, old)
    code = code.replace(old, new)
  ns = {"MOD": MOD, "INF": TINF}
  exec(code, ns)
  return ns["mat_pow"]


def best_walks(w, e, best):
  """Optimal walks of exactly e edges by a DP over the length;
  None: no edge or no walk."""
  n = len(w)
  d = [[0 if i == j else None for j in range(n)] for i in range(n)]
  for _ in range(e):
    nd = [[None] * n for _ in range(n)]
    for s in range(n):
      for u in range(n):
        for v in range(n):
          if d[s][u] is not None and w[u][v] is not None:
            c = d[s][u] + w[u][v]
            nd[s][v] = c if nd[s][v] is None else best(nd[s][v], c)
    d = nd
  return d


for best, sign in ((min, 1), (max, -1)):  # Tweak: (min, +), (max, +)
  tpow = tropical_pow(best.__name__)
  for _ in range(iters(400)):
    n, e, hi = rnd(1, 5), rnd(0, 8), 3 if rnd(0, 1) else 10**9
    w = [
      [rnd(-hi, hi) if rnd(0, 2) else None for _ in range(n)] for _ in range(n)
    ]
    m = [[sign * TINF if v is None else v for v in row] for row in w]
    got = [
      [None if sign * v >= TINF // 2 else v for v in row] for row in tpow(m, e)
    ]
    check_eq(got, best_walks(w, e, best))

check_eq(mat_pow([], 5), [])
check_eq(mat_pow([[0] * 3 for _ in range(3)], 0), naive_pow([[0] * 3] * 3, 0))
check_eq(mat_mul([[1, 2]], [[3], [4]]), [[11]])


def speed(n):
  m = [[rnd(0, MOD - 1) for _ in range(n)] for _ in range(n)]
  with Timer(f"{n} x {n} to the power 1e18"):
    big = mat_pow(m, 10**18)
  check_eq(mat_mul(big, m), mat_mul(m, big))


speed(64)
print("OK linalg/matrix")
