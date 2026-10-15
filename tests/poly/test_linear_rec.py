from testlib import Timer, check, check_eq, iters, load, rnd, rnd_list

lib = load("poly/linear-rec")
MOD = lib.MOD


def fits(s, L):
  """Does some recurrence of order L fit s? Gaussian elimination."""
  rows = [[s[i - 1 - j] for j in range(L)] + [s[i]] for i in range(L, len(s))]
  r = 0
  for col in range(L):
    p = next((i for i in range(r, len(rows)) if rows[i][col]), None)
    if p is None:
      continue
    rows[r], rows[p] = rows[p], rows[r]
    inv = pow(rows[r][col], MOD - 2, MOD)
    rows[r] = [x * inv % MOD for x in rows[r]]
    for i, row in enumerate(rows):
      if i != r and row[col]:
        f = row[col]
        rows[i] = [(x - f * y) % MOD for x, y in zip(row, rows[r])]
    r += 1
  return all(row[L] == 0 for row in rows[r:])


def min_order(s):
  L = 0
  while not fits(s, L):
    L += 1
  return L


def extend(s, c, n):
  s = s[:]
  while len(s) < n:
    x = 0
    for j, y in enumerate(c):
      x = (x + y * s[-1 - j]) % MOD
    s.append(x)
  return s


def generates(s, c):
  return extend(s[: len(c)], c, len(s)) == s if len(c) <= len(s) else False


def matrix_kth(s, c, k):
  """s_k from the companion matrix power, independent of linear_rec."""
  L = len(c)
  if k < L:
    return s[k]

  def mul(A, B):
    return [
      [sum(A[i][t] * B[t][j] for t in range(L)) % MOD for j in range(L)]
      for i in range(L)
    ]

  M = [c[:]] + [[int(j == i - 1) for j in range(L)] for i in range(1, L)]
  P = [[int(i == j) for j in range(L)] for i in range(L)]
  e = k - L + 1
  while e:
    if e & 1:
      P = mul(P, M)
    M = mul(M, M)
    e >>= 1
  return sum(P[0][j] * s[L - 1 - j] for j in range(L)) % MOD


for it in range(iters(800)):
  n = rnd(0, 12)
  s = rnd_list(n, 0, 2 if it % 2 else MOD - 1)
  c = lib.berlekamp_massey(s)
  check_eq(len(c), min_order(s), f"s={s}")
  check(generates(s, c))

for it in range(iters(300)):
  L = rnd(0, 8)
  c = rnd_list(L, 0, MOD - 1 if it % 3 else 2)
  s = extend(rnd_list(L, 0, MOD - 1), c, 2 * L + rnd(0, 4))
  got = lib.berlekamp_massey(s)
  check(len(got) <= L and generates(s, got))
  full = extend(s, c, 300)
  check(generates(full, got))
  for t in range(20):
    k = t if t < 10 else rnd(0, 299)
    check_eq(lib.linear_rec(s, got, k), full[k])
  k = rnd(0, 10**18)
  check_eq(lib.linear_rec(s, got, k), matrix_kth(s, got, k))
  check_eq(lib.linear_rec(s, c, k), matrix_kth(s, c, k))
  # Hand-written input may use negative representatives.
  sn = [x - MOD * rnd(0, 1) for x in s]
  cn = [x - MOD * rnd(0, 1) for x in c]
  check_eq(lib.berlekamp_massey(sn), got)
  for _ in range(5):
    kk = rnd(0, 299)
    check_eq(lib.linear_rec(sn, cn, kk), full[kk])
  check_eq(lib.linear_rec(sn, cn, k), matrix_kth(s, c, k))

check_eq(lib.berlekamp_massey([]), [])
check_eq(lib.berlekamp_massey([0, 0, 0]), [])
check_eq(len(lib.berlekamp_massey([0, 0, 1])), 3)
check_eq(lib.berlekamp_massey([1, 1, 2, 3, 5, 8]), [1, 1])
pw = [(pow(2, k, MOD) + pow(3, k, MOD)) % MOD for k in range(10)]
check_eq(lib.berlekamp_massey(pw), [5, MOD - 6])  # orders add
check_eq(lib.linear_rec([], [], 10**18), 0)
check_eq(lib.linear_rec([7], [0], 0), 7)
check_eq(lib.linear_rec([7], [0], 1), 0)
check_eq(lib.linear_rec([0, 1], [1, 1], 90), 2880067194370816120 % MOD)
check_eq(lib.linear_rec([0, 1], [2, -1], 10**18), 10**18 % MOD)
check_eq(lib.linear_rec([1, 1], [-1, 0], 10), MOD - 1)
check_eq(lib.linear_rec([-3], [1], 0), MOD - 3)

with Timer("berlekamp_massey, 10^4 terms of order 4995"):
  L = 4995
  c = rnd_list(L, 0, MOD - 1)
  s = extend(rnd_list(L, 0, MOD - 1), c, 10**4)
  check(lib.berlekamp_massey(s) == c)
with Timer("linear_rec, order 500, k = 10^18"):
  L = 500
  c, s = c[:L], s[:L]
  period = [0] * (L - 1) + [1]  # s_i = s_(i - L)
  check_eq(lib.linear_rec(s, period, 10**18), s[10**18 % L])
check_eq(lib.linear_rec(s, c, 2999), extend(s, c, 3000)[2999])
print("OK poly/linear-rec")
