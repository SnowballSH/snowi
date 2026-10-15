from testlib import Timer, check_eq, iters, load, rnd

lib = load("dp/knuth")


def brute(n, w):
  """O(n^3) interval DP, dp[i][i + 1] = 0."""
  dp = [[0] * (n + 1) for _ in range(n + 1)]
  for length in range(2, n + 1):
    for l in range(n - length + 1):
      r = l + length
      best = min(dp[l][k] + dp[k][r] for k in range(l + 1, r))
      dp[l][r] = best + w(l, r)
  return dp


def random_w(n, c):
  """f(sum) with f convex increasing, plus pair weights: the
  quadrangle inequality and monotone under inclusion."""
  lin, quad = rnd(0, 3), rnd(0, 2)
  p = [[0 if rnd(0, 2) else rnd(0, c) for _ in range(n)] for _ in range(n)]
  pairs = [[0] * (n + 1) for _ in range(n + 1)]
  for l in range(n + 1):
    for r in range(l + 1, n + 1):
      pairs[l][r] = pairs[l][r - 1] + sum(p[x][r - 1] for x in range(l, r - 1))
  pre = [0]
  for _ in range(n):
    pre.append(pre[-1] + rnd(0, c))

  def w(l, r):
    s = pre[r] - pre[l]
    return lin * s + quad * s * s + pairs[l][r]

  return w


def stress(rounds, maxn, c):
  for _ in range(iters(rounds)):
    n = rnd(0, maxn)
    w = random_w(n, c)
    check_eq(lib.knuth(n, w), brute(n, w))


stress(3000, 7, 3)  # tiny values: many ties
stress(1000, 14, 20)
stress(100, 40, 10**6)
check_eq(lib.knuth(0, lambda l, r: 1), [[0]])
check_eq(lib.knuth(1, lambda l, r: 1)[0][1], 0)
check_eq(lib.knuth(2, lambda l, r: 7)[0][2], 7)

for _ in range(iters(300)):  # optimal BST via n + 1 boundaries
  n = rnd(0, 9)
  pre = [0]
  for _ in range(n):
    pre.append(pre[-1] + rnd(0, 10))
  e = [[0] * (n + 1) for _ in range(n + 1)]  # keys [l, r), roots k
  for length in range(1, n + 1):
    for l in range(n - length + 1):
      r = l + length
      e[l][r] = min(e[l][k] + e[k + 1][r] for k in range(l, r))
      e[l][r] += pre[r] - pre[l]
  dp = lib.knuth(n + 1, lambda l, r, pre=pre: pre[r - 1] - pre[l])
  check_eq(dp[0][n + 1], e[0][n])


def product_w(n):
  """For large n: f(sum) plus pairs x < y inside with weight
  p[x] * q[y], in O(1) from prefix sums."""
  pre, P, Q, QP = [0], [0], [0], [0]
  for _ in range(n):
    p, q = rnd(0, 100), rnd(0, 100)
    pre.append(pre[-1] + rnd(0, 1000))
    QP.append(QP[-1] + q * P[-1])
    P.append(P[-1] + p)
    Q.append(Q[-1] + q)

  def w(l, r):
    s = pre[r] - pre[l]
    return s * s + (QP[r] - QP[l]) - P[l] * (Q[r] - Q[l])

  return w


def units(length):
  """Cost of merging `length` unit piles: a balanced tree."""
  q = length.bit_length() - 1
  return length * q + 2 * (length - (1 << q))


n = 2000
with Timer("n = 2000, unit piles"):
  dp = lib.knuth(n, lambda l, r: r - l)
for l in range(n):
  row = dp[l]
  for r in range(l + 1, n + 1):
    if row[r] != units(r - l):
      check_eq(row[r], units(r - l), (l, r))
w = product_w(n)
with Timer("n = 2000, random w"):
  dp = lib.knuth(n, w)
small = brute(80, w)
for l in range(81):
  check_eq(dp[l][l + 1 : 81], small[l][l + 1 : 81])
print("OK dp/knuth")
