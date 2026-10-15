from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("dp/dc-opt")
BIG = 1 << 62


def brute(old, cost):
  """new[i] = min over j < i of old[j] + cost(j, i), capped at BIG."""
  res = [BIG] * len(old)
  for i in range(1, len(old)):
    for j in range(i):
      res[i] = min(res[i], old[j] + cost(j, i))
  return res


def random_cost(n, c):
  """f(S_i - S_j) with f convex, plus weights of pairs inside [j, i),
  plus terms in i alone and in j alone: the quadrangle inequality."""
  quad = rnd(0, 2)
  w = [[rnd(0, 1) * rnd(0, c) for _ in range(n)] for _ in range(n)]
  pairs = [[0] * (n + 1) for _ in range(n + 1)]
  for j in range(n + 1):
    for i in range(j + 1, n + 1):
      pairs[j][i] = pairs[j][i - 1] + sum(w[x][i - 1] for x in range(j, i - 1))
  pre = [0]
  for _ in range(n):
    pre.append(pre[-1] + rnd(0, c))
  gi = [rnd(-c, c) for _ in range(n + 1)]
  hj = [rnd(-c, c) for _ in range(n + 1)]

  def cost(j, i):
    s = pre[i] - pre[j]
    return quad * s * s + pairs[j][i] + gi[i] + hj[j]

  return cost


def stress(rounds, maxn, c):
  for _ in range(iters(rounds)):
    n = rnd(0, maxn)
    cost = random_cost(n, c)
    dp = [BIG] * (n + 1)
    if rnd(0, 3):
      dp[0] = 0
    else:
      dp = [rnd(-c, c) if rnd(0, 1) else BIG for _ in range(n + 1)]
    for _ in range(rnd(1, n + 2)):
      want = brute(dp, cost)
      dp = lib.dc_layer(dp, cost)
      check_eq(dp, want)


def stress_smallest_opt(rounds, maxn):
  """Regression for the Precondition: no quadrangle inequality, only
  the smallest optimal j non-decreasing in i. Larger j often tie with
  it, so taking the largest optimal j (`<=`) must fail."""
  for _ in range(iters(rounds)):
    n, o = rnd(1, maxn), 0
    old = [rnd(-5, 5) for _ in range(n + 1)]
    w = [[0] * (n + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
      o = rnd(o, i - 1)  # the smallest optimal j for i
      best = rnd(-5, 5)
      for j in range(i):
        v = best + (rnd(1, 3) if j < o else rnd(0, 2) if j > o else 0)
        w[j][i] = v - old[j]
    cost = lambda j, i, w=w: w[j][i]
    check_eq(lib.dc_layer(old, cost), brute(old, cost), (old, w))


stress(3000, 6, 3)  # tiny values: many ties
stress(1500, 14, 10)
stress(200, 30, 10**6)
stress_smallest_opt(2000, 12)

check_eq(lib.dc_layer([5], lambda j, i: 0), [BIG])
dp = [0, BIG, BIG, BIG]
for _ in range(5):  # unreachable states stay capped
  dp = lib.dc_layer(dp, lambda j, i: (i - j) << 60)
check_eq(dp, [BIG] * 4)

# Large: one layer with squared group sums, n = 2e5.
n = 200_000
pre = [0]
for _ in range(n):
  pre.append(pre[-1] + rnd(0, 1000))
old = [rnd(0, 10**12) for _ in range(n + 1)]
calls = [0]


def sq(j, i):
  calls[0] += 1
  s = pre[i] - pre[j]
  return s * s


with Timer("2e5, one layer"):
  res = lib.dc_layer(old, sq)
check(calls[0] <= (n + 1) * 2 * 19, calls[0])
for _ in range(30):
  i = rnd(1, n)
  check_eq(res[i], min(old[j] + (pre[i] - pre[j]) ** 2 for j in range(i)))
print("OK dp/dc-opt")
