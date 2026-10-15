from itertools import product

from testlib import Timer, check_eq, iters, load, rnd, rnd_list

lib = load("dp/wqs")
BIG = 1 << 62


def brute_split(a):
  """Exactly k non-empty groups, minimum sum of squared group sums."""
  n = len(a)
  F = [BIG] * (n + 1)
  if not n:
    F[0] = 0
  for mask in range(1 << max(n - 1, 0)):  # bit i: cut after a[i]
    total = s = 0
    for i in range(n):
      s += a[i]
      if i == n - 1 or mask >> i & 1:
        total += s * s
        s = 0
    if n:
      k = bin(mask).count("1") + 1
      F[k] = min(F[k], total)
  return F


def solve_split(a, lam):  # (value, count) tuples: ties to fewer
  pre = [0]
  for x in a:
    pre.append(pre[-1] + x)
  dp = [(0, 0)] + [(BIG, 0)] * len(a)
  for i in range(1, len(a) + 1):
    for j in range(i):
      s = pre[i] - pre[j]
      dp[i] = min(dp[i], (dp[j][0] + s * s + lam, dp[j][1] + 1))
  return dp[-1]


def brute_nonadj(a):
  """Exactly k pairwise non-adjacent elements, largest sum."""
  n = len(a)
  F = [-BIG] * (n + 1)
  for mask in range(1 << n):
    if not mask & mask >> 1:
      k = bin(mask).count("1")
      F[k] = max(F[k], sum(a[i] for i in range(n) if mask >> i & 1))
  return F


def solve_nonadj(a, lam):  # minimise -sum + lam * count
  skip, take = (0, 0), (BIG, 0)
  for x in a:
    t = (skip[0] - x + lam, skip[1] + 1)
    skip, take = min(skip, take), t
  return min(skip, take)


def brute_subarrays(a):
  """Exactly k disjoint non-empty subarrays, largest total."""
  n = len(a)
  F = [-BIG] * (n + 1)
  for labels in product(range(3), repeat=n):  # out, open, extend
    if all(t != 2 or (i and labels[i - 1]) for i, t in enumerate(labels)):
      k = labels.count(1)
      F[k] = max(F[k], sum(x for x, t in zip(a, labels) if t))
  return F


def solve_subarrays(a, lam):
  out, inside = (0, 0), (BIG, 0)
  for x in a:
    best = min(out, inside)
    opened = (best[0] - x + lam, best[1] + 1)
    out, inside = best, min(opened, (inside[0] - x, inside[1]))
  return min(out, inside)


for it in range(iters(800)):
  a = rnd_list(rnd(0, 10), 0, 1 if it % 2 else 5)
  F = brute_split(a)
  C = (sum(a) + 1) ** 2
  for k in range(1 if a else 0, len(a) + 1):
    check_eq(lib.wqs(k, -C, C, lambda lam, a=a: solve_split(a, lam)), F[k], a)

for it in range(iters(800)):
  a = rnd_list(rnd(0, 12), 0 if it % 2 else -3, 3)
  F = brute_nonadj(a)
  C = sum(map(abs, a))
  solve = lambda lam, a=a: solve_nonadj(a, lam)
  for k in range((len(a) + 1) // 2 + 1):
    check_eq(-lib.wqs(k, -C, C, solve), F[k], (a, k))
    check_eq(-lib.wqs(k, 0, C, solve), max(F[: k + 1]), (a, k))  # at most

for it in range(iters(200)):
  a = rnd_list(rnd(0, 8), -3, 0 if it % 2 else 3)
  F = brute_subarrays(a)
  C = sum(map(abs, a)) + 1
  for k in range(len(a) + 1):
    got = -lib.wqs(k, -C, C, lambda lam, a=a: solve_subarrays(a, lam))
    check_eq(got, F[k], (a, k))

check_eq(lib.wqs(0, -5, 5, lambda lam: (0, 0)), 0)  # no items at all


def solve_packed(a, lam):
  """Contest-speed oracle: (value, count) packed as value * B + count."""
  B = 1 << 20
  skip, take = 0, BIG
  for x in a:
    t = skip + (lam - x) * B + 1
    skip, take = (skip if skip < take else take), t
  best = skip if skip < take else take
  return best >> 20, best & (B - 1)


n = 2000  # moderate n against an O(n k) DP
a = rnd_list(n, 0, 10**6)
m = n // 2 + 2
no, yes = [0] + [-BIG] * m, [-BIG] * (m + 1)  # last skipped / taken
for x in a:
  nno = [max(no[c], yes[c]) for c in range(m + 1)]
  nyes = [-BIG] + [no[c] + x if no[c] > -BIG else -BIG for c in range(m)]
  no, yes = nno, nyes
C = sum(a)
for _ in range(10):
  k = rnd(0, (n + 1) // 2)
  check_eq(
    -lib.wqs(k, -C, C, lambda lam: solve_packed(a, lam)), max(no[k], yes[k])
  )

# Regression: a = 0, M, 0, ..., M, 0. The last slope is the whole
# sum (k = m + 1 forces every 0), so lam * count passes 2^63.
m, M = 100_000, 10**9
alt = [0, M] * m + [0]
with Timer("2e5 alternating, big ints, 3 searches"):
  for k, want in ((1, M), (m, m * M), (m + 1, 0)):
    got = -lib.wqs(k, -sum(alt), sum(alt), lambda lam: solve_packed(alt, lam))
    check_eq(got, want, k)

n = 200_000
same = [7] * n  # straight F: every count ties at lam = 7
with Timer("2e5 elements, 5 searches"):
  for k in (0, 1, 2, n // 4, n // 2):
    got = -lib.wqs(k, -7 * n, 7 * n, lambda lam: solve_packed(same, lam))
    check_eq(got, 7 * k)
rand = rnd_list(n, 0, 10**9)
with Timer("2e5 random elements, 1 search"):
  lib.wqs(n // 3, -sum(rand), sum(rand), lambda lam: solve_packed(rand, lam))
print("OK dp/wqs")
