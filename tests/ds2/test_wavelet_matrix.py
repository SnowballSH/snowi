from itertools import accumulate

from testlib import Timer, check_eq, iters, load, rnd, rnd_list

lib = load("ds2/wavelet-matrix")


def brute_less(a, l, r, x):
  return sum(1 for v in a[l:r] if v < x)


def brute_count(a, l, r, lo, hi):  # none when lo >= hi
  return sum(1 for v in a[l:r] if lo <= v < hi)


class SumSmallest:  # the Tweak: sum of the k smallest in a[l:r]
  def __init__(self, a):
    self.wm = wm = lib.WaveletMatrix(a)
    self.ps = [None] * wm.lg  # prefix sums after level h's partition
    for h in range(wm.lg - 1, -1, -1):
      a = [x for x in a if not x >> h & 1] + [x for x in a if x >> h & 1]
      self.ps[h] = list(accumulate(a, initial=0))

  def sum(self, l, r, k):  # 0 <= k <= r - l
    wm, m = self.wm, self.wm.n + 1
    res = s = 0
    for h in range(wm.lg - 1, -1, -1):
      l1, r1 = wm.ones[h * m + l], wm.ones[h * m + r]
      z = r - l - (r1 - l1)
      if k < z:
        l, r = l - l1, r - r1
      else:
        k -= z
        res |= 1 << h
        s += self.ps[h][r - r1] - self.ps[h][l - l1]
        l, r = wm.zeros[h] + l1, wm.zeros[h] + r1
    return s + k * res  # the k copies of the k-th value left


def stress_sum(rounds, maxn, hi):
  for _ in range(iters(rounds)):
    n = rnd(0, maxn)
    a = rnd_list(n, 0, hi)
    ss = SumSmallest(a)
    for _ in range(25):
      l = rnd(0, n)
      r = rnd(l, n)
      k = rnd(0, r - l)
      check_eq(ss.sum(l, r, k), sum(sorted(a[l:r])[:k]))


def stress(rounds, maxn, hi):
  for _ in range(iters(rounds)):
    n = rnd(0, maxn)
    a = rnd_list(n, 0, hi)
    wm = lib.WaveletMatrix(a)
    for _ in range(25):
      l = rnd(0, n)
      r = rnd(l, n)
      if l < r:
        k = rnd(0, r - l - 1)
        check_eq(wm.kth(l, r, k), sorted(a[l:r])[k])
      x, y = rnd(-2, hi + 2), rnd(-2, hi + 2)
      if n and rnd(0, 3) == 0:
        x = a[rnd(0, n - 1)]
      check_eq(wm.count_less(l, r, x), brute_less(a, l, r, x))
      check_eq(wm.count(l, r, x, y), brute_count(a, l, r, x, y))
      check_eq(wm.rank(x, r), a[:r].count(x))


stress(1500, 10, 0)  # all zeros: no levels
stress(1500, 20, 3)  # many duplicates
stress(600, 70, 1000)
stress(200, 200, 2**62 - 1)  # 62 levels
stress_sum(1500, 20, 3)
stress_sum(300, 100, 2**62 - 1)
check_eq(SumSmallest([1, 1]).sum(0, 2, 1), 1)  # Tweak on duplicates
check_eq(lib.WaveletMatrix([5, 1]).count(0, 2, 5, 1), 0)  # lo > hi

wm = lib.WaveletMatrix([2**62 - 1, 0, 2**61])  # extremes for x, k
check_eq(wm.kth(0, 3, 2), 2**62 - 1)
check_eq(wm.count_less(0, 3, 2**63 - 1), 3)
check_eq(wm.count_less(0, 3, -(2**63)), 0)
check_eq((wm.count_less(0, 3, 2**61), wm.rank(0, 3)), (1, 1))


def answer(wm, queries):  # timed parts run in functions, as in main()
  total = 0
  for l, r, k, x in queries:
    total += wm.kth(l, r, k) + wm.count_less(l, r, x)
  return total


n = 200_000
a = rnd_list(n, 0, 10**9)
queries = []
for _ in range(n):
  l = rnd(0, n - 1)
  r = rnd(l + 1, n)
  queries.append((l, r, rnd(0, r - l - 1), rnd(0, 10**9)))
with Timer("build at n = 2e5, values < 1e9"):
  wm = lib.WaveletMatrix(a)
with Timer("2e5 kth + 2e5 count_less"):
  answer(wm, queries)
for l, r, k, x in queries[:30]:
  check_eq(wm.kth(l, r, k), sorted(a[l:r])[k])
  check_eq(wm.count_less(l, r, x), brute_less(a, l, r, x))
print("OK ds2/wavelet-matrix")
