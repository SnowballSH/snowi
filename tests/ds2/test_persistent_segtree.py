from bisect import bisect_left

from testlib import Timer, check, check_eq, iters, load, rnd, rnd_list

lib = load("ds2/persistent-segtree")


def brute_kth(d, k):  # smallest i whose count over [0, i] exceeds k
  pre = 0
  for i, x in enumerate(d):
    pre += x
    if pre > k:
      return i
  return -1


def stress_versions(rounds, c):
  for _ in range(iters(rounds)):
    n = rnd(1, 12)
    st = lib.PersistentSegTree(n)
    root, arr = [0], [[0] * n]
    for _ in range(40):
      v, kind = rnd(0, len(root) - 1), rnd(0, 2)
      if kind == 0:
        i, d = rnd(0, n - 1), rnd(-c, c)
        root.append(st.add(root[v], i, d))
        arr.append(arr[v][:])
        arr[-1][i] += d
      elif kind == 1:
        l = rnd(0, n)
        r = rnd(l, n)
        check_eq(st.query(root[v], l, r), sum(arr[v][l:r]))
      else:  # b = version v plus some positive adds
        b, d = root[v], [0] * n
        for _ in range(rnd(1, 4)):
          i, x = rnd(0, n - 1), rnd(1, 3)
          b = st.add(b, i, x)
          d[i] += x
        k = rnd(0, sum(d) - 1)
        check_eq(st.kth(root[v], b, k), brute_kth(d, k))
    for v in range(len(root)):
      for i in range(n):
        check_eq(st.query(root[v], i, i + 1), arr[v][i])
    check_eq(st.query(0, 0, n), 0)


class RangeKth:  # the Example: k-th smallest, count < y in a[l:r]
  def __init__(self, a):
    self.vals = vals = sorted(set(a))
    self.st = st = lib.PersistentSegTree(len(vals))
    self.root = root = [0]
    for x in a:
      root.append(st.add(root[-1], bisect_left(vals, x), 1))

  def kth(self, l, r, k):
    return self.vals[self.st.kth(self.root[l], self.root[r], k)]

  def less(self, l, r, y):
    c, st = bisect_left(self.vals, y), self.st
    return st.query(self.root[r], 0, c) - st.query(self.root[l], 0, c)


def stress_range_kth(rounds, lo, hi):
  for _ in range(iters(rounds)):
    n = rnd(1, 15)
    a = rnd_list(n, lo, hi)
    rk = RangeKth(a)
    for _ in range(20):
      l = rnd(0, n - 1)
      r = rnd(l + 1, n)
      b = sorted(a[l:r])
      k = rnd(0, r - l - 1)
      check_eq(rk.kth(l, r, k), b[k])
      y = rnd(lo - 1, hi + 1)
      check_eq(rk.less(l, r, y), bisect_left(b, y))


stress_versions(1500, 5)
stress_versions(300, 10**17)
stress_range_kth(1500, 0, 5)  # many duplicates
stress_range_kth(500, -(10**18), 10**18)

st = lib.PersistentSegTree(6)  # the Usage, i = 4, j = 1
v1 = st.add(0, 4, 5)
v2 = st.add(v1, 1, 3)
d = [st.query(v2, p, p + 1) - st.query(v1, p, p + 1) for p in range(6)]
check(min(d) >= 0)  # kth needs v2 - v1 >= 0
for k in range(3):
  check_eq(st.kth(v1, v2, k), brute_kth(d, k))
check_eq(st.query(v1, 0, 6), 5)

st = lib.PersistentSegTree(1)  # n = 1: a single leaf
v = st.add(st.add(0, 0, 2), 0, 3)
check_eq((st.query(v, 0, 1), st.kth(0, v, 4)), (5, 0))


def build(a):  # timed parts run in functions, as in main()
  return RangeKth(a)


def answer(rk, queries):
  total = 0
  for l, r, k in queries:
    total += rk.kth(l, r, k) + rk.less(l, r, k)
  return total


n = 200_000
a = rnd_list(n, 0, 10**9)
queries = []
for _ in range(n):
  l = rnd(0, n - 1)
  r = rnd(l + 1, n)
  queries.append((l, r, rnd(0, r - l - 1)))
with Timer("2e5 adds (prefix versions) at n = 2e5"):
  rk = build(a)
check(len(rk.st.S) <= 1 + n * 19)
with Timer("2e5 kth + 4e5 prefix sums"):
  answer(rk, queries)
for l, r, k in queries[:20]:
  check_eq(rk.kth(l, r, k), sorted(a[l:r])[k])
print("OK ds2/persistent-segtree")
