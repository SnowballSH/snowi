from bisect import bisect_left

from testlib import Timer, check, check_eq, iters, load, rnd, rnd_list

lib = load("dp/lis")


def brute(a):
  """Length of the longest strictly increasing subsequence, O(n^2)."""
  best, length = 0, [1] * len(a)
  for i in range(len(a)):
    for j in range(i):
      if a[j] < a[i]:
        length[i] = max(length[i], length[j] + 1)
    best = max(best, length[i])
  return best


def check_valid(a, idx):
  for t, i in enumerate(idx):
    check(0 <= i < len(a))
    if t:
      check(idx[t - 1] < i and a[idx[t - 1]] < a[i], (a, idx))


def fenwick_length(a):
  """Independent O(n log n) length: prefix max over value ranks."""
  vals = sorted(set(a))
  t, best = [0] * (len(vals) + 1), 0
  for x in a:
    r, cur = bisect_left(vals, x), 0
    i = r
    while i > 0:
      cur = max(cur, t[i])
      i &= i - 1
    cur += 1
    i = r + 1
    while i <= len(vals):
      t[i] = max(t[i], cur)
      i += i & -i
    best = max(best, cur)
  return best


for it in range(iters(15000)):
  n = rnd(0, 12)
  c = (2, 8, 10**12)[it % 3]
  a = rnd_list(n, -c, c)
  idx = lib.lis(a)
  check_valid(a, idx)
  check_eq(len(idx), brute(a), a)

for _ in range(iters(1000)):  # tuples compare lexicographically
  a = [(rnd(0, 2), rnd(0, 2)) for _ in range(rnd(0, 10))]
  idx = lib.lis(a)
  check_valid(a, idx)
  check_eq(len(idx), brute(a), a)

check_eq(lib.lis([]), [])
check_eq(lib.lis([7]), [0])
check_eq(lib.lis([5, 5, 5]), [2])
check_eq(lib.lis([3, 2, 1]), [2])
check_eq(lib.lis([-(2**63), 2**63 - 1]), [0, 1])
check_eq(lib.lis([3, 1, 4, 1, 5, 9, 2, 6]), [1, 2, 4, 7])

with Timer("2e5 random values"):
  a = rnd_list(200_000, 0, 10**6)
  idx = lib.lis(a)
check_valid(a, idx)
check_eq(len(idx), fenwick_length(a))
check_eq(len(lib.lis(list(range(200_000)))), 200_000)
print("OK dp/lis")
