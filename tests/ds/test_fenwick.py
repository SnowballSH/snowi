from testlib import Timer, check_eq, iters, load, rnd

lib = load("ds/fenwick")


def brute_lower_bound(a, s):
  if s <= 0:
    return -1
  pre = 0
  for i, x in enumerate(a):
    pre += x
    if pre >= s:
      return i
  return len(a)


for _ in range(iters(2000)):
  n = rnd(1, 20)
  a = [0] * n
  fw = lib.Fenwick(n)
  for _ in range(30):
    i, x = rnd(0, n - 1), rnd(0, 10)
    a[i] += x
    fw.add(i, x)
    l = rnd(0, n)
    r = rnd(l, n)
    check_eq(fw.query(r) - fw.query(l), sum(a[l:r]))
    s = rnd(-2, sum(a) + 2)
    check_eq(fw.lower_bound(s), brute_lower_bound(a, s))

with Timer("2e5 adds + 2e5 queries"):
  n = 200_000
  fw = lib.Fenwick(n)
  for i in range(n):
    fw.add(i, 1)
  total = 0
  for i in range(n):
    total += fw.query(i)
  check_eq(total, n * (n - 1) // 2)
print("OK ds/fenwick")
