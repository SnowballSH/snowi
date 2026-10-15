from testlib import Timer, check_eq, iters, load, rnd

lib = load("ds2/treap")


def stress(rounds, c):
  for _ in range(iters(rounds)):
    tr, a = lib.Treap(), []
    for _ in range(60):
      n, kind = len(a), rnd(0, 6)
      l, r = sorted((rnd(0, n), rnd(0, n)))
      if kind <= 1 or n == 0:
        v = rnd(-c, c)
        tr.insert(l, v)
        a.insert(l, v)
      elif kind == 2:
        i = rnd(0, n - 1)
        tr.erase(i)
        del a[i]
      elif kind == 3:
        tr.reverse(l, r)
        a[l:r] = a[l:r][::-1]
      elif kind == 4:
        check_eq(tr.query(l, r), sum(a[l:r]))
        check_eq(tr.prefix(r), sum(a[:r]))
      elif kind == 5:
        i = rnd(0, n - 1)
        check_eq(tr.get(i), a[i])
      else:  # Usage: move the first k elements to the end
        x, y = tr.split(tr.root, l)
        tr.root = tr.merge(y, x)
        a = a[l:] + a[:l]
      check_eq(tr.size(), len(a))
    check_eq(tr.to_list(), a)


def cut_paste(rounds):  # cut [l, r) and paste at p of the rest
  for _ in range(iters(rounds)):
    n = rnd(0, 20)
    tr, a = lib.Treap(), [rnd(0, 9) for _ in range(n)]
    for i, v in enumerate(a):
      tr.insert(i, v)
    l = rnd(0, n)
    r = rnd(l, n)
    p = rnd(0, n - (r - l))
    x0, bc = tr.split(tr.root, l)
    b, c = tr.split(bc, r - l)
    x, y = tr.split(tr.merge(x0, c), p)
    tr.root = tr.merge(tr.merge(x, b), y)
    rest = a[:l] + a[r:]
    check_eq(tr.to_list(), rest[:p] + a[l:r] + rest[p:])


stress(1500, 9)
stress(500, 10**18)
cut_paste(1000)

tr = lib.Treap()  # empty ranges and a single element
check_eq((tr.query(0, 0), tr.size(), tr.to_list()), (0, 0, []))
tr.reverse(0, 0)
tr.insert(0, 5)
tr.reverse(0, 1)
tr.reverse(1, 1)
check_eq((tr.get(0), tr.query(0, 1), tr.to_list()), (5, 5, [5]))
tr.erase(0)
check_eq((tr.size(), tr.to_list()), (0, []))


def inserts(tr, pos):  # timed parts run in functions, as in main()
  for i, p in enumerate(pos):
    tr.insert(p % (i + 1), i)


def reversals(tr, ranges):
  for l, r in ranges:
    tr.reverse(l, r)


def reads(tr, ranges):
  total = 0
  for l, r in ranges:
    total += tr.query(l, r) + tr.get(l)
  return total


n = 100_000
pos = [rnd(0, n) for _ in range(n)]
ranges = []
for _ in range(2 * n):
  l = rnd(0, n - 1)
  ranges.append((l, rnd(l + 1, n)))
tr = lib.Treap()
with Timer("1e5 inserts at random positions"):
  inserts(tr, pos)
with Timer("1e5 reversals at n = 1e5"):
  reversals(tr, ranges[:n])
with Timer("2e5 range sums + 2e5 gets at n = 1e5"):
  reads(tr, ranges)
for l, r in ranges[:2000]:  # a reversal keeps the sum, swaps ends
  s, first, last = tr.query(l, r), tr.get(l), tr.get(r - 1)
  tr.reverse(l, r)
  check_eq((tr.query(l, r), tr.get(l), tr.get(r - 1)), (s, last, first))
check_eq(sorted(tr.to_list()), list(range(n)))
check_eq(tr.query(0, n), n * (n - 1) // 2)


def mixed(tr, ops):
  total = 0
  for t, i, j in ops:
    if t == 0:
      tr.insert(i % (tr.size() + 1), j)
    elif t == 1:
      tr.erase(i % tr.size())
    elif t == 2:
      l, r = sorted((i % tr.size(), j % tr.size()))
      tr.reverse(l, r + 1)
    else:
      l, r = sorted((i % tr.size(), j % tr.size()))
      total += tr.query(l, r + 1)
  return total


ops = [(rnd(0, 3), rnd(0, 10**9), rnd(0, 10**9)) for _ in range(2 * n)]
with Timer("2e5 mixed ops at n = 1e5 (insert, erase, reverse, sum)"):
  mixed(tr, ops)
print("OK ds2/treap")
