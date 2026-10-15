from testlib import Timer, check_eq, iters, load, rnd

lib = load("ds2/binary-trie")


def check_all(tr, ys, x, bits):
  full = 1 << bits
  n = len(ys)
  check_eq(tr.size(), n)
  xs = sorted(y ^ x for y in ys)
  srt = sorted(ys)
  for k in range(n):
    check_eq(tr.kth(k, x), xs[k])
    check_eq(tr.kth(k), srt[k])
  if n:
    check_eq((tr.min_xor(x), tr.max_xor(x)), (xs[0], xs[-1]))
  for k in (0, 1, x, full - 1, full, full + 5, rnd(0, full)):
    check_eq(tr.count_less(x, k), sum(v < k for v in xs))
  check_eq(tr.count_less(x, 1), ys.count(x))


def stress(rounds, bits, hi):
  for _ in range(iters(rounds)):
    tr, ys = lib.BinaryTrie(bits), []
    for _ in range(30):
      kind = rnd(0, 3)
      if kind <= 1 or not ys:
        y, c = rnd(0, hi), 1 if kind == 0 else rnd(1, 3)
        tr.add(y, c)
        ys += [y] * c
      else:
        y = ys.pop(rnd(0, len(ys) - 1))
        tr.add(y, -1)
      check_all(tr, ys, rnd(0, hi), bits)
    for y in ys:  # empty again, nodes stay
      tr.add(y, -1)
    check_all(tr, [], rnd(0, hi), bits)
    tr.add(hi)
    check_all(tr, [hi], rnd(0, hi), bits)


stress(1000, 4, 15)
stress(600, 30, 2**30 - 1)
stress(400, 60, 2**60 - 1)
stress(200, 1, 1)

tr = lib.BinaryTrie(30)  # extremes
tr.add(0)
tr.add(2**30 - 1)
check_eq((tr.max_xor(0), tr.min_xor(2**30 - 1)), (2**30 - 1, 0))
check_eq((tr.count_less(0, 2**30), tr.count_less(0, 2**30 - 1)), (2, 1))


def inserts(tr, ys):  # timed parts run in functions, as in main()
  for y in ys:
    tr.add(y)


def queries(tr, xs):
  total = 0
  for x in xs:
    total += tr.max_xor(x) + tr.count_less(x, x)
  return total


n = 200_000
ys = [rnd(0, 2**30 - 1) for _ in range(n)]
xs = [rnd(0, 2**30 - 1) for _ in range(n)]
tr = lib.BinaryTrie(30)
with Timer("2e5 inserts, B = 30"):
  inserts(tr, ys)
with Timer("2e5 max_xor + 2e5 count_less, B = 30"):
  queries(tr, xs)
for x in xs[:20]:
  check_eq(tr.max_xor(x), max(x ^ y for y in ys))
  check_eq(tr.count_less(x, x), sum((x ^ y) < x for y in ys))
srt = sorted(ys)
for k in range(0, n, 1009):
  check_eq(tr.kth(k), srt[k])
print("OK ds2/binary-trie")
