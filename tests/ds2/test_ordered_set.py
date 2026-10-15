from bisect import bisect_left, bisect_right, insort

from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("ds2/ordered-set")
SortedList = lib.SortedList


def check_all(s, ref, values):
  check_eq(len(s), len(ref))
  check_eq(list(s), ref)
  for x in values:
    check_eq(x in s, x in ref)
    check_eq(s.count(x), ref.count(x))
    check_eq(s.bisect_left(x), bisect_left(ref, x))
    check_eq(s.bisect_right(x), bisect_right(ref, x))
  for k in range(-len(ref), len(ref)):
    check_eq(s[k], ref[k])


def stress(rounds, gen, ratio):
  SortedList.RATIO = ratio  # tiny ratios force many splits
  for _ in range(iters(rounds)):
    init = [gen() for _ in range(rnd(0, 12))]
    s, ref = SortedList(init), sorted(init)
    for _ in range(60):
      op, x = rnd(0, 9), gen()
      if op <= 3:
        s.add(x)
        insort(ref, x)
      elif op <= 6:
        if ref and rnd(0, 1):
          x = rng.choice(ref)
        check_eq(s.discard(x), x in ref)
        if x in ref:
          ref.remove(x)
      elif op == 7 and ref:
        k = rnd(-len(ref), len(ref) - 1)
        check_eq(s.pop(k), ref.pop(k))
      elif op == 8:
        try:
          s.remove(x)
          check(x in ref, "remove of an absent value")
          ref.remove(x)
        except ValueError:
          check(x not in ref, "remove raised for a present value")
      elif op == 9 and ref:
        k = rnd(-len(ref), len(ref) - 1)
        check_eq(s[k], ref[k])
      check_eq(len(s), len(ref))
    check_all(s, ref, [gen() for _ in range(5)])
    while ref:
      check_eq(s.pop(0), ref.pop(0))
    check_eq(len(s), 0)
    check_all(s, ref, [gen()])
  SortedList.RATIO = 32


stress(1000, lambda: rnd(0, 8), 1)
stress(800, lambda: rnd(-50, 50), 2)
stress(300, lambda: rnd(-50, 50), 32)
stress(200, lambda: (rnd(0, 3), rng.choice("ab")), 1)

# Edge cases: empty, a single value, extremes, out-of-range k.
s = SortedList()
check_eq((len(s), list(s), 5 in s, s.discard(5)), (0, [], False, False))
check_eq((s.bisect_left(5), s.bisect_right(5), s.count(5)), (0, 0, 0))
for bad in (0, -1):
  try:
    s[bad]
    check(False, "index into an empty list")
  except IndexError:
    pass
big = 10**18
s = SortedList([big, -big, 0, big])
check_eq(list(s), [-big, 0, big, big])
check_eq((s.count(big), s[-1], s.bisect_right(big)), (2, big, 4))


def mixed(s, ops, n):  # in a function, like contest code in main()
  total = 0
  for t, x in ops:
    if t == 0:
      s.add(x)
    elif t == 1:
      s.discard(s[x % len(s)])
    elif t == 2:
      total += s.bisect_left(x)
    elif t == 3:
      total += s[x % len(s)]
    else:
      total += x in s
  return total


def churn(s, ops):  # add and discard keep n fixed
  total = 0
  for t, x in ops:
    if t <= 1:
      s.add(x)
      s.discard(x)
    elif t == 2:
      total += s.bisect_left(x)
    elif t == 3:
      total += s[x % len(s)]
    else:
      total += x in s
  return total


n = 200_000
values = [rnd(0, 10**9) for _ in range(n)]
ops = [(rnd(0, 4), rnd(0, 10**9)) for _ in range(2 * n)]
with Timer("2e5 adds then 2e5 mixed ops, from empty"):
  s = SortedList()
  for x in values:
    s.add(x)
  mixed(s, ops[:n], n)
with Timer("4e5 mixed ops at n = 2e5"):
  churn(s, ops)
ref = sorted(values)
for t, x in ops[:n]:
  if t == 0:
    insort(ref, x)
  elif t == 1:
    del ref[x % len(ref)]
check_eq(list(s), ref)
# The sqrt(n) bound: buckets split, so none outgrows RATIO * count.
check(max(map(len, s.a)) <= s.RATIO * len(s.a), "bucket too large")
print("OK ds2/ordered-set")
