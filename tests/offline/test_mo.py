from testlib import Timer, check, check_eq, iters, load, rnd, rnd_list

lib = load("offline/mo")


def run(a, qs, V, hil):
  """Runs mo with callbacks that check the window: every add or
  remove happens at an end, and answer(q) sees exactly qs[q]. An
  empty window has no position: removing its last element from
  either end looks the same."""
  n = len(a)
  win = [0, 0]  # lo, hi
  cnt, state = [0] * V, [0, 0, 0]  # distinct, sum of cnt^2, moves
  res = [None] * len(qs)

  def add(i):
    lo, hi = win
    check(0 <= i < n)
    if lo == hi:
      win[:] = [i, i + 1]
    else:
      check(i in (lo - 1, hi))
      win[:] = [lo, hi + 1] if i == hi else [lo - 1, hi]
    state[0] += cnt[a[i]] == 0
    state[1] += 2 * cnt[a[i]] + 1
    state[2] += 1
    cnt[a[i]] += 1

  def remove(i):
    lo, hi = win
    check(lo < hi and i in (lo, hi - 1))
    win[:] = [lo + 1, hi] if i == lo else [lo, hi - 1]
    cnt[a[i]] -= 1
    state[0] -= cnt[a[i]] == 0
    state[1] -= 2 * cnt[a[i]] + 1
    state[2] += 1

  def answer(q):
    l, r = qs[q]
    lo, hi = win
    check(l == r if lo == hi else (lo, hi) == (l, r))
    check(res[q] is None)  # answered once
    res[q] = (state[0], state[1])

  lib.mo(n, qs, add, remove, answer, hil)
  return res, state[2]


def brute(a, l, r):
  c = {}
  for v in a[l:r]:
    c[v] = c.get(v, 0) + 1
  return len(c), sum(k * k for k in c.values())


def stress(rounds, maxn, maxv):
  for _ in range(iters(rounds)):
    n, q = rnd(0, maxn), rnd(0, 30)
    a = rnd_list(n, 0, maxv - 1)
    qs = [tuple(sorted((rnd(0, n), rnd(0, n)))) for _ in range(q)]
    for hil in (False, True):
      res, _ = run(a, qs, maxv, hil)
      check_eq(res, [brute(a, l, r) for l, r in qs])


# The first 4^k cells of the curve fill the 2^k square at the
# origin, and consecutive cells are neighbours.
side = 1 << 5
at = [None] * (side * side)
for x in range(side):
  for y in range(side):
    d = lib.hilbert(x, y)
    check(d < side * side and at[d] is None)
    at[d] = (x, y)
for d in range(1, side * side):
  (x0, y0), (x1, y1) = at[d - 1], at[d]
  check_eq(abs(x0 - x1) + abs(y0 - y1), 1)
check_eq(lib.hilbert(2**21 - 1, 0), 2**42 - 1)  # last cell

stress(1500, 12, 3)
stress(600, 40, 10)
stress(100, 300, 1000)
check_eq(run([], [(0, 0), (0, 0)], 1, False), ([(0, 0), (0, 0)], 0))
check_eq(run([5], [], 6, True), ([], 0))


def distinct(a, qs, V, hil):  # Usage, timed
  cnt, res, cur = [0] * V, [0] * len(qs), [0]

  def add(i):
    cur[0] += cnt[a[i]] == 0
    cnt[a[i]] += 1

  def remove(i):
    cnt[a[i]] -= 1
    cur[0] -= cnt[a[i]] == 0

  def answer(q):
    res[q] = cur[0]

  lib.mo(len(a), qs, add, remove, answer, hil)
  return res


n, V = 100_000, 1000
a = rnd_list(n, 0, V - 1)
qs = []
for _ in range(n):
  l = rnd(0, n - 1)
  qs.append((l, rnd(l + 1, n)))
with Timer("distinct values, n = q = 1e5, blocks"):
  got = distinct(a, qs, V, False)
with Timer("distinct values, n = q = 1e5, Hilbert"):
  check_eq(distinct(a, qs, V, True), got)
for j in range(0, n, 2000):
  l, r = qs[j]
  check_eq(got[j], len(set(a[l:r])))
# moves stay near n sqrt(q) for both orders
m = 20_000
for hil in (False, True):
  _, moves = run(a[:m], [(l // 5, r // 5) for l, r in qs[:m]], V, hil)
  check(moves < 1.2 * m * m**0.5, f"{moves} moves")
print("OK offline/mo")
