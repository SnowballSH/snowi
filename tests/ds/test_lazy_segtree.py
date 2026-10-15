from functools import reduce
from operator import add

from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("ds/lazy-segtree")
MOD = 1_000_003
INF = 10**18


# Range affine on (sum, len) nodes: tags compose non-commutatively.
def sum_op(x, y):
  return (x[0] + y[0]) % MOD, x[1] + y[1]


def affine_map(f, s):
  return (f[0] * s[0] + f[1] * s[1]) % MOD, s[1]


def affine_comp(f, g):  # g, then f
  return f[0] * g[0] % MOD, (f[0] * g[1] + f[1]) % MOD


# Runs of equal values (first, last, runs) under range add: a
# non-commutative op. The identity has runs = 0.
def run_op(x, y):
  if not x[2]:
    return y
  if not y[2]:
    return x
  return x[0], y[1], x[2] + y[2] - (x[1] == y[0])


def run_map(f, s):
  return (s[0] + f, s[1] + f, s[2]) if s[2] else s


# Range assign (None = no tag), range max.
def assign_comp(f, g):
  return g if f is None else f


def assign_map(f, s):
  return s if f is None else f


def brute_max_right(a, l, op, e, pred):
  acc, r = e, l
  while r < len(a) and pred(op(acc, a[r])):
    acc = op(acc, a[r])
    r += 1
  return r


def brute_min_left(a, r, op, e, pred):
  acc, l = e, r
  while l > 0 and pred(op(a[l - 1], acc)):
    acc = op(a[l - 1], acc)
    l -= 1
  return l


def stress(monoid, leaf, tag, preds, rounds):
  op, e, mapping = monoid[:3]
  for _ in range(iters(rounds)):
    n = rnd(0, 18)
    a = [leaf() for _ in range(n)]
    st = lib.LazySegTree(a, *monoid)
    for _ in range(25):
      if n and rnd(0, 2) == 0:
        i = rnd(0, n - 1)
        a[i] = leaf()
        st.set(i, a[i])
      else:
        l = rnd(0, n)
        r = rnd(l, n)
        f = tag()
        a[l:r] = [mapping(f, s) for s in a[l:r]]
        st.apply(l, r, f)
      l = rnd(0, n)
      r = rnd(l, n)
      check_eq(st.query(l, r), reduce(op, a[l:r], e))
      check_eq(st.query_all(), reduce(op, a, e))
      if n:
        i = rnd(0, n - 1)
        check_eq(st.get(i), a[i])
      i = rnd(0, n)
      for pred in preds():
        want = brute_max_right(a, i, op, e, pred)
        check_eq(st.max_right(i, pred), want)
        want = brute_min_left(a, i, op, e, pred)
        check_eq(st.min_left(i, pred), want)


def min_preds():
  x = rnd(-8, 8)
  return [lambda v: v > x, lambda v: v >= x]


def max_preds():
  x = rnd(-8, 8)
  return [lambda v: v < x]


def run_preds():
  k = rnd(1, 4)
  return [lambda s: s[2] <= k]


def no_preds():
  return []


# (sum, len) packed into one int, as in the snippet's Tweak.
M = (1 << 20) - 1


def packed_map(f, s):
  return s + (f * (s & M) << 20)


def packed_preds():  # values >= 0, so sums grow with the range
  x = rnd(0, 40)
  return [lambda s: s >> 20 <= x]


def rnd_assign():
  return rnd(-5, 5) if rnd(0, 3) else None


stress(
  (min, INF, add, add, 0),
  lambda: rnd(-5, 5),
  lambda: rnd(-3, 3),
  min_preds,
  800,
)
stress(
  (sum_op, (0, 0), affine_map, affine_comp, (1, 0)),
  lambda: (rnd(0, 6), 1),
  lambda: (rnd(0, 3), rnd(0, 3)),
  no_preds,
  500,
)
stress(
  (run_op, (0, 0, 0), run_map, add, 0),
  lambda: (lambda v: (v, v, 1))(rnd(0, 2)),
  lambda: rnd(-1, 1),
  run_preds,
  500,
)
stress(
  (max, -INF, assign_map, assign_comp, None),
  lambda: rnd(-5, 5),
  rnd_assign,
  max_preds,
  500,
)
stress(
  (add, 0, packed_map, add, 0),
  lambda: (rnd(0, 6) << 20) + 1,
  lambda: rnd(0, 3),
  packed_preds,
  300,
)
stress(
  (add, 0, packed_map, add, 0),
  lambda: (rnd(-6, 6) << 20) + 1,
  lambda: rnd(-3, 3),
  no_preds,
  300,
)

# The Tweak's bound: with |sum| < 2^43 every packed value is a machine
# int (PyPy makes larger ones slow big ints) and stays exact.
for _ in range(iters(200)):
  n = rnd(1, 18)
  a = [(rnd(-(2**37), 2**37) << 20) + 1 for _ in range(n)]
  st = lib.LazySegTree(a, add, 0, packed_map, add, 0)
  for _ in range(25):
    l = rnd(0, n)
    r = rnd(l, n)
    f = rnd(-(2**30), 2**30)
    a[l:r] = [packed_map(f, s) for s in a[l:r]]
    st.apply(l, r, f)
    l = rnd(0, n)
    r = rnd(l, n)
    check_eq(st.query(l, r), sum(a[l:r]))
  check_eq([st.get(i) for i in range(n)], a)
  check(all(-(2**63) <= x < 2**63 for x in st.d))
  check(abs(st.query_all() >> 20) < 2**43)

with Timer("2e5 range adds + 2e5 range mins"):
  n = 200_000
  st = lib.LazySegTree([0] * n, min, INF, add, add, 0)
  ops = [(rnd(0, n - 1), rnd(0, n - 1), rnd(-5, 5)) for _ in range(n)]
  for l, r, f in ops:
    if l > r:
      l, r = r, l
    st.apply(l, r + 1, f)
    st.query(l, r + 1)
  diff = [0] * (n + 1)
  for l, r, f in ops:
    if l > r:
      l, r = r, l
    diff[l] += f
    diff[r + 1] -= f
  run = 0
  for i in range(n):
    run += diff[i]
    diff[i] = run
  check_eq([st.get(i) for i in range(n)], diff[:n])
  check_eq(st.query_all(), min(diff[:n]))
print("OK ds/lazy-segtree")
