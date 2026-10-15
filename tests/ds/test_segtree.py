from functools import reduce
from operator import add

from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("ds/segtree")
MOD = 1_000_003


def compose(f, g):  # x -> b x + c: f first, then g
  return g[0] * f[0] % MOD, (g[0] * f[1] + g[1]) % MOD


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


def stress(op, e, gen, preds, rounds):
  for _ in range(iters(rounds)):
    n = rnd(0, 18)
    a = [gen() for _ in range(n)]
    st = lib.SegTree(a, op, e)
    for _ in range(25):
      if n and rnd(0, 2) == 0:
        i = rnd(0, n - 1)
        a[i] = gen()
        st.set(i, a[i])
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


def max_preds():
  x = rnd(-6, 6)
  return [lambda v: v < x, lambda v: v <= x]


def sum_preds():  # values >= 0, so prefix sums are monotone
  x = rnd(-1, 30)
  return [lambda v: v <= x]


def str_preds():
  k = rnd(0, 8)
  return [lambda s: "ab" not in s, lambda s: s.count("b") <= k]


def no_preds():
  return []


stress(max, -(10**18), lambda: rnd(-5, 5), max_preds, 1200)
stress(add, 0, lambda: rnd(0, 5), sum_preds, 800)
stress(add, "", lambda: rng.choice("ab"), str_preds, 800)
stress(compose, (1, 0), lambda: (rnd(0, 5), rnd(0, 5)), no_preds, 400)

# The Tweak's packing, (sum << 20) + count under add: with sums below
# 2^43 every node stays a machine int and unpacks exactly.
for _ in range(iters(200)):
  n = rnd(1, 18)
  a = [(rnd(-(2**38), 2**38) << 20) + 1 for _ in range(n)]
  st = lib.SegTree(a, add, 0)
  for _ in range(20):
    i = rnd(0, n - 1)
    a[i] = (rnd(-(2**38), 2**38) << 20) + 1
    st.set(i, a[i])
    l = rnd(0, n)
    r = rnd(l, n)
    got = st.query(l, r)
    check_eq(
      (got >> 20, got & (2**20 - 1)), (sum(x >> 20 for x in a[l:r]), r - l)
    )
  check(all(-(2**63) <= x < 2**63 for x in st.t))

with Timer("2e5 sets + 2e5 queries + 2e5 max_right, max"):
  n = 200_000
  a = [rnd(0, 10**9) for _ in range(n)]
  st = lib.SegTree(a, max, -1)
  for i in range(n):
    st.set(i, i)
  total = 0
  for i in range(n):
    total += st.query(i // 2, i + 1)
  check_eq(total, n * (n - 1) // 2)
  for i in range(n):
    total += st.max_right(i, lambda v: v < n - 1)
  check_eq(total, n * (n - 1) // 2 + (n - 1) * n)
print("OK ds/segtree")
