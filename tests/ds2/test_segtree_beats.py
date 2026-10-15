from testlib import Timer, check_eq, iters, load, rnd, rnd_list

lib = load("ds2/segtree-beats")
INF = lib.Beats.INF


def brute_query(b):
  return (sum(b), max(b), min(b)) if b else (0, -INF, INF)


def check_root(st, a):  # every field of the root node
  if not a:
    return
  b = sorted(a)
  mx, mn = b[-1], b[0]
  below = [v for v in b if v < mx]
  above = [v for v in b if v > mn]
  got = (st.MX[1], st.CMX[1], st.MX2[1], st.MN[1], st.CMN[1], st.MN2[1])
  want = (mx, b.count(mx), below[-1] if below else -INF)
  want += (mn, b.count(mn), above[0] if above else INF)
  check_eq(got, want)


def stress(rounds, maxn, ops, c, d):
  for _ in range(iters(rounds)):
    n = rnd(0, maxn)
    a = rnd_list(n, -c, c)
    st = lib.Beats(a)
    for _ in range(ops):
      l, r = sorted((rnd(0, n), rnd(0, n)))
      kind, x = rnd(0, 4), rnd(-c, c)
      if kind == 0:
        st.chmin(l, r, x)
        a[l:r] = [min(v, x) for v in a[l:r]]
      elif kind == 1:
        st.chmax(l, r, x)
        a[l:r] = [max(v, x) for v in a[l:r]]
      elif kind == 2:
        x = rnd(-d, d)
        st.add(l, r, x)
        a[l:r] = [v + x for v in a[l:r]]
      elif kind == 3:  # Tweak: assign
        st.chmin(l, r, x)
        st.chmax(l, r, x)
        a[l:r] = [x] * (r - l)
      else:
        check_eq(st.query(l, r), brute_query(a[l:r]))
      check_root(st, a)
    check_eq([st.query(i, i + 1)[0] for i in range(n)], a)


stress(1200, 12, 40, 3, 2)  # few distinct values
stress(600, 16, 40, 10, 10)
stress(100, 64, 150, 1000, 300)
stress(200, 8, 40, 2 * 10**17, 10**15)  # sums up to 1.6e18
stress(1, 500, 2000, 10**9, 10**6)  # medium

st = lib.Beats([])  # n = 0
check_eq(st.query(0, 0), (0, -INF, INF))
st.chmin(0, 0, 1)
st.add(0, 0, 1)
st = lib.Beats([7])  # n = 1
st.chmax(0, 1, 9)
st.add(0, 1, -2)
st.chmin(0, 1, 3)
check_eq(st.query(0, 1), (3, 3, 3))


def run(st, ops):  # timed parts run in functions, as in main()
  total = 0
  for kind, l, r, x in ops:
    if kind == 0:
      st.chmin(l, r, x)
    elif kind == 1:
      st.chmax(l, r, x)
    elif kind == 2:
      st.add(l, r, x - 500_000_000)
    else:
      total += st.query(l, r)[0]
  return total


def random_ops(n, q, kinds):
  ops = []
  for _ in range(q):
    l = rnd(0, n - 1)
    ops.append((kinds[rnd(0, 3)], l, rnd(l + 1, n), rnd(0, 10**9)))
  return ops


n = 100_000
a = rnd_list(n, 0, 10**9)
ops = random_ops(n, n, [0, 3, 0, 3])
with Timer("1e5 chmin + sum at n = 1e5"):
  st = lib.Beats(a)
  run(st, ops)
for l, r in [(0, n), (5, 77_777), (n - 1, n)]:
  b = []
  for i in range(l, r):  # every element, read one at a time
    b.append(st.query(i, i + 1)[0])
  check_eq(st.query(l, r), brute_query(b))
m = 50_000  # the documented PyPy limit; warm here, about 2.8 s cold
ops = random_ops(m, m, [0, 1, 2, 3])
with Timer("5e4 mixed chmin, chmax, add, sum at n = 5e4"):
  st = lib.Beats(a[:m])
  run(st, ops)
print("OK ds2/segtree-beats")
