from testlib import Timer, check_eq, iters, load, rnd, rng

lib = load("ds/li-chao")
INF = lib.LiChao.INF


def brute(segs, x):
  y = INF
  for k, m, l, r in segs:
    if l <= x < r:
      y = min(y, k * x + m)
  return y


def rnd_xs(n, c, distinct):
  xs = sorted(rnd(-c, c) for _ in range(n))
  return sorted(set(xs)) if distinct else xs


def stress(rounds, cx, ck, cm):
  for _ in range(iters(rounds)):
    xs = rnd_xs(rnd(1, 14), cx, rnd(0, 3) > 0)
    lc = lib.LiChao(xs)
    segs = []
    for _ in range(40):
      k, m = rnd(-ck, ck), rnd(-cm, cm)
      if rnd(0, 1):
        segs.append((k, m, -(10**30), 10**30))
        lc.add_line(k, m)
      else:
        l, r = rnd(-cx - 2, cx + 2), rnd(-cx - 2, cx + 2)
        if l > r and rnd(0, 3):  # else empty
          l, r = r, l
        segs.append((k, m, l, r))
        lc.add_segment(k, m, l, r)
      x = rng.choice(xs)
      check_eq(lc.query(x), brute(segs, x))
    for x in xs:
      check_eq(lc.query(x), brute(segs, x))


stress(2000, 8, 3, 10)
stress(600, 30, 10, 100)
stress(600, 10**9, 10**9, 10**18)

lc = lib.LiChao([7])  # one point, no line yet
check_eq(lc.query(7), INF)
lc.add_segment(1, 1, 8, 9)
check_eq(lc.query(7), INF)
lc.add_segment(2, -3, 7, 8)
check_eq(lc.query(7), 11)

# Maximum by negation (the Tweak).
lc = lib.LiChao(list(range(-5, 6)))
for k, m in [(1, 0), (-1, 0), (0, 2)]:
  lc.add_line(-k, -m)
want = [max(x, -x, 2) for x in range(-5, 6)]
check_eq([-lc.query(x) for x in range(-5, 6)], want)

n = 200_000
with Timer("2e5 lines + 2e5 queries, n = 2e5"):
  # Tangents -2 t x + t^2 of -x^2 in random order give -x^2.
  ts = list(range(n))
  rng.shuffle(ts)
  lc = lib.LiChao(list(range(n)))
  for t in ts:
    lc.add_line(-2 * t, t * t)
  for x in range(n):
    if lc.query(x) != -x * x:
      check_eq(lc.query(x), -x * x)
segs = []
for _ in range(100_000):
  l = rnd(0, n - 1)
  r = rnd(l + 1, n)
  segs.append((rnd(-(10**6), 10**6), rnd(-(10**12), 10**12), l, r))
with Timer("1e5 segments + 1e5 queries, n = 2e5"):
  lc = lib.LiChao(list(range(n)))
  for k, m, l, r in segs:
    lc.add_segment(k, m, l, r)
  got = [lc.query(x) for x in range(0, n, 2)]
for x in range(0, n, 998):
  check_eq(got[x // 2], brute(segs, x))
print("OK ds/li-chao")
