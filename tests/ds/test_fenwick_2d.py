from testlib import Timer, check_eq, iters, load, rnd

lib = load("ds/fenwick-2d")

for _ in range(iters(1500)):
  n, m = rnd(0, 8), rnd(0, 8)
  a = [[0] * m for _ in range(n)]
  fw = lib.Fenwick2D(n, m)
  big = 10 if rnd(0, 3) else 10**18
  for _ in range(30):
    if n and m:
      x, y, v = rnd(0, n - 1), rnd(0, m - 1), rnd(-big, big)
      a[x][y] += v
      fw.add(x, y, v)
    x1 = rnd(0, n)
    x2 = rnd(x1, n)
    y1 = rnd(0, m)
    y2 = rnd(y1, m)
    want = sum(sum(row[y1:y2]) for row in a[x1:x2])
    check_eq(fw.rect(x1, y1, x2, y2), want)
    check_eq(fw.query(x2, y2), sum(sum(row[:y2]) for row in a[:x2]))

# The Tweak's rectangle add: corners at x2 = n or y2 = m are no-ops.
for _ in range(iters(800)):
  n, m = rnd(0, 8), rnd(0, 8)
  a = [[0] * m for _ in range(n)]
  fw = lib.Fenwick2D(n, m)
  for _ in range(20):
    x1, y1 = rnd(0, n), rnd(0, m)
    x2, y2, v = rnd(x1, n), rnd(y1, m), rnd(-10, 10)
    fw.add(x1, y1, v)
    fw.add(x2, y2, v)
    fw.add(x1, y2, -v)
    fw.add(x2, y1, -v)
    for i in range(x1, x2):
      for j in range(y1, y2):
        a[i][j] += v
    got = [[fw.query(i + 1, j + 1) for j in range(m)] for i in range(n)]
    check_eq(got, a)

n, m, q = 1000, 1000, 200_000
adds = [(rnd(0, n - 1), rnd(0, m - 1), rnd(-(10**9), 10**9)) for _ in range(q)]
rects = []
for _ in range(q):
  x1, y1 = rnd(0, n), rnd(0, m)
  rects.append((x1, y1, rnd(x1, n), rnd(y1, m)))
with Timer("2e5 adds + 2e5 rectangle sums on 1000 x 1000"):
  fw = lib.Fenwick2D(n, m)
  for x, y, v in adds:
    fw.add(x, y, v)
  got = [fw.rect(*r) for r in rects]
pre = [[0] * (m + 1) for _ in range(n + 1)]
for x, y, v in adds:
  pre[x + 1][y + 1] += v
for i in range(n):
  for j in range(m):
    pre[i + 1][j + 1] += pre[i][j + 1] + pre[i + 1][j] - pre[i][j]
for (x1, y1, x2, y2), g in zip(rects, got):
  want = pre[x2][y2] - pre[x1][y2] - pre[x2][y1] + pre[x1][y1]
  check_eq(g, want)
print("OK ds/fenwick-2d")
