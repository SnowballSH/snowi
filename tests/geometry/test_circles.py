import math

from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("geometry/circles")


def rnd_pt(c):
  return (rng.uniform(-c, c), rng.uniform(-c, c))


def polar(r, t):
  return (r * math.cos(t), r * math.sin(t))


def add(a, b):
  return (a[0] + b[0], a[1] + b[1])


def dist(a, b):
  return math.hypot(a[0] - b[0], a[1] - b[1])


def close(u, v, tol=1e-9):
  return abs(u - v) <= tol * max(1.0, abs(v))


def near(p, q, tol=1e-9):
  return dist(p, q) <= tol * max(1.0, math.hypot(*q))


def same_set(got, want, tol=1e-7):
  return len(got) == len(want) and all(
    any(near(g, w, tol) for g in got) for w in want
  )


def ternary(lo, hi, f):
  for _ in range(200):
    m1, m2 = lo + (hi - lo) / 3, hi - (hi - lo) / 3
    if f(m1) < f(m2):
      hi = m2
    else:
      lo = m1
  return f((lo + hi) / 2)


def by_angle(c, r, t0, x):
  """Points at angle t0 +- acos(x) on circle (c, r)."""
  if not -1 <= x <= 1:
    return []
  al = math.acos(x)
  return [add(c, polar(r, t0 + al)), add(c, polar(r, t0 - al))]


def angle(a, b):
  return math.atan2(b[1] - a[1], b[0] - a[0])


def check_line(c, r, a, b):
  v = lib.circle_line(c, r, a, b)
  ux, uy = b[0] - a[0], b[1] - a[1]
  dmin = ternary(-1e4, 1e4, lambda t: dist((a[0] + ux * t, a[1] + uy * t), c))
  if dmin > r + 1e-6:
    check_eq(v, [])
  if dmin < r - 1e-6:
    check_eq(len(v), 2)
  for q in v:
    check(close(dist(q, c), r, 1e-6))
    check(
      abs(ux * (q[1] - a[1]) - uy * (q[0] - a[0])) / math.hypot(ux, uy) < 1e-6
    )
  if len(v) == 2:
    t0 = (v[0][0] - a[0]) * ux + (v[0][1] - a[1]) * uy
    t1 = (v[1][0] - a[0]) * ux + (v[1][1] - a[1]) * uy
    check(t0 < t1)


def check_inter(a, r1, b, r2):
  v = lib.circle_inter(a, r1, b, r2)
  d = dist(a, b)
  if d > r1 + r2 + 1e-6 or d < abs(r1 - r2) - 1e-6:
    check_eq(v, [])
  if abs(r1 - r2) + 1e-6 < d < r1 + r2 - 1e-6:
    x = (r1 * r1 + d * d - r2 * r2) / (2 * r1 * d)
    check(same_set(v, by_angle(a, r1, angle(a, b), x)))
    check(
      (b[0] - a[0]) * (v[0][1] - a[1]) - (b[1] - a[1]) * (v[0][0] - a[0]) > 0
    )
  for q in v:
    check(close(dist(q, a), r1, 1e-6) and close(dist(q, b), r2, 1e-6))


def check_tangents(c1, r1, c2, r2):
  t = lib.tangents(c1, r1, c2, r2)
  d, dr = dist(c1, c2), r1 - r2
  if d < abs(dr) - 1e-6:
    check_eq(t, [])
  if d > abs(dr) + 1e-6:
    want = by_angle(c1, r1, angle(c1, c2), dr / d)
    check(same_set([p for p, _ in t], want))
  for p, q in t:
    check(close(dist(p, c1), r1, 1e-6) and close(dist(q, c2), abs(r2), 1e-6))
    L = dist(p, q)
    if L < 1e-6:
      continue
    ux, uy = (q[0] - p[0]) / L, (q[1] - p[1]) / L
    check(abs(ux * (p[0] - c1[0]) + uy * (p[1] - c1[1])) < 1e-6)
    check(abs(ux * (q[0] - c2[0]) + uy * (q[1] - c2[1])) < 1e-6)
    s1 = ux * (c1[1] - p[1]) - uy * (c1[0] - p[0])
    s2 = ux * (c2[1] - p[1]) - uy * (c2[0] - p[0])
    if r2:
      check(s1 * s2 > 0 if r2 > 0 else s1 * s2 < 0)


def brute_mec(ps):
  best = None

  def consider(o, r):
    nonlocal best
    if all(dist(p, o) <= r + 1e-9 for p in ps) and (
      best is None or r < best[1]
    ):
      best = (o, r)

  consider(ps[0], 0.0)
  n = len(ps)
  for i in range(n):
    for j in range(i + 1, n):
      a, b = ps[i], ps[j]
      o = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
      consider(o, dist(o, a))
      for k in range(j + 1, n):
        c = ps[k]
        a1, b1 = 2 * (b[0] - a[0]), 2 * (b[1] - a[1])
        e1 = b[0] ** 2 + b[1] ** 2 - a[0] ** 2 - a[1] ** 2
        a2, b2 = 2 * (c[0] - a[0]), 2 * (c[1] - a[1])
        e2 = c[0] ** 2 + c[1] ** 2 - a[0] ** 2 - a[1] ** 2
        det = a1 * b2 - a2 * b1
        if abs(det) < 1e-12:
          continue
        q = ((e1 * b2 - e2 * b1) / det, (a1 * e2 - a2 * e1) / det)
        consider(q, dist(q, a))
  return best


for _ in range(iters(2000)):
  c, a, b = rnd_pt(10), rnd_pt(10), rnd_pt(10)
  r, t = rng.uniform(0.1, 10), rng.uniform(0, 2 * math.pi)
  check_line(c, r, a, b)
  T, (dx, dy) = add(c, polar(r, t)), polar(1, t)
  s1, s2 = rng.uniform(0.5, 5), rng.uniform(0.5, 5)
  v = lib.circle_line(
    c, r, (T[0] + dy * s1, T[1] - dx * s1), (T[0] - dy * s2, T[1] + dx * s2)
  )
  check(len(v) == 1 and near(v[0], T, 1e-7))

  r1, r2 = rng.uniform(0.1, 10), rng.uniform(0.1, 10)
  check_inter(a, r1, b, r2)
  check_inter(a, r1, add(a, polar(r1 + r2 + rng.uniform(-1, 1), t)), r2)
  ext, inn = add(a, polar(r1 + r2, t)), add(a, polar(abs(r1 - r2), t))
  v = lib.circle_inter(a, r1, ext, r2)
  check(len(v) == 1 and near(v[0], add(a, polar(r1, t)), 1e-7))
  v = lib.circle_inter(a, r1, inn, r2)
  check(
    len(v) == 1 and near(v[0], add(a, polar(-r1 if r1 < r2 else r1, t)), 1e-7)
  )
  check_eq(lib.circle_inter(a, r1, a, r2), [])

  check_tangents(a, r1, b, r2)
  check_tangents(a, r1, b, -r2)
  check_tangents(a, r1, b, 0.0)
  check_tangents(a, r1, add(a, polar(abs(r1 - r2) + rng.uniform(-1, 1), t)), r2)
  check_eq(len(lib.tangents(a, r1, ext, r2)), 2)
  check_eq(len(lib.tangents(a, r1, ext, -r2)), 1)
  check_eq(len(lib.tangents(a, r1, inn, r2)), 1)
  check_eq(len(lib.tangents(a, r1, add(a, polar(r1, t)), 0.0)), 1)
  check_eq(lib.tangents(a, r1, a, r2), [])

  p, q, s = rnd_pt(10), rnd_pt(10), rnd_pt(10)
  if abs((q[0] - p[0]) * (s[1] - p[1]) - (q[1] - p[1]) * (s[0] - p[0])) > 1e-2:
    o, rad = lib.circumcircle(p, q, s)
    check(all(close(dist(x, o), rad, 1e-7) for x in (p, q, s)))

check_eq(lib.circle_line((0, 0), 5, (5, -1), (5, 1)), [(5.0, 0.0)])
check_eq(lib.circle_line((0, 0), 5, (6, -1), (6, 1)), [])
check(same_set(lib.circle_inter((0, 0), 5, (6, 0), 5), [(3, -4), (3, 4)]))

for it in range(iters(800)):
  n = rnd(1, 8)
  if it % 2:
    ps = [rnd_pt(10) for _ in range(n)]
  else:
    ps = [(rnd(-3, 3), rnd(-3, 3)) for _ in range(n)]
  if it % 5 == 0:
    ps.append(ps[0])
  ctr, rad = lib.mec(ps)
  bc, br = brute_mec(ps)
  check(close(rad, br, 1e-7), f"{ps}: {rad} vs {br}")
  check(near(ctr, bc, 1e-6))
  check(all(dist(p, ctr) <= rad * (1 + 1e-7) + 1e-12 for p in ps))

# Far from the origin the precision must follow the spread of the
# points, not their position. Reproducer: ZeroDivisionError at 1e9
# (a circle through a collinear triple after a rounded test).
repro = [(0, 2), (1, 0), (2, 1), (0, 2), (0, 2)]
for it in range(iters(1500)):
  fixed = it < iters(300)
  B = (10**9, 10**9) if fixed else (rnd(-(10**9), 10**9), rnd(-(10**9), 10**9))
  spread = (1, 3, 10, 100)[it % 4]
  rel = repro
  if not fixed:
    m = rnd(1, 8)
    if it % 3:
      rel = [(rnd(0, spread), rnd(0, spread)) for _ in range(m)]
    else:
      rel = [(rng.uniform(0, spread), rng.uniform(0, spread)) for _ in range(m)]
    if it % 2:
      rel.append(rel[rnd(0, m - 1)])
  ps = [(B[0] + x, B[1] + y) for x, y in rel]
  rel = [(x - B[0], y - B[1]) for x, y in ps]  # exact (Sterbenz)
  ctr, rad = lib.mec(ps)
  bc, br = brute_mec(rel)
  check(close(rad, br, 1e-7), f"{ps}: {rad} vs {br}")
  check(near((ctr[0] - B[0], ctr[1] - B[1]), bc, 1e-6))
  # The returned center is rounded to about 1e-7 at 1e9.
  check(all(dist(p, ctr) <= rad * (1 + 1e-7) + 1e-6 for p in ps))

n = 200_000
ps = [
  polar(rng.uniform(0, 1e3) if i % 10 else 1e3, rng.uniform(0, 2 * math.pi))
  for i in range(n)
]
with Timer("mec of 2e5 points"):
  ctr, rad = lib.mec(ps)
check(999 < rad <= 1e3 * (1 + 1e-9))
check(all(dist(p, ctr) <= rad * (1 + 1e-7) for p in ps))
with Timer("2e5 each of circle_line, circle_inter, tangents"):
  total = 0
  for i in range(n):
    p, q = ps[i - 1], ps[i]
    total += len(lib.circle_line(p, 1.0, q, (0, 0)))
    total += len(lib.circle_inter(p, 3.0, q, 4.0))
    total += len(lib.tangents(p, 3.0, q, -4.0))
print("OK geometry/circles")
