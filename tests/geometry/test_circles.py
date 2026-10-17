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


def in_poly(v, p):
  inside = False
  for i in range(len(v)):
    a, b = v[i - 1], v[i]
    if (a[1] > p[1]) != (b[1] > p[1]):
      x = a[0] + (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1])
      inside ^= p[0] < x
  return inside


def signed_area(v):
  return sum(cross(v[i - 1], v[i]) for i in range(len(v))) / 2


def cross(a, b):
  return a[0] * b[1] - a[1] * b[0]


def green_area(c, r, v):
  """Disc (c, r) inside a simple polygon by Green's theorem, another
  method than circle_poly: polygon edges cut at the circle (pieces in
  the disc), then circle arcs between crossings (arcs in the
  polygon). Generic input only."""
  sign = 1
  if signed_area(v) < 0:
    v, sign = v[::-1], -1
  v = [(x - c[0], y - c[1]) for x, y in v]
  total, angles = 0.0, []
  for i in range(len(v)):
    a, b = v[i - 1], v[i]
    d = (b[0] - a[0], b[1] - a[1])
    A = d[0] ** 2 + d[1] ** 2
    B = 2 * (a[0] * d[0] + a[1] * d[1])
    C = a[0] ** 2 + a[1] ** 2 - r * r
    ts, disc = [0.0, 1.0], B * B - 4 * A * C
    if disc > 0:
      for sg in (-1, 1):
        t = (-B + sg * math.sqrt(disc)) / (2 * A)
        if 0 < t < 1:
          ts.append(t)
          angles.append(math.atan2(a[1] + d[1] * t, a[0] + d[0] * t))
    ts.sort()
    for t0, t1 in zip(ts, ts[1:]):
      p = (a[0] + d[0] * t0, a[1] + d[1] * t0)
      q = (a[0] + d[0] * t1, a[1] + d[1] * t1)
      if math.hypot((p[0] + q[0]) / 2, (p[1] + q[1]) / 2) < r:
        total += cross(p, q) / 2
  if not angles:
    return sign * (total + (math.pi * r * r if in_poly(v, (r, 0)) else 0))
  angles.sort()
  for i, lo in enumerate(angles):
    hi = angles[i + 1] if i + 1 < len(angles) else angles[0] + 2 * math.pi
    if in_poly(v, polar(r, (lo + hi) / 2)):
      total += r * r * (hi - lo) / 2
  return sign * total


def clip(v, w):
  """v clipped to the left of every edge of the convex polygon w."""
  for k in range(len(w)):
    s, e = w[k - 1], w[k]
    ex, ey, out = e[0] - s[0], e[1] - s[1], []
    for i in range(len(v)):
      prev, cur = v[i - 1], v[i]
      a = ex * (cur[1] - s[1]) - ey * (cur[0] - s[0])
      b = ex * (prev[1] - s[1]) - ey * (prev[0] - s[0])
      if (a < 0) != (b < 0):
        f = a / (a - b)
        out.append(
          (cur[0] + (prev[0] - cur[0]) * f, cur[1] + (prev[1] - cur[1]) * f)
        )
      if a >= 0:
        out.append(cur)
    v = out
  return v


def ngon_bounds(c, r, v, m=128):
  """Areas inside regular m-gons inscribed in and around the circle."""
  ins = [add(c, polar(r, 2 * math.pi * k / m)) for k in range(m)]
  R = r / math.cos(math.pi / m)
  out = [add(c, polar(R, 2 * math.pi * k / m)) for k in range(m)]
  x, y = signed_area(clip(v, ins)), signed_area(clip(v, out))
  return min(x, y), max(x, y)


def rnd_star(o, n, rmax):
  """Star-shaped (so simple) counter-clockwise polygon around o."""
  while True:
    t = sorted(rng.uniform(0, 2 * math.pi) for _ in range(n))
    gaps = [t[0] + 2 * math.pi - t[-1]] + [y - x for x, y in zip(t, t[1:])]
    if max(gaps) < 3:
      return [add(o, polar(rng.uniform(0.05, rmax), x)) for x in t]


def simple(v):
  n = len(v)
  if abs(signed_area(v)) < 1e-9 or any(v[i - 1] == v[i] for i in range(n)):
    return False
  for i in range(n):
    for j in range(i + 2, n):
      if i == 0 and j == n - 1:
        continue
      a, b, p, q = v[i], v[(i + 1) % n], v[j], v[(j + 1) % n]
      d1 = (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])
      d2 = (b[0] - a[0]) * (q[1] - a[1]) - (b[1] - a[1]) * (q[0] - a[0])
      d3 = (q[0] - p[0]) * (a[1] - p[1]) - (q[1] - p[1]) * (a[0] - p[0])
      d4 = (q[0] - p[0]) * (b[1] - p[1]) - (q[1] - p[1]) * (b[0] - p[0])
      if d1 * d2 <= 0 and d3 * d4 <= 0:
        return False
  return True


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

# circle_poly: generic polygons against Green's theorem.
for it in range(iters(2000)):
  c, v = rnd_pt(5), rnd_star(rnd_pt(5), rnd(3, 9), 6)
  if it % 2:
    v.reverse()
  r = (0.3, 1, 3, 8, 20)[it % 5] * rng.uniform(0.5, 1)
  got, want = lib.circle_poly(c, r, v), green_area(c, r, v)
  check(close(got, want, 1e-7), f"{c} {r} {v}: {got} vs {want}")
  check(close(lib.circle_poly(c, r, v[::-1]), -want, 1e-7))

# Degenerate integer cases (vertices on the circle, tangent edges, the
# center on a vertex or an edge): inside regular-polygon bounds and
# close to the area for a slightly larger radius.
degenerate = 0
for it in range(iters(400)):
  if it % 4 == 0:
    v = [(0, 0), (6, 0), (6, 6), (0, 6)]
  else:
    v = [(rnd(-5, 5), rnd(-5, 5)) for _ in range(rnd(3, 6))]
  if it % 8 == 4:
    v.reverse()
  if not simple(v):
    continue
  mid = ((v[0][0] + v[1][0]) / 2, (v[0][1] + v[1][1]) / 2)
  c = (rnd(-3, 3), rnd(-3, 3)) if it % 2 else (v[0], mid, (3, 3))[it // 2 % 3]
  r = rnd(1, 5)
  lo, hi = ngon_bounds(c, r, v)
  got = lib.circle_poly(c, r, v)
  check(lo - 1e-9 <= got <= hi + 1e-9, f"{c} {r} {v}: {got}")
  check(close(got, green_area(c, r + 1e-8, v), 1e-6))
  degenerate += 1
check(degenerate > iters(400) // 4)

sq = [(0, 0), (2, 0), (2, 2), (0, 2)]
check(close(lib.circle_poly((0, 0), 1, sq), math.pi / 4, 1e-12))
check(close(lib.circle_poly((1, 1), 1, sq), math.pi, 1e-12))
check(close(lib.circle_poly((1, 1), 5, sq), 4, 1e-12))
check(close(lib.circle_poly((9, 9), 1, sq), 0, 1e-12))
check(close(lib.circle_poly((1, 0), 1, sq), math.pi / 2, 1e-12))

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
with Timer("circle_poly of 2e5 vertices"):
  lib.circle_poly((1, 2), 500.0, ps)
print("OK geometry/circles")
