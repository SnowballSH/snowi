import math
from fractions import Fraction

from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("geometry/lines")


def rnd_pt(c):
  return (rnd(-c, c), rnd(-c, c))


def close(u, v, tol=1e-9):
  return abs(u - v) <= tol * max(1.0, abs(v))


def brute_side(s, e, p):
  """Orientation from angles, without cross products."""
  if p == s:
    return 0
  t = math.sin(
    math.atan2(p[1] - s[1], p[0] - s[0]) - math.atan2(e[1] - s[1], e[0] - s[0])
  )
  return 0 if abs(t) < 1e-9 else (1 if t > 0 else -1)


def eq(a, b):
  """Line through a and b as A x + B y = C."""
  A, B = b[1] - a[1], a[0] - b[0]
  return A, B, A * a[0] + B * a[1]


def brute_line(a, b, c, d):
  """Cramer's rule on the line equations, exact with Fractions."""
  A1, B1, C1 = eq(a, b)
  A2, B2, C2 = eq(c, d)
  det = A1 * B2 - A2 * B1
  if det == 0:
    same = A1 * C2 == A2 * C1 and B1 * C2 == B2 * C1
    return (-1 if same else 0), None
  return 1, (Fraction(C1 * B2 - C2 * B1, det), Fraction(A1 * C2 - A2 * C1, det))


def in_box(q, a, b):
  return min(a[0], b[0]) <= q[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= q[
    1
  ] <= max(a[1], b[1])


def exact_on(s, e, p):
  if s == e:
    return p == s
  A, B, C = eq(s, e)
  return A * p[0] + B * p[1] == C and in_box(p, s, e)


def brute_seg(a, b, c, d):
  if a != b and c != d:
    kind, q = brute_line(a, b, c, d)
    if kind == 1:
      return [q] if in_box(q, a, b) and in_box(q, c, d) else []
  return sorted(
    {p for p in (a, b, c, d) if exact_on(a, b, p) and exact_on(c, d, p)}
  )


def same_pt(p, q):
  return close(p[0], float(q[0]), 1e-12) and close(p[1], float(q[1]), 1e-12)


def ternary(lo, hi, f):
  for _ in range(200):
    m1, m2 = lo + (hi - lo) / 3, hi - (hi - lo) / 3
    if f(m1) < f(m2):
      hi = m2
    else:
      lo = m1
  return f((lo + hi) / 2)


def check_seg(a, b, c, d):
  got, want = lib.seg_inter(a, b, c, d), brute_seg(a, b, c, d)
  check_eq(len(got), len(want), f"seg_inter{(a, b, c, d)}")
  for p, q in zip(got, want):
    check(same_pt(p, q), f"seg_inter{(a, b, c, d)}: {got} vs {want}")


def check_all(a, b, c, d, p):
  if a != b:
    check_eq(lib.side_of(a, b, p), brute_side(a, b, p))
  check_eq(lib.on_segment(a, b, p), exact_on(a, b, p))
  check_seg(a, b, c, d)

  def dist_at(t):
    return math.hypot(
      a[0] + (b[0] - a[0]) * t - p[0], a[1] + (b[1] - a[1]) * t - p[1]
    )

  check(close(lib.seg_dist(a, b, p), ternary(0, 1, dist_at), 1e-7))
  if a == b:
    return
  ld = lib.line_dist(a, b, p)
  check(close(abs(ld), ternary(-200, 200, dist_at), 1e-7))
  check_eq(1 if ld > 1e-9 else (-1 if ld < -1e-9 else 0), brute_side(a, b, p))
  q, r = lib.line_proj(a, b, p), lib.line_proj(a, b, p, True)
  vx, vy = b[0] - a[0], b[1] - a[1]
  check(abs(vx * (q[1] - a[1]) - vy * (q[0] - a[0])) < 1e-9)
  check(abs(vx * (p[0] - q[0]) + vy * (p[1] - q[1])) < 1e-9)
  check(close(math.dist(p, q), abs(ld)))
  check(close(r[0], 2 * q[0] - p[0]) and close(r[1], 2 * q[1] - p[1]))
  check(close(math.dist(r, a), math.dist(p, a)))
  if c == d:
    return
  kind, x = lib.line_inter(a, b, c, d)
  bk, bq = brute_line(a, b, c, d)
  check_eq(kind, bk)
  if kind == 1:
    check(same_pt(x, bq))


for it in range(iters(8000)):
  C = 4 if it % 3 else 30
  a, b, c, d, p = (rnd_pt(C) for _ in range(5))
  if it % 5 == 0:
    k1, k2 = rnd(-2, 2), rnd(-2, 2)
    c = (a[0] + (b[0] - a[0]) * k1, a[1] + (b[1] - a[1]) * k1)
    d = (a[0] + (b[0] - a[0]) * k2, a[1] + (b[1] - a[1]) * k2)
  if it % 7 == 0:
    k = rnd(-1, 2)
    p = (a[0] + (b[0] - a[0]) * k, a[1] + (b[1] - a[1]) * k)
  check_all(a, b, c, d, p)

# Hand-made: T-junction, touching ends, overlap, disjoint collinear,
# identical, parallel, point segments.
check_eq(lib.seg_inter((0, 0), (4, 0), (2, 0), (2, 2)), [(2, 0)])
check_eq(lib.seg_inter((0, 0), (2, 0), (2, 0), (4, 0)), [(2, 0)])
check_eq(lib.seg_inter((0, 0), (1, 0), (2, 0), (4, 0)), [])
check_eq(lib.seg_inter((0, 0), (4, 0), (4, 0), (0, 0)), [(0, 0), (4, 0)])
check_eq(lib.seg_inter((0, 0), (2, 0), (0, 2), (2, 2)), [])
check_eq(lib.seg_inter((1, 1), (1, 1), (0, 0), (2, 2)), [(1, 1)])
check_eq(lib.seg_inter((1, 1), (1, 1), (0, 0), (2, 0)), [])
check_eq(lib.line_inter((0, 0), (2, 0), (4, 0), (0, 0)), (-1, None))
check(lib.on_segment((0, 0), (0, 0), (0, 0)))
check(not lib.on_segment((0, 0), (0, 0), (1, 1)))
check_eq(lib.seg_dist((0, 0), (0, 0), (3, 4)), 5.0)
check_eq(lib.line_dist((0, 0), (2, 0), (7, -3)), -3.0)

# Extremes: Python ints are exact at any size.
C = 10**9
for it in range(iters(3000)):
  a, b, c, d = (rnd_pt(C) for _ in range(4))
  if it % 2:
    v, a = rnd_pt(1000), rnd_pt(C - 2 * 10**8)
    ks = [rnd(-(10**5), 10**5) for _ in range(3)]
    b, c, d = [(a[0] + v[0] * k, a[1] + v[1] * k) for k in ks]
  if it % 3 == 0:
    c, d = (a[0], b[1]), (b[0], a[1])
  check_seg(a, b, c, d)
  check_eq(lib.on_segment(a, b, c), exact_on(a, b, c))


def far_error(s1, e1, s2, e2):
  """Max coordinate error of line_inter against exact Fractions."""
  kind, x = lib.line_inter(s1, e1, s2, e2)
  bk, bq = brute_line(
    *((Fraction(p[0]), Fraction(p[1])) for p in (s1, e1, s2, e2))
  )
  check_eq(kind, bk)
  if kind != 1:
    return 0
  if all(type(c) is int for p in (s1, e1, s2, e2) for c in p):
    check_eq(x, (float(bq[0]), float(bq[1])))  # correctly rounded
  return max(abs(Fraction(x[0]) - bq[0]), abs(Fraction(x[1]) - bq[1]))


# line_inter at 1e9 with the point far from s1 compared with
# |e1 - s1|: the weights of s1 and e1 are huge and opposite. Integer
# input is exact; float input was off by about 3600 (reproducer).
s2, e2 = (-810232897, -796596951), (731982243, 720915984)
for h in (0, 0.5):
  s1, e1 = (999999655 + h, 999999162), (999999658 + h, 999999165)
  check(far_error(s1, e1, s2, e2) < 1e-4)
for it in range(iters(6000)):
  s1, s2, e2 = rnd_pt(C), rnd_pt(C), rnd_pt(C)
  e1 = rnd_pt(C) if it % 2 else (s1[0] + rnd(-10, 10), s1[1] + rnd(-10, 10))
  if it % 3:
    s1, e1, s2, e2 = ((p[0] + 0.5, p[1] - 0.25) for p in (s1, e1, s2, e2))
  if s1 == e1 or s2 == e2:
    continue
  bq = brute_line(*((Fraction(p[0]), Fraction(p[1])) for p in (s1, e1, s2, e2)))
  if bq[0] == 1 and max(abs(bq[1][0]), abs(bq[1][1])) > C:
    continue
  check(far_error(s1, e1, s2, e2) < 1e-12 * C, f"{(s1, e1, s2, e2)}")
check(lib.on_segment((-C, -C), (C, C), (C - 1, C - 1)))
check_eq(lib.side_of((-C, -C), (C, C - 1), (C - 1, C - 2)), -1)

segs = [
  tuple((rng.random(), rng.random()) for _ in range(4)) for _ in range(200_000)
]
with Timer("2e5 seg_inter + seg_dist on floats"):
  total = 0
  for a, b, c, d in segs:
    total += len(lib.seg_inter(a, b, c, d)) + lib.seg_dist(a, b, c)
print("OK geometry/lines")
