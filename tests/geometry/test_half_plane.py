import math
from fractions import Fraction

from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("geometry/half-plane")


def rnd_pt(c):
  return (rnd(-c, c), rnd(-c, c))


def cr(o, a, b):
  return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def box(big):
  c = [(-big, -big), (big, -big), (big, big), (-big, big)]
  return [(c[i - 1], c[i]) for i in range(4)]


def meet(a, b):
  """Exact meeting point of two lines, or None if parallel."""
  (s1, e1), (s2, e2) = a, b
  d1 = (e1[0] - s1[0], e1[1] - s1[1])
  d2 = (e2[0] - s2[0], e2[1] - s2[1])
  den = d1[0] * d2[1] - d1[1] * d2[0]
  if den == 0:
    return None
  t = Fraction((s2[0] - s1[0]) * d2[1] - (s2[1] - s1[1]) * d2[0], den)
  return (s1[0] + d1[0] * t, s1[1] + d1[1] * t)


def exact_region(hs, big):
  """Vertices of the region: the hull of feasible line meetings."""
  hs = hs + box(big)
  pts = set()
  for i in range(len(hs)):
    for j in range(i):
      q = meet(hs[i], hs[j])
      if q is not None and all(cr(s, e, q) >= 0 for s, e in hs):
        pts.add(q)
  pts = sorted(pts)
  if len(pts) < 3:
    return []
  h = []
  for chain in (pts, pts[::-1]):
    s = len(h)
    for p in chain:
      while len(h) >= s + 2 and cr(h[-2], h[-1], p) <= 0:
        h.pop()
      h.append(p)
    h.pop()
  return h if len(h) >= 3 else []


def same_region(got, want, msg):
  """Same vertices, counter-clockwise, up to the starting vertex."""
  check_eq(len(got), len(want), msg)
  if not got:
    return 0
  w0 = want[0]
  k = min(range(len(got)), key=lambda i: math.dist(got[i], w0))
  for i, w in enumerate(want):
    g = got[(k + i) % len(got)]
    check(abs(g[0] - w[0]) < 1e-6 and abs(g[1] - w[1]) < 1e-6, msg)
  return len(got)


def check_region(hs, big, own_box=False):
  got = lib.half_plane(hs + box(big)) if own_box else lib.half_plane(hs, big)
  return same_region(got, exact_region(hs, big), str(hs))


def random_planes(n, c):
  o = rnd_pt(c)  # most half-planes contain o: non-empty results
  hs = []
  while len(hs) < n:
    s, e = rnd_pt(c), rnd_pt(c)
    if s == e:
      continue
    if rnd(0, 4) and cr(s, e, o) < 0:
      s, e = e, s
    hs.append((s, e))
  return hs


nonempty, total = 0, iters(2500)
for it in range(total):
  n, C = rnd(0, 9), rnd(1, 5)
  big = rnd(C, 3 * C)
  hs = random_planes(n, C)
  if it % 5 == 1 and n:  # repeated and reversed half-planes
    for _ in range(3):
      s, e = hs[rnd(0, n - 1)]
      if rnd(0, 1):
        s, e = e, s
      if rnd(0, 1):
        e = (2 * e[0] - s[0], 2 * e[1] - s[1])
      hs.append((s, e))
  if it % 5 == 2:  # all lines through one point: a cone or a point
    o = rnd_pt(C)
    hs = [(s, (s[0] + 1, s[1] + 1) if s == o else o) for s, _ in hs]
  if it % 5 == 3:  # few directions: many parallel lines, strips
    for i, (s, _) in enumerate(hs):
      k = rnd(-2, 2)
      hs[i] = (s, (s[0] + rnd(0, 1) * k, s[1] + rnd(-1, 1) * k))
  hs = [(s, e) for s, e in hs if s != e]
  nonempty += check_region(hs, big, it % 3 == 0) > 0
check(total // 3 < nonempty < total * 9 // 10)

# The default box at 1e9, where doubles resolve only about 1e-7.
# Lines y = x and y = -x pass through its corners; their points,
# integer or decimal, must not split a corner into two vertices.
for it in range(iters(1500)):
  hs = random_planes(rnd(0, 4), 5)
  hd = list(hs)
  for _ in range(rnd(1, 3)):
    a, b, sg = rnd(-5000, 5000), rnd(-5000, 5000), rnd(0, 1) * 2 - 1
    if a == b:
      continue
    k = 1 if a < b else -1
    hs.append(((0, 0), (k, sg * k)))  # the same line
    if it % 2 == 0:  # decimal points on the same line
      a, b = a / 1000, b / 1000
    hd.append(((a, sg * a), (b, sg * b)))
  same_region(lib.half_plane(hd), exact_region(hs, 10**9), str(hd))

o, X, Y = (0, 0), (1, 0), (0, 1)
check_eq(
  lib.half_plane([], 1), [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
)
check_region([(o, X), (Y, (-1, 1))], 3)  # strip 0 <= y <= 1
check_region([(o, X), (o, (-1, 0))], 3)  # the line y = 0: empty
check_region([(o, X), (X, Y), (Y, o), (o, (1, 1))], 3)
check_region([(o, (0, -1)), (o, Y), (o, X), (o, (-1, 0))], 3)  # a point
up = lib.half_plane([((0, 0), (1, 0))])  # y >= 0, default box
check_eq(up, [(-1e9, 0.0), (1e9, 0.0), (1e9, 1e9), (-1e9, 1e9)])
# x + y >= 0 cuts the default box through two corners: a triangle.
check_eq(len(lib.half_plane([((-3, 3), (-1, 1))])), 3)
check_eq(len(lib.half_plane([((3.48, 3.48), (-3.141, -3.141))])), 3)
L = [(0, 0), (4, 0), (4, 1), (1, 1), (1, 4), (0, 4)]
ker = lib.half_plane([(L[i - 1], L[i]) for i in range(6)])
check_eq(ker, [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
# A direction with y = -0.0 (atan2 gives -pi) and a parallel one.
z = [((1.0, 0.0), (0.0, -0.0)), ((1, 1), (0, 1)), ((0, -1), (1, -1))]
check_eq(lib.half_plane(z, 3.0), [(-3, -1), (3, -1), (3, 0), (-3, 0)])

# Large: tangents of a circle at jittered angles (all are edges),
# area from the angles; then random half-planes containing the disc.
n, R = 200_000, 1000.0
ang = [2 * math.pi * (i + rng.uniform(0.1, 0.9)) / n for i in range(n)]
hs = []
for a in ang:
  t = (R * math.cos(a), R * math.sin(a))
  hs.append((t, (t[0] - math.sin(a), t[1] + math.cos(a))))
rng.shuffle(hs)
with Timer("half-plane intersection of 2e5 tangents"):
  poly = lib.half_plane(hs)
check_eq(len(poly), n)
area = sum(
  poly[i - 1][0] * poly[i][1] - poly[i - 1][1] * poly[i][0] for i in range(n)
)
want = sum(
  R * R * math.tan(((ang[(i + 1) % n] - ang[i]) % (2 * math.pi)) / 2)
  for i in range(n)
)
check(abs(area / 2 - want) < 1e-6 * want)
hs = []
for _ in range(n):
  a, r = rng.uniform(0, 2 * math.pi), rng.uniform(R, 3 * R)
  t = (r * math.cos(a), r * math.sin(a))
  hs.append((t, (t[0] - math.sin(a), t[1] + math.cos(a))))
with Timer("half-plane intersection of 2e5 random"):
  poly = lib.half_plane(hs)
k = len(poly)
check(all(cr(poly[i - 2], poly[i - 1], poly[i]) > 0 for i in range(k)))
area = sum(
  poly[i - 1][0] * poly[i][1] - poly[i - 1][1] * poly[i][0] for i in range(k)
)
check(area / 2 > math.pi * R * R)
for s, e in hs[:2000]:
  d = math.hypot(e[0] - s[0], e[1] - s[1])
  check(all(cr(s, e, p) / d > -1e-6 for p in poly))
print("OK geometry/half-plane")
