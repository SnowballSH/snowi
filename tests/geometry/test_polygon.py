import math

from testlib import Timer, check, check_eq, iters, load, rnd

lib = load("geometry/polygon")


def rnd_pt(c):
  return (rnd(-c, c), rnd(-c, c))


def close(u, v, tol=1e-7):
  return abs(u - v) <= tol * max(1.0, abs(v))


def exact_on(a, b, p):
  A, B = b[1] - a[1], a[0] - b[0]
  return (
    A * p[0] + B * p[1] == A * a[0] + B * a[1]
    and min(a[0], b[0]) <= p[0] <= max(a[0], b[0])
    and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])
  )


def brute_in(v, p):
  """1 / 0 / -1 by winding number (angle sum), boundary exact."""
  n = len(v)
  if any(exact_on(v[i - 1], v[i], p) for i in range(n)):
    return 0
  w = 0.0
  for i in range(n):
    a, b = v[i - 1], v[i]
    d = math.atan2(b[1] - p[1], b[0] - p[0]) - math.atan2(
      a[1] - p[1], a[0] - p[0]
    )
    while d > math.pi:
      d -= 2 * math.pi
    while d < -math.pi:
      d += 2 * math.pi
    w += d
  return 1 if abs(w) > math.pi else -1


def is_simple(v):
  n = len(v)
  for i in range(n):
    for j in range(i + 1, n):
      r = lib.seg_inter(v[i], v[(i + 1) % n], v[j], v[(j + 1) % n])
      adjacent = j == i + 1 or (i == 0 and j == n - 1)
      if len(r) > (1 if adjacent else 0):
        return False
  return True


def make_ccw(v):
  k = v.index(min(v))
  a, b, c = v[k - 1], v[k], v[(k + 1) % len(v)]
  if lib.cross(lib.sub(b, a), lib.sub(c, b)) < 0:
    v = v[::-1]
  return v


def star_polygon(n, C):
  while True:
    c = rnd_pt(C // 2)
    pts = [p for p in (rnd_pt(C) for _ in range(n)) if p != c]
    w = [(p[0] - c[0], p[1] - c[1], p) for p in pts]
    lib.angle_sort(w)
    v = [q[2] for q in w]
    m = len(v)
    if m >= 3 and all(
      lib.cross(lib.sub(v[i - 1], c), lib.sub(v[i], c)) > 0 for i in range(m)
    ):
      return v


def monotone_polygon(n, C):
  while True:
    pts = sorted(set(rnd_pt(C) for _ in range(n)))
    if len(pts) < 3:
      continue
    lo, hi = [pts[0]], []
    for p in pts[1:-1]:
      (lo if rnd(0, 1) else hi).append(p)
    v = lo + [pts[-1]] + hi[::-1]
    if is_simple(v):
      return make_ccw(v)


def slab_brute(v, cut=None):
  """Area and moments of v (left of cut = (s, e) if given) by exact
  Gauss quadrature over vertical slabs."""
  n = len(v)
  xs = [p[0] for p in v]
  if cut:
    s, e = cut
    xs.append(s[0])
    for i in range(n):
      a, b = v[i - 1], v[i]
      sa, sb = lib.cross3(s, e, a), lib.cross3(s, e, b)
      if sa * sb < 0:
        xs.append(a[0] + (b[0] - a[0]) * sa / (sa - sb))
  xs.sort()
  area = mx = my = 0.0
  for k in range(len(xs) - 1):
    w, m = xs[k + 1] - xs[k], (xs[k] + xs[k + 1]) / 2
    if w < 1e-12:
      continue
    for g in (-1, 1):
      x = m + g * w / (2 * math.sqrt(3))
      if cut and s[0] == e[0] and (x > s[0] if e[1] > s[1] else x < s[0]):
        continue
      ys = sorted(
        a[1] + (b[1] - a[1]) * (x - a[0]) / (b[0] - a[0])
        for a, b in ((v[i - 1], v[i]) for i in range(n))
        if (a[0] < x) != (b[0] < x)
      )
      length = mom = 0.0
      for i in range(0, len(ys) - 1, 2):
        lo, hi = ys[i], ys[i + 1]
        if cut and s[0] != e[0]:
          yl = s[1] + (e[1] - s[1]) * (x - s[0]) / (e[0] - s[0])
          if e[0] > s[0]:
            lo = max(lo, yl)
          else:
            hi = min(hi, yl)
        if hi > lo:
          length += hi - lo
          mom += (hi * hi - lo * lo) / 2
      area += w / 2 * length
      mx += w / 2 * x * length
      my += w / 2 * mom
  return area, mx, my


def check_polygon(v):
  C = max(max(abs(p[0]), abs(p[1])) for p in v)
  inside = boundary = 0
  for x in range(-C - 1, C + 2):
    for y in range(-C - 1, C + 2):
      want = brute_in(v, (x, y))
      check_eq(lib.in_polygon(v, (x, y)), want, f"{v} {(x, y)}")
      inside += want == 1
      boundary += want == 0
  check_eq(lib.area2(v), 2 * inside + boundary - 2)  # Pick
  rv = v[::-1]
  check_eq(lib.area2(rv), -lib.area2(v))
  area, mx, my = slab_brute(v)
  check(close(area, lib.area2(v) / 2, 1e-9))
  g, rg = lib.centroid(v), lib.centroid(rv)
  check(close(g[0], mx / area) and close(g[1], my / area))
  check(close(rg[0], g[0]) and close(rg[1], g[1]))
  for k in range(4):
    s, e = rnd_pt(C + 1), rnd_pt(C + 1)
    if k == 1:
      s = v[rnd(0, len(v) - 1)]
    if k == 2:
      s, e = v[rnd(0, len(v) - 1)], v[rnd(0, len(v) - 1)]
    if k == 3:
      s, e = v[0], v[1]
    if s == e:
      continue
    ca = slab_brute(v, (s, e))[0]
    for poly, sign in ((v, 1), (rv, -1)):
      cut = lib.polygon_cut(poly, s, e)
      check(close(lib.area2(cut) / 2 * sign, ca, 1e-9), f"{poly} {s} {e}")
      for q in cut:
        check(lib.cross3(s, e, q) >= -1e-9)


def check_far(v):
  """v moved by up to 1e9 as floats (like polygon_cut output):
  precision must follow its size, not its position. Cut points there
  are rounded to about 6e-8, which a sliver's centroid amplifies, so
  the cut is compared with the same points moved back to 0."""
  bx, by = rnd(-(10**9), 10**9) + 0.5, rnd(-(10**9), 10**9) - 0.25
  far = [(x + bx, y + by) for x, y in v]
  check_eq(lib.area2(far), lib.area2(v))  # differences are exact
  g, gf = lib.centroid(v), lib.centroid(far)
  check(close(gf[0] - bx, g[0], 1e-6) and close(gf[1] - by, g[1], 1e-6))
  s, e = rnd_pt(8), rnd_pt(8)
  if s == e:
    return
  ca = slab_brute(v, (s, e))[0]
  shift = ((s[0] + bx, s[1] + by), (e[0] + bx, e[1] + by))
  cut = lib.polygon_cut(far, *shift)
  check(close(lib.area2(cut) / 2, ca, 1e-5), f"{v} {s} {e}")
  back = [(x - bx, y - by) for x, y in cut]
  check(close(lib.area2(cut), lib.area2(back), 1e-9))
  if ca > 1e-3:
    g, gb = lib.centroid(cut), lib.centroid(back)
    check(close(g[0] - bx, gb[0], 1e-6) and close(g[1] - by, gb[1], 1e-6))


for it in range(iters(400)):
  n, C = rnd(3, 10), rnd(2, 7)
  v = star_polygon(n, C) if it % 2 else monotone_polygon(n, C)
  check_polygon(v)
  for _ in range(5):
    check_far(v)

# Reproducer (float input): area2 128 instead of 36.
B = 987654321.5
tri = [(B, B + 7), (B + 3, B + 1), (B + 5, B + 9)]
check_eq(lib.area2(tri), 36.0)
g = lib.centroid(tri)
check(abs(g[0] - (B + 8 / 3)) < 1e-6 and abs(g[1] - (B + 17 / 3)) < 1e-6)
check_eq(lib.centroid([(1, 1), (4, 1), (1, 4)]), (2.0, 2.0))
for _ in range(iters(2000)):  # cut a rectangle at 1e9 in half
  ox, oy = rnd(-(10**9), 10**9), rnd(-(10**9), 10**9)
  w, h, mx = rnd(1, 100), rnd(1, 100), ox + 0.5 * rnd(1, 199)
  rect = [(ox, oy), (ox + w, oy), (ox + w, oy + h), (ox, oy + h)]
  left = lib.polygon_cut(rect, (mx, 0), (mx, 1))
  lw = min(mx - ox, w)
  check_eq(lib.area2(left), 2 * lw * h)
  g = lib.centroid(left)
  check(close(g[0] - ox, lw / 2, 1e-6) and close(g[1] - oy, h / 2, 1e-6))

comb = [(0, 0), (6, 0), (6, 4), (5, 4), (5, 1), (4, 1)]
comb += [(4, 4), (3, 4), (3, 1), (2, 1), (2, 4), (0, 4)]
line_pts = [(-3, -3), (0, -3), (3, -3), (3, 0)]
line_pts += [(3, 3), (0, 3), (-3, 3), (-3, 0)]
spiral = [(0, 0), (7, 0), (7, 7), (0, 7), (0, 2), (5, 2), (5, 5), (2, 5)]
spiral += [(2, 4), (4, 4), (4, 3), (1, 3), (1, 6), (6, 6), (6, 1), (0, 1)]
for v in (comb, line_pts, spiral):
  check(is_simple(v))
  check_polygon(v)
check_eq(lib.area2(comb), 2 * 18)  # 6 x 4 minus two 1 x 3 notches
check_eq(lib.centroid(lib.polygon_cut(comb, (1, 0), (1, 1))), (0.5, 2.0))
check_eq(len(lib.polygon_cut(comb, (0, 0), (0, 1))), 2)

check_eq(lib.area2([]), 0)
check_eq(lib.in_polygon([], (0, 0)), -1)
check_eq(lib.in_polygon([(1, 1)], (1, 1)), 0)
check_eq(lib.in_polygon([(1, 1)], (0, 1)), -1)
check_eq(lib.polygon_cut([], (0, 0), (1, 0)), [])

# Extremes: Python ints are exact; triangles with 1e9 coordinates.
C = 10**9
for it in range(iters(3000)):
  tri = [rnd_pt(C) for _ in range(3)]
  if lib.area2(tri) < 0:
    tri[0], tri[1] = tri[1], tri[0]
  if lib.area2(tri) == 0:
    continue
  p = rnd_pt(C)
  if it % 3 == 0:
    p = ((tri[0][0] + tri[1][0]) // 2, (tri[0][1] + tri[1][1]) // 2)
  want = 1
  for i in range(3):
    a, b = tri[i - 1], tri[i]
    c = lib.cross3(a, b, p)
    if c < 0:
      want = -1
    elif c == 0 and want == 1:
      want = 0 if exact_on(a, b, p) else -1
  check_eq(lib.in_polygon(tri, p), want)

n = 200_000
big = []
for i in range(n):
  t, r = 2 * math.pi * i / n, rnd(5 * 10**8, C)
  big.append((round(r * math.cos(t)), round(r * math.sin(t))))
queries = [rnd_pt(C) for _ in range(5)] + big[:5]
with Timer("2e5-gon area2 + centroid + cut + 10 in_polygon"):
  a2 = lib.area2(big)
  g = lib.centroid(big)
  cut = lib.polygon_cut(big, (0, -1), (0, 1))
  res = [lib.in_polygon(big, q) for q in queries]
check(a2 > 0 and abs(g[0]) < 5e7 and abs(g[1]) < 5e7)
check(len(cut) > n // 3)
check_eq(res[5:], [0] * 5)
check_eq(res[:2], [brute_in(big, q) for q in queries[:2]])
print("OK geometry/polygon")
