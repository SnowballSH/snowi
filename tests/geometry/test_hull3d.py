import math

from testlib import Timer, check, check_eq, iters, load, rnd, rng

lib = load("geometry/hull3d")


def sub(u, v):
  return (u[0] - v[0], u[1] - v[1], u[2] - v[2])


def cross(u, v):
  return (
    u[1] * v[2] - u[2] * v[1],
    u[2] * v[0] - u[0] * v[2],
    u[0] * v[1] - u[1] * v[0],
  )


def dot(u, v):
  return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def check_hull(p, faces):
  """Outward faces with all points behind, each directed edge once
  and its reverse once, V - E + F = 2. Returns 6 times the volume."""
  check(len(faces) >= 4)
  edges, verts, vol6 = {}, set(), 0
  for a, b, c in faces:
    check(len({a, b, c}) == 3)
    q = cross(sub(p[b], p[a]), sub(p[c], p[a]))
    check(q != (0, 0, 0))
    check(all(dot(q, sub(r, p[a])) <= 0 for r in p))
    for e in ((a, b), (b, c), (c, a)):
      edges[e] = edges.get(e, 0) + 1
    verts |= {a, b, c}
    vol6 += dot(p[a], cross(p[b], p[c]))
  check(all(k == 1 and (y, x) in edges for (x, y), k in edges.items()))
  check_eq(len(verts) - len(edges) // 2 + len(faces), 2)
  check(vol6 > 0)
  return vol6


def brute_vol6(p):
  """Hull of points in general position: the facets are the triples
  with all other points strictly on one side."""
  n, vol6 = len(p), 0.0
  for i in range(n):
    for j in range(i + 1, n):
      for k in range(j + 1, n):
        q = cross(sub(p[j], p[i]), sub(p[k], p[i]))
        s = [dot(q, sub(p[m], p[i])) for m in range(n) if m not in (i, j, k)]
        t = dot(p[i], cross(p[j], p[k]))
        if all(x < 0 for x in s):
          vol6 += t
        if all(x > 0 for x in s):
          vol6 -= t
  return vol6


def coplanar(p):
  for a in p:
    for b in p:
      for c in p:
        q = cross(sub(b, a), sub(c, a))
        if any(dot(q, sub(d, a)) for d in p):
          return False
  return True


def cross2(o, a, b):
  return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def hull2_area2(pts):
  """Twice the area of the 2D convex hull, by monotone chain."""
  h, st = sorted(set(pts)), []
  for _ in range(2):
    base = len(st)
    for q in h:
      while len(st) >= base + 2 and cross2(st[-2], st[-1], q) <= 0:
        st.pop()
      st.append(q)
    st.pop()
    h.reverse()
  return sum(cross2((0, 0), st[i - 1], st[i]) for i in range(len(st)))


def circle_det(a, b, c, r):
  """Negative iff r is strictly inside the circumcircle of the
  clockwise triangle a, b, c."""
  m = []
  for q in (a, b, c):
    x, y = q[0] - r[0], q[1] - r[1]
    m.append((x, y, x * x + y * y))
  return (
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
    - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
    + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
  )


cube = [(x, y, z) for x in range(3) for y in range(3) for z in range(3)]
rng.shuffle(cube)
check_eq(check_hull(cube, lib.hull3d(cube)), 48)
for it in range(iters(1500)):
  n, c = rnd(4, 12 if it % 2 else 30), (1, 2, 5, 100)[it % 4]
  p = [(rnd(-c, c), rnd(-c, c), rnd(-c, c)) for _ in range(n)]
  if it % 5 == 0:
    p += p[:2]
  if n <= 12 and coplanar(p):
    continue
  vol6 = check_hull(p, lib.hull3d(p))
  if n <= 12:  # perturbed copies are in general position
    e = 1e-7
    pd = [tuple(x + rng.uniform(-e, e) for x in q) for q in p]
    check(abs(brute_vol6(pd) - vol6) < 1e-3 * max(1, c * c))
for _ in range(iters(30)):  # largest coordinates: below 2^63
  c = 500_000
  p = [(rnd(-c, c), rnd(-c, c), rnd(-c, c)) for _ in range(100)]
  check_hull(p, lib.hull3d(p))

# Delaunay by lifting: lower faces, empty circumcircles, and the
# triangles tile the 2D hull.
for it in range(iters(1500)):
  n, c = rnd(3, 25), (2, 4, 20_000)[it % 3]
  p = []
  for _ in range(n):
    x, y = rnd(-c, c), rnd(-c, c)
    p.append((x, y, x * x + y * y))
  if coplanar(p):  # all on a line or circle
    continue
  area2 = 0
  for a, b, cc in lib.hull3d(p):
    q = cross(sub(p[b], p[a]), sub(p[cc], p[a]))
    if q[2] >= 0:
      continue
    area2 -= q[2]  # clockwise in the plane
    check(all(circle_det(p[a], p[b], p[cc], r) >= 0 for r in p))
  check_eq(area2, hull2_area2([q[:2] for q in p]))

n = 2000
p = []
for _ in range(n):  # near a sphere: most points are on the hull
  t, z = rng.uniform(0, 2 * math.pi), rng.uniform(-1, 1)
  r = math.sqrt(1 - z * z)
  p.append(
    (int(1e5 * r * math.cos(t)), int(1e5 * r * math.sin(t)), int(1e5 * z))
  )
with Timer("hull3d of 2000 points near a sphere"):
  faces = lib.hull3d(p)
check_hull(p, faces)
print("OK geometry/hull3d")
