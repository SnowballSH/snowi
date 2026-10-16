"""
Title: 3D convex hull
Description: Triangles `(a, b, c)` (indices) covering the surface of
  the convex hull, counter-clockwise seen from outside: the normal
  $(p_b - p_a) times (p_c - p_a)$ points out. A face with more than
  three points comes as several triangles.
Precondition: Integer tuples `(x, y, z, ...)`, not all coplanar;
  duplicates are fine. Exact (Python integers), and below $2^63$, so
  fast, for $|x|, |y|, |z| <= 5 dot 10^5$.
Time: $O(n^2)$ time and $n^2$ bytes; 1 s for 2000 points on a
  sphere.
Usage:
  faces = hull3d(pts)
Tweak: Delaunay triangulation of points in the plane, not all on
  one line or circle: `hull3d` of $(x, y, x^2 + y^2)$, keeping the
  faces whose normal has $z < 0$ (clockwise in the plane); fast for
  $|x|, |y| <= 2 dot 10^4$. Its edges contain the Euclidean MST and
  every nearest neighbour; Voronoi vertices are the circumcenters.
Example:
  >>> b = (0, 1)
  >>> cube = [(x, y, z) for x in b for y in b for z in b]
  >>> len(hull3d(cube + [(0, 0, 1), (1, 1, 1)]))  # 2 per square
  12
Notes: Incremental: a new point removes the faces it sees strictly
  and joins the horizon edges to itself. A new face keeps the
  direction of the horizon edge from the removed face, so it is
  oriented correctly even when it is coplanar with its neighbour;
  hence exact arithmetic handles coplanar points. Cocircular points
  in the Delaunay lift give an arbitrary triangulation of their
  polygon. Surface area: half the sum of the normals' lengths.
Source: Incremental construction as in kactl 3dHull (CC0),
  orientation by edge direction.
Status: stress-tested against brute force
"""


def hull3d(p):
  n = len(p)
  X, Y, Z = zip(*(q[:3] for q in p))

  def normal(a, b, c):
    ux, uy, uz = X[b] - X[a], Y[b] - Y[a], Z[b] - Z[a]
    vx, vy, vz = X[c] - X[a], Y[c] - Y[a], Z[c] - Z[a]
    nx, ny = uy * vz - uz * vy, uz * vx - ux * vz
    return nx, ny, ux * vy - uy * vx

  def side(f, i):  # > 0: face f sees point i
    a = f[0]
    dx, dy, dz = X[i] - X[a], Y[i] - Y[a], Z[i] - Z[a]
    return f[3] * dx + f[4] * dy + f[5] * dz

  a, b = 0, 1
  while (X[a], Y[a], Z[a]) == (X[b], Y[b], Z[b]):
    b += 1
  c = b + 1
  while normal(a, b, c) == (0, 0, 0):
    c += 1
  d = c + 1
  while side((a, b, c) + normal(a, b, c), d) == 0:
    d += 1
  if side((a, b, c) + normal(a, b, c), d) > 0:
    b, c = c, b
  faces = []
  E = bytearray(n * n)  # E[x * n + y]: edge x -> y is live

  def add(x, y, z):
    faces.append((x, y, z) + normal(x, y, z))
    E[x * n + y] = E[y * n + z] = E[z * n + x] = 1

  add(a, b, c)
  add(a, d, b)
  add(b, d, c)
  add(c, d, a)
  for i in range(n):
    if i in (a, b, c, d):
      continue
    keep, vis = [], []
    for f in faces:
      (vis if side(f, i) > 0 else keep).append(f)
    if not vis:
      continue
    for x, y, z, *_ in vis:
      E[x * n + y] = E[y * n + z] = E[z * n + x] = 0
    faces[:] = keep
    for f in keep:  # horizon edges: the face across is gone
      for x, y in ((f[0], f[1]), (f[1], f[2]), (f[2], f[0])):
        if not E[y * n + x]:
          add(y, x, i)
  return [f[:3] for f in faces]
