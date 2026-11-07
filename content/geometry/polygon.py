"""
Title: Polygon
Description: `area2` is twice the signed area, positive when the
  vertices go counter-clockwise. `centroid` is the center of mass of
  the enclosed area. `in_polygon` returns 1 if `p` is strictly
  inside, 0 if it is on the boundary and -1 if outside.
  `polygon_cut` keeps the part of the polygon on the left of the
  directed line `s -> e` (the line itself included).
Precondition: Vertices are tuples in order, either orientation,
  first not repeated at the end. Integer input keeps `area2` and
  `in_polygon` exact; `centroid` (nonzero area) is exact up to its
  final division. Self-intersecting polygons: `in_polygon` uses the
  even-odd rule.
Time: $O(n)$.
Usage:
  a2 = area2(poly)  # area = abs(a2) / 2
  where = in_polygon(poly, p)  # 1 in, 0 edge, -1 out
  gx, gy = centroid(poly)
  left = polygon_cut(poly, s, e)
Tweak: For non-integer floats, make the boundary test in
  `in_polygon` a tolerance: `seg_dist(a, b, p) < eps`. Exact cut
  points: `Fraction` coordinates.
Example:
  >>> sq = [(0, 0), (4, 0), (4, 4), (0, 4)]
  >>> area2(sq), centroid(sq), in_polygon(sq, (4, 2))
  (32, (2.0, 2.0), 0)
  >>> polygon_cut(sq, (1, 0), (3, 4))  # left of the line
  [(0, 0), (1.0, 0.0), (3.0, 4.0), (0, 4)]
  >>> b = sum(math.gcd(p[0] - q[0], p[1] - q[1])
  ...         for p, q in zip(sq, sq[1:] + sq[:1]))
  >>> (area2(sq) - b + 2) // 2  # lattice points inside (Pick)
  9
Notes: Pick's theorem: a simple lattice polygon has area
  $A = I + B/2 - 1$, where $I$ counts the lattice points strictly
  inside and $B = sum gcd(|Delta x|, |Delta y|)$ over the edges
  counts those on the boundary. Input given clockwise flips the sign
  of `area2`; check it before using code that assumes
  counter-clockwise order. `polygon_cut` of a non-convex polygon can
  return several pieces joined by zero-width bridges along the line;
  its area is still right. It is the step of Sutherland-Hodgman
  clipping: cutting by every edge of a convex polygon intersects
  with it in $O(n m)$. `area2` and `centroid` sum relative to
  `v[0]`, so a small float polygon far from the origin keeps its
  precision.
Uses: geometry/lines
Source: kactl PolygonArea, PolygonCenter, InsidePolygon,
  PolygonCut (CC0).
Status: stress-tested against brute force
"""


def area2(v):
  a = 0
  for i in range(1, len(v) - 1):
    a += cross3(v[0], v[i], v[i + 1])
  return a


def centroid(v):
  ox, oy = v[0][0], v[0][1]
  cx = cy = a = 0
  for i in range(1, len(v) - 1):
    px, py = v[i][0] - ox, v[i][1] - oy
    qx, qy = v[i + 1][0] - ox, v[i + 1][1] - oy
    c = px * qy - py * qx
    cx += (px + qx) * c
    cy += (py + qy) * c
    a += 3 * c
  return ((ox * a + cx) / a, (oy * a + cy) / a)


def in_polygon(v, p):
  inside, py = False, p[1]
  for i in range(len(v)):
    a, b = v[i - 1], v[i]
    if on_segment(a, b, p):
      return 0
    if ((py < a[1]) - (py < b[1])) * cross3(p, a, b) > 0:
      inside = not inside
  return 1 if inside else -1


def polygon_cut(v, s, e):
  res = []
  for i in range(len(v)):
    prev, cur = v[i - 1], v[i]
    a, b = cross3(s, e, cur), cross3(s, e, prev)
    if a * b < 0:
      t = a / (a - b)
      x = cur[0] + (prev[0] - cur[0]) * t
      res.append((x, cur[1] + (prev[1] - cur[1]) * t))
    if a >= 0:
      res.append(cur)
  return res
