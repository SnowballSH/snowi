/**
 * Title: Polygon
 * Description: `area2` is twice the signed area, positive when the
 *   vertices go counter-clockwise. `centroid` is the center of mass
 *   of the enclosed area. `in_polygon` returns 1 if `p` is strictly
 *   inside, 0 if it is on the boundary and -1 if outside.
 *   `polygon_cut` keeps the part of the polygon on the left of the
 *   directed line `s -> e` (the line itself included).
 * Precondition: Vertices in order, either orientation, first not
 *   repeated at the end. `area2` and `in_polygon` are exact for
 *   `P<ll>`. `centroid` (nonzero area) and `polygon_cut` need
 *   floating `T`. Self-intersecting polygons: `in_polygon` uses the
 *   even-odd rule.
 * Time: $O(n)$.
 * Usage:
 *   ll a2 = area2(poly); // area = abs(a2) / 2.0
 *   int where = in_polygon(poly, p); // 1 in, 0 edge, -1 out
 *   P<double> g = centroid(pd);
 *   vector<P<double>> left = polygon_cut(pd, s, e);
 * Tweak: Coordinates near $10^9$ with many vertices: accumulate
 *   `area2` in `__int128`. For non-integer doubles, make the boundary
 *   test in `in_polygon` a tolerance: `seg_dist(a, b, p) < eps`.
 * Example:
 *   // Area of a convex polygon inside the strip between two
 *   // parallel lines: cut twice, ends swapped for the other side.
 *   auto mid = polygon_cut(polygon_cut(pd, s1, e1), e2, s2);
 *   double area = area2(mid) / 2;
 *   // Lattice points strictly inside (Pick), integer poly:
 *   ll b = 0;
 *   rep(i, 0, n) {
 *     P<ll> d = poly[(i + 1) % n] - poly[i];
 *     b += gcd(abs(d.x), abs(d.y));
 *   }
 *   ll inside = (abs(area2(poly)) - b + 2) / 2;
 * Notes: Pick's theorem: a polygon with lattice vertices has area
 *   $A = I + B/2 - 1$, where $I$ counts the lattice points strictly
 *   inside and $B = sum gcd(|Delta x|, |Delta y|)$ over the edges
 *   counts those on the boundary. Input given clockwise flips the
 *   sign of `area2`; check it before using code that assumes
 *   counter-clockwise order. `polygon_cut` of a non-convex polygon
 *   can return several pieces joined by zero-width bridges along
 *   the line; its area is still right. It is the step of
 *   Sutherland-Hodgman clipping: cutting by every edge of a convex
 *   polygon intersects with it in $O(n m)$. `area2` and `centroid`
 *   sum relative to `v[0]`, so a small polygon far from the origin
 *   keeps its precision in `P<double>`.
 * Source: kactl PolygonArea, PolygonCenter, InsidePolygon,
 *   PolygonCut (CC0).
 * Status: stress-tested against brute force
 */
#pragma once
#include "lines.hpp"

template<class T> T area2(const vector<P<T>>& v) {
  T a = 0;
  rep(i, 1, sz(v) - 1) a += v[0].cross(v[i], v[i + 1]);
  return a;
}
template<class T> P<T> centroid(const vector<P<T>>& v) {
  P<T> c{};
  T a = 0;
  rep(i, 1, sz(v) - 1) {
    P<T> p = v[i] - v[0], q = v[i + 1] - v[0];
    c = c + (p + q) * p.cross(q), a += p.cross(q);
  }
  return v[0] + c / a / 3;
}
template<class T>
int in_polygon(const vector<P<T>>& v, P<T> p) {
  int in = 0;
  rep(i, 0, sz(v)) {
    P<T> a = v[i], b = v[(i + 1) % sz(v)];
    if (on_segment(a, b, p)) return 0;
    in ^= ((p.y < a.y) - (p.y < b.y)) * p.cross(a, b) > 0;
  }
  return in ? 1 : -1;
}
template<class T> vector<P<T>>
polygon_cut(const vector<P<T>>& v, P<T> s, P<T> e) {
  vector<P<T>> res;
  rep(i, 0, sz(v)) {
    P<T> cur = v[i], prev = v[(i + sz(v) - 1) % sz(v)];
    T a = s.cross(e, cur), b = s.cross(e, prev);
    if (sgn(a) * sgn(b) < 0)
      res.push_back(cur + (prev - cur) * (a / (a - b)));
    if (a >= 0) res.push_back(cur);
  }
  return res;
}
