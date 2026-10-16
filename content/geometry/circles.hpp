/**
 * Title: Circles
 * Description: `circle_line` meets circle `(c, r)` with line `ab`:
 *   no point, one (tangent) or two, in the direction `a -> b`.
 *   `circle_inter` meets two circles: no point, one (touching) or
 *   two, the first left of the line from `a` to `b`. `tangents`
 *   gives the common tangents as pairs of touch points (on circle 1,
 *   on circle 2): the outer ones, the inner ones if `r2` is negated,
 *   and the tangents from point `c2` if `r2 = 0`. `circumcircle`
 *   returns the center and radius of the circle through three
 *   points, `mec` those of the smallest circle containing all
 *   points (Welzl, randomised). `circle_poly` is the signed area of
 *   the circle's intersection with a simple polygon (positive if
 *   counter-clockwise).
 * Precondition: Floating `T`, radii of type `T` too (`5.0`, not
 *   `5`). Lengths within `eps` $= 10^(-9)$ count as equal (fine for
 *   coordinates up to about $10^4$). Concentric circles give no
 *   point and no tangent. `circumcircle`: points not collinear.
 *   `mec`: at least one point. `circle_poly`: no two consecutive
 *   vertices equal.
 * Time: $O(1)$; `mec` expected $O(n)$, `circle_poly` $O(n)$.
 * Usage:
 *   vector<P<double>> v = circle_line(c, r, a, b);
 *   vector<P<double>> w = circle_inter(c1, r1, c2, r2);
 *   auto t = tangents(c1, r1, c2, r2); // t[i].first on circle 1
 *   auto [o, rad] = circumcircle(a, b, c);
 *   auto [center, radius] = mec(pts);
 *   double area = circle_poly(c, r, poly);
 * Tweak: For larger coordinates scale `eps` (the noise is about
 *   $10^(-13)$ times the largest coordinate); `mec` has no `eps`.
 *   More precision: `P<long double>` with `long double` radii.
 * Example:
 *   // Touch points of the tangents from p to circle (c, r):
 *   for (auto [t, q] : tangents(c, r, p, 0.0))
 *     printf("%.9f %.9f\n", t.x, t.y);
 *   // Inner tangents (they cross between the circles):
 *   auto inner = tangents(c1, r1, c2, -r2);
 *   // Radius of the circle through p, q, s (not collinear):
 *   double R = circumcircle(p, q, s).second;
 * Notes: `mec` works relative to its first point, so its precision
 *   follows the spread of the points, not their distance from the
 *   origin; it counts points up to `r * 1e-9` outside as inside.
 *   With integer input, decide the case exactly in integers
 *   (circles touch iff $d^2 = (r_1 plus.minus r_2)^2$; line `ab` is
 *   tangent iff $((b - a) times (c - a))^2 = r^2 |b - a|^2$) and
 *   use these functions only for the points. The minimum enclosing
 *   circle hides behind "place one station to minimise the farthest
 *   distance". For the largest circle inside a convex polygon,
 *   binary search the radius, push every edge inward by it and test
 *   whether the half-plane intersection is empty. Circumradius from
 *   side lengths $x, y, z$ and area $S$: $R = x y z / (4 S)$.
 * Source: kactl CircleLine, CircleIntersection, CircleTangents,
 *   Circumcircle, MinimumEnclosingCircle, CirclePolygonIntersection
 *   (CC0).
 * Status: stress-tested against brute force
 */
#pragma once
#include "point.hpp"

template<class T>
vector<P<T>> circle_line(P<T> c, T r, P<T> a, P<T> b) {
  const T eps = 1e-9;
  P<T> ab = b - a, p = a + ab * ((c - a).dot(ab) / ab.dist2());
  T d = (c - p).dist();
  if (d > r + eps) return {};
  if (d > r - eps) return {p};
  P<T> h = ab.unit() * sqrt(r * r - d * d);
  return {p - h, p + h};
}
template<class T>
vector<P<T>> circle_inter(P<T> a, T r1, P<T> b, T r2) {
  const T eps = 1e-9;
  P<T> v = b - a;
  T d = v.dist(), hi = r1 + r2, lo = abs(r1 - r2);
  if (d < eps || d > hi + eps || d < lo - eps) return {};
  P<T> m = a + v * ((d * d + r1 * r1 - r2 * r2) / (2 * d * d));
  if (d > hi - eps || d < lo + eps) return {m};
  T h = sqrt((d - lo) * (d + lo) * (hi - d) * (hi + d));
  P<T> n = v.perp() * (h / (2 * d * d));
  return {m + n, m - n};
}
template<class T> vector<pair<P<T>, P<T>>>
tangents(P<T> c1, T r1, P<T> c2, T r2) {
  const T eps = 1e-9;
  P<T> d = c2 - c1;
  T dr = r1 - r2, d2 = d.dist2(), len = sqrt(d2);
  if (len < eps || len < abs(dr) - eps) return {};
  T h = len < abs(dr) + eps ? 0 : sqrt(d2 - dr * dr);
  vector<pair<P<T>, P<T>>> out;
  for (int s : {-1, 1}) {
    P<T> v = (d * dr + d.perp() * h * s) / d2;
    out.push_back({c1 + v * r1, c2 + v * r2});
  }
  if (h == 0) out.pop_back();
  return out;
}
template<class T>
pair<P<T>, T> circumcircle(P<T> a, P<T> b, P<T> c) {
  P<T> u = c - a, v = b - a, w = u * v.dist2() - v * u.dist2();
  P<T> o = a + w.perp() / u.cross(v) / 2;
  return {o, (o - a).dist()};
}
template<class T> pair<P<T>, T> mec(vector<P<T>> ps) {
  shuffle(all(ps), mt19937(time(0)));
  P<T> base = ps[0], o{};
  for (P<T>& p : ps) p = p - base;
  T r = 0, slack = 1 + 1e-9;
  auto out = [&](P<T> p) { return (o - p).dist() > r * slack; };
  rep(i, 0, sz(ps))
    if (out(ps[i])) {
      o = ps[i], r = 0;
      rep(j, 0, i)
        if (out(ps[j])) {
          o = (ps[i] + ps[j]) / 2, r = (o - ps[i]).dist();
          rep(k, 0, j)
            if (out(ps[k]))
              tie(o, r) = circumcircle(ps[i], ps[j], ps[k]);
        }
    }
  return {base + o, r};
}
template<class T>
T circle_poly(P<T> c, T r, const vector<P<T>>& poly) {
  auto arc = [&](P<T> p, P<T> q) { // sector of the turn p -> q
    return atan2(p.cross(q), p.dot(q)) * r * r / 2;
  };
  T sum = 0;
  rep(i, 0, sz(poly)) {
    P<T> p = poly[i] - c, q = poly[(i + 1) % sz(poly)] - c;
    P<T> d = q - p;
    T a = d.dot(p) / d.dist2(), s = 1, t = 0;
    T h = a * a - (p.dist2() - r * r) / d.dist2();
    if (h > 0) // the chord p + d [s, t] lies inside
      h = sqrt(h), s = max<T>(0, -a - h), t = min<T>(1, -a + h);
    if (t < 0 || s >= 1) sum += arc(p, q);
    else {
      P<T> u = p + d * s, w = p + d * t;
      sum += arc(p, u) + u.cross(w) / 2 + arc(w, q);
    }
  }
  return sum;
}
