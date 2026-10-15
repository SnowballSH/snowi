/**
 * Title: Lines and segments
 * Description: `side_of(s, e, p)` is 1 if `p` is left of the
 *   directed line `s -> e`, 0 on it and -1 right of it.
 *   `on_segment` tests the closed segment. `line_inter` returns
 *   `{1, point}`, or `{0, _}` for parallel lines and `{-1, _}` for
 *   the same line. `seg_inter` intersects two closed segments: no
 *   point, one point, or the two ends of a collinear overlap,
 *   sorted. `line_dist` is the signed distance from `p` to line
 *   `ab` (positive on the left), `seg_dist` the distance to segment
 *   `ab`. `line_proj` projects `p` onto line `ab`, or reflects it
 *   across the line when `refl` is true.
 * Precondition: The two points defining a line differ. `side_of`,
 *   `on_segment` and the case taken by `seg_inter` and `seg_dist`
 *   are exact for `P<ll>` with $|x|, |y| <= 10^9$. Computed points
 *   need floating `T`; `P<long double>` keeps those cases exact for
 *   such integer input (it holds integers below $2^64$).
 * Time: $O(1)$.
 * Usage:
 *   int s = side_of(a, b, p); bool on = on_segment(a, b, p);
 *   auto [kind, x] = line_inter(a, b, c, d);
 *   vector<P<double>> v = seg_inter(a, b, c, d); // sz(v) <= 2
 *   double d1 = seg_dist(a, b, p), d2 = line_dist(a, b, p);
 *   P<double> q = line_proj(a, b, p), r = line_proj(a, b, p, 1);
 * Tweak: Exact yes/no for crossing `P<ll>` segments: in the first
 *   branch of `seg_inter` return `{a}` instead of the point. For
 *   non-integer doubles, the `== 0` in `on_segment` should become a
 *   tolerance: `seg_dist(s, e, p) < eps`.
 * Example:
 *   // Where do segments ab and cd meet?
 *   vector<P<double>> v = seg_inter(a, b, c, d);
 *   if (v.empty()) puts("none");
 *   else if (sz(v) == 1) printf("%.9f %.9f\n", v[0].x, v[0].y);
 *   else puts("they overlap from v[0] to v[1]");
 *   // Are p and q strictly on the same side of line ab?
 *   bool same = side_of(a, b, p) * side_of(a, b, q) > 0;
 *   // Mirror image of p in line ab (light reflection):
 *   P<double> img = line_proj(a, b, p, true);
 * Notes: Two segments cross properly iff the ends of each lie
 *   strictly on opposite sides of the other (`side_of` products are
 *   -1); everything else is a touching or collinear case, which
 *   `on_segment` decides. Read whether touching counts. "Mirror",
 *   "reflection", "unfold the path" call for `line_proj(.., true)`;
 *   "distance to a wall or road" is `seg_dist`. With `P<ll>`,
 *   `line_inter` would truncate (and may overflow) its point:
 *   convert to `P<long double>` first, which keeps the cross
 *   products of integers up to $10^9$ exact. `P<double>` rounds
 *   them, so the point may be off by about $10^(-5)$ there.
 * Source: kactl SegmentIntersection, LineIntersection,
 *   SegmentDistance, LineProjectionReflection (CC0).
 * Status: stress-tested against brute force
 */
#pragma once
#include "point.hpp"

template<class T> int side_of(P<T> s, P<T> e, P<T> p) {
  return sgn(s.cross(e, p));
}
template<class T> bool on_segment(P<T> s, P<T> e, P<T> p) {
  return p.cross(s, e) == 0 && (s - p).dot(e - p) <= 0;
}
template<class T>
pair<int, P<T>> line_inter(P<T> s1, P<T> e1, P<T> s2, P<T> e2) {
  T d = (e1 - s1).cross(e2 - s2);
  if (d == 0) return {-(s1.cross(e1, s2) == 0), {}};
  T q = s2.cross(e2, s1);
  return {1, (s1 * d + (e1 - s1) * q) / d};
}
template<class T>
vector<P<T>> seg_inter(P<T> a, P<T> b, P<T> c, P<T> d) {
  T oa = c.cross(d, a), ob = c.cross(d, b), oc = a.cross(b, c),
    od = a.cross(b, d);
  if (sgn(oa) * sgn(ob) < 0 && sgn(oc) * sgn(od) < 0)
    return {(a * ob - b * oa) / (ob - oa)};
  set<P<T>> s;
  if (on_segment(c, d, a)) s.insert(a);
  if (on_segment(c, d, b)) s.insert(b);
  if (on_segment(a, b, c)) s.insert(c);
  if (on_segment(a, b, d)) s.insert(d);
  return {all(s)};
}
template<class T> auto line_dist(P<T> a, P<T> b, P<T> p) {
  return a.cross(b, p) / (b - a).dist();
}
template<class T> auto seg_dist(P<T> a, P<T> b, P<T> p) {
  if ((p - a).dot(b - a) <= 0) return (p - a).dist();
  if ((p - b).dot(a - b) <= 0) return (p - b).dist();
  return abs(line_dist(a, b, p));
}
template<class T>
P<T> line_proj(P<T> a, P<T> b, P<T> p, bool refl = false) {
  P<T> v = b - a;
  return p - v.perp() * (1 + refl) * a.cross(b, p) / v.dist2();
}
