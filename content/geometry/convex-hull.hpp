/**
 * Title: Convex hull
 * Description: `convex_hull` returns the vertices of the convex hull
 *   counter-clockwise, starting from the smallest point (by `x`,
 *   then `y`), without collinear points (Andrew's monotone chain).
 *   `in_hull` returns 1 if `p` is strictly inside such a hull, 0 if
 *   on its boundary and -1 if outside. `hull_diameter` returns the
 *   farthest pair of points (rotating calipers).
 * Precondition: `in_hull` and `hull_diameter` take a hull as
 *   returned by `convex_hull` (non-empty for `hull_diameter`). All
 *   exact for `P<ll>` with $|x|, |y| <= 10^9$.
 * Time: $O(n log n)$ for the hull, $O(log n)$ per `in_hull` query,
 *   $O(n)$ for `hull_diameter`.
 * Usage:
 *   vector<P<ll>> h = convex_hull(pts); // 1 or 2 if degenerate
 *   int where = in_hull(h, p); // 1 in, 0 edge, -1 out
 *   auto [a, b] = hull_diameter(h); ll d2 = (a - b).dist2();
 * Tweak: To keep points in the middle of hull edges, pop on `< 0`
 *   instead of `<= 0`; if all points are collinear the middle ones
 *   then appear twice, so return the sorted `pts` in that case.
 *   `in_hull` and `hull_diameter` need the hull without them.
 * Example:
 *   // Fence around all trees: length, area, membership queries.
 *   vector<P<ll>> h = convex_hull(pts);
 *   double len = 0;
 *   rep(i, 0, sz(h)) len += (h[(i + 1) % sz(h)] - h[i]).dist();
 *   ll a2 = area2(h); // from geometry/polygon
 *   for (P<ll> q : qs) puts(in_hull(h, q) >= 0 ? "in" : "out");
 *   // Minimum width: for each hull edge, the farthest vertex from
 *   // its line; that vertex only moves forward, as in the code.
 * Notes: The answer only ever involves hull vertices when the
 *   objective is a maximum of a linear function (max of
 *   $a x + b y$ over points), a farthest pair, a smallest enclosing
 *   shape or "mixtures" (any weighted average of options). Points
 *   given as subarray or window data become hull points after a
 *   rewrite: "add $x$ to every element" turns each candidate into a
 *   line, and the upper envelope of lines is a hull of points
 *   (slope, value). Decide whether collinear points belong to the
 *   hull and whether touching counts before coding. Points that
 *   arrive sorted by `x`, then `y`, without duplicates, can skip
 *   the sort (each chain is $O(n)$); sorted by `x` alone is not
 *   enough.
 * Source: kactl ConvexHull, PointInsideHull, HullDiameter (CC0);
 *   UCup S4-S5 Training Packet 20.
 * Status: stress-tested against brute force
 */
#pragma once
#include "lines.hpp"

template<class T> vector<P<T>> convex_hull(vector<P<T>> pts) {
  sort(all(pts)), pts.erase(unique(all(pts)), pts.end());
  if (sz(pts) <= 2) return pts;
  vector<P<T>> h(2 * sz(pts));
  int k = 0;
  rep(it, 0, 2) {
    for (int s = k; P<T> p : pts) {
      while (k >= s + 2 && h[k - 2].cross(h[k - 1], p) <= 0)
        k--;
      h[k++] = p;
    }
    k--, reverse(all(pts));
  }
  return {h.begin(), h.begin() + k};
}
template<class T> int in_hull(const vector<P<T>>& h, P<T> p) {
  int n = sz(h), lo = 1, hi = n - 1;
  if (n < 3) return n && on_segment(h[0], h[n - 1], p) ? 0 : -1;
  T a = h[0].cross(h[1], p), b = h[0].cross(h[n - 1], p);
  if (a < 0 || b > 0) return -1;
  while (hi - lo > 1) {
    int mid = (lo + hi) / 2;
    (h[0].cross(h[mid], p) >= 0 ? lo : hi) = mid;
  }
  T c = h[lo].cross(h[hi], p);
  return c < 0 ? -1 : c == 0 || a == 0 || b == 0 ? 0 : 1;
}
template<class T>
array<P<T>, 2> hull_diameter(const vector<P<T>>& h) {
  int n = sz(h);
  pair<T, array<P<T>, 2>> res{0, {h[0], h[0]}};
  for (int i = 0, j = 1 % n; i < n; i++) {
    P<T> a = h[i], b = h[(i + 1) % n];
    while ((b - a).cross(h[(j + 1) % n] - h[j]) > 0)
      j = (j + 1) % n;
    for (P<T> p : {a, b})
      res = max(res, {(p - h[j]).dist2(), {p, h[j]}});
  }
  return res.second;
}
