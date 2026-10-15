/**
 * Title: Half-plane intersection
 * Description: Intersects the half-planes left of the directed
 *   lines `s -> e` inside the box $[-"big", "big"]^2$ and returns
 *   the vertices of the result counter-clockwise. Sorts the lines
 *   by angle and keeps the useful ones in a deque.
 * Precondition: Floating `T`; `s != e` in each half-plane. Points
 *   within `eps` $= 10^(-9)$ of a line, plus $10^(-15) (|x| + |y|)$
 *   for rounding far out, count as on it: the result is empty when
 *   the intersection is empty or has no interior (a point or a
 *   segment), and has no repeated or collinear vertices. Directions
 *   closer than `eps` count as parallel. Fine for coordinates up to
 *   about $10^4$.
 * Time: $O(n log n)$.
 * Usage:
 *   vector<array<P<double>, 2>> hs; // keep left of s -> e
 *   hs.push_back({s, e});
 *   vector<P<double>> poly = half_plane(hs); // empty if none
 * Tweak: Scale `eps` with larger coordinates; `P<long double>`
 *   with a `long double` box gives more precision. The default box
 *   only keeps unbounded results finite; pass a smaller `big` when
 *   the problem bounds the plane.
 * Example:
 *   // Kernel of a simple polygon (where a guard sees everything):
 *   // the polygon counter-clockwise, each edge keeps its left.
 *   rep(i, 0, n) hs.push_back({poly[i], poly[(i + 1) % n]});
 *   auto ker = half_plane(hs); // area2(ker) / 2 is its area
 *   // Largest circle inside a convex polygon: binary search r,
 *   // shift each edge left by r and test for a non-empty result.
 *   P<double> sh = (e - s).perp().unit() * r;
 *   hs2.push_back({s + sh, e + sh});
 * Notes: Decide whether touching counts. This code drops
 *   results without interior, so to find a point on the closed
 *   side of every line (a single point may be the answer),
 *   first move each line outward by a little more than `eps`. A
 *   vertex on the box means the true intersection is unbounded.
 *   Doubles near the box at $10^9$ are only good to about
 *   $10^(-7)$: without the rounding term in `out`, a line
 *   through a box corner would split it into two vertices.
 *   Linear programs in two variables are half-plane
 *   intersections: maximise the objective over the vertices. "A
 *   point that sees every wall" (polygon kernel), "the region
 *   closer to site $a$ than to all others" (a Voronoi cell:
 *   keep the $a$ side of each perpendicular bisector) and "the
 *   largest circle or offset that fits" are the usual
 *   disguises. Integer input with exact answers needs a
 *   rational version: every test is the sign of a degree-4
 *   polynomial in the coordinates.
 * Source: Zhu Zeyuan's sort-and-deque algorithm; cp-algorithms
 *   Half-plane intersection (idea).
 * Status: stress-tested against an exact brute force
 */
#pragma once
#include "lines.hpp"

template<class T> vector<P<T>>
half_plane(vector<array<P<T>, 2>> hs, T big = 1e9) {
  const T eps = 1e-9;
  P<T> c{-big, -big}; // box corners: c, c.perp(), ...
  rep(i, 0, 4) hs.push_back({c, c.perp()}), c = c.perp();
  auto ang = [](auto& h) { // + 0 turns a y of -0.0 into 0.0
    P<T> d = h[1] - h[0];
    return atan2(d.y + 0, d.x);
  };
  sort(all(hs),
       [&](auto& a, auto& b) { return ang(a) < ang(b); });
  auto dir = [](auto& h) { return (h[1] - h[0]).unit(); };
  auto at = [](auto& a, auto& b) {
    return line_inter(a[0], a[1], b[0], b[1]).second;
  };
  auto out = [&](auto& h, P<T> p) { // eps, + rounding at |p|
    T tol = eps + 1e-15 * (abs(p.x) + abs(p.y));
    return line_dist(h[0], h[1], p) < tol;
  };
  deque<array<P<T>, 2>> dq;
  for (auto& h : hs) {
    while (sz(dq) > 1 && out(h, at(dq.end()[-2], dq.back())))
      dq.pop_back();
    while (sz(dq) > 1 && out(h, at(dq[0], dq[1])))
      dq.pop_front();
    if (sz(dq) && abs(dir(h).cross(dir(dq.back()))) < eps) {
      if (dir(h).dot(dir(dq.back())) < 0) return {};
      if (!out(h, dq.back()[0])) continue;
      dq.pop_back();
    }
    dq.push_back(h);
  }
  while (sz(dq) > 2 && out(dq[0], at(dq.end()[-2], dq.back())))
    dq.pop_back();
  while (sz(dq) > 2 && out(dq.back(), at(dq[0], dq[1])))
    dq.pop_front();
  if (sz(dq) < 3) return {};
  vector<P<T>> res;
  rep(i, 0, sz(dq))
    res.push_back(at(dq[i], dq[(i + 1) % sz(dq)]));
  return res;
}
