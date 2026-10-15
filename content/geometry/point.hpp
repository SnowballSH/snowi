/**
 * Title: Point and vector
 * Description: 2D point or vector with coordinates of type `T`.
 *   `a.cross(b)` is positive when `b` is counter-clockwise of `a`,
 *   and `o.cross(a, b)` is positive when `o -> a -> b` turns left.
 *   `angle_less` orders vectors by polar angle, starting at the
 *   positive x-axis and going counter-clockwise; it is exact for
 *   integers. `sgn` is the sign as -1, 0 or 1.
 * Precondition: `P<ll>` is exact while products of two coordinate
 *   differences fit: $|x|, |y| <= 10^9$ gives cross products up to
 *   $8 dot 10^18$. `unit` and `rotate` need floating `T`.
 *   `angle_less`: no zero vectors; equal directions tie.
 * Time: $O(1)$.
 * Usage:
 *   P<ll> o{}, a{3, 1}, b(x, y); // (x, y) also converts types
 *   ll c = a.cross(b), turn = o.cross(a, b), d2 = (a - b).dist2();
 *   sort(all(vecs), angle_less<ll>);
 * Tweak: Use `P<double>` or `P<long double>` for real coordinates
 *   (`dist` and `angle` then return that type). To compare sums or
 *   products of cross products with coordinates near $10^9$, use
 *   `P<__int128>`.
 * Example:
 *   // Sort points around c by angle, ties by distance:
 *   sort(all(v), [&](P<ll> a, P<ll> b) {
 *     P<ll> u = a - c, w = b - c;
 *     if (angle_less(u, w) || angle_less(w, u))
 *       return angle_less(u, w);
 *     return u.dist2() < w.dist2();
 *   });
 *   // Same direction: u.cross(w) == 0 && u.dot(w) > 0
 * Notes: The sign of a cross product answers "left or right turn?",
 *   "which side of a line?" and "clockwise or not?"; with integer
 *   input it never needs an epsilon. Compare lengths with `dist2`
 *   to stay in integers. A cross product of differences of
 *   coordinates up to $10^9$ is just below the `ll` limit, so the
 *   sum of two of them overflows. Never sort by `atan2` when
 *   coordinates exceed about $10^7$: distinct directions can round
 *   to the same double. `angle` is in $(-pi, pi]$, except that a
 *   floating `y` of `-0.0` (as in `P<double>{1, 0} * -1.0`) gives
 *   $-pi$; `half` still treats that vector as angle $pi$. `perp`
 *   turns a quarter turn counter-clockwise; `half` is true for
 *   angles in $[pi, 2 pi)$.
 * Source: kactl Point (CC0); half-plane angle order is folklore.
 * Status: stress-tested against brute force
 */
#pragma once

template<class T> int sgn(T x) { return (x > 0) - (x < 0); }
template<class T> struct P {
  T x = 0, y = 0;
  P operator+(P p) const { return {x + p.x, y + p.y}; }
  P operator-(P p) const { return {x - p.x, y - p.y}; }
  P operator*(T d) const { return {x * d, y * d}; }
  P operator/(T d) const { return {x / d, y / d}; }
  auto operator<=>(const P&) const = default; // by x, then y
  T dot(P p) const { return x * p.x + y * p.y; }
  T cross(P p) const { return x * p.y - y * p.x; }
  T cross(P a, P b) const {
    return (a - *this).cross(b - *this);
  }
  T dist2() const { return x * x + y * y; }
  auto dist() const { return sqrt(dist2()); }
  auto angle() const { return atan2(y, x); }
  P unit() const { return *this / dist(); }
  P perp() const { return {-y, x}; }
  P rotate(T a) const {
    return {x * cos(a) - y * sin(a), x * sin(a) + y * cos(a)};
  }
  bool half() const { return y < 0 || (y == 0 && x < 0); }
};
template<class T> bool angle_less(P<T> a, P<T> b) {
  if (a.half() != b.half()) return b.half();
  return a.cross(b) > 0;
}
