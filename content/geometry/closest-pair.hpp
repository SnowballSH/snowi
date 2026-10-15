/**
 * Title: Closest pair of points
 * Description: Returns two points of `v` at the smallest distance.
 *   Sweeps the points by `y`, keeping those less than the best
 *   distance below the sweep line in a set ordered by `x`.
 * Precondition: At least two points; equal points are allowed
 *   (distance 0). Exact for `P<ll>` with $|x|, |y| <= 10^9$:
 *   squared distances stay below $8 dot 10^18$.
 * Time: $O(n log n)$.
 * Usage:
 *   auto [a, b] = closest_pair(pts); ll d2 = (a - b).dist2();
 * Tweak: `P<double>` works as is. The window is `sqrt(best)`, plus
 *   1 for integer `T`, where truncation would make it too short.
 *   For indices, find `a` in the input, then `b` at another index
 *   (the two may be equal points).
 * Example:
 *   // Tricky function: min over i != j of (i - j)^2 + (S_i - S_j)^2
 *   // with prefix sums S is a closest pair of the points (i, S_i).
 *   vector<P<ll>> pts;
 *   rep(i, 0, n) pts.push_back({i, pre[i]});
 *   auto [a, b] = closest_pair(pts);
 *   printf("%lld\n", (a - b).dist2());
 * Notes: "Minimum distance between any two" of $10^5$ or more
 *   points, or "minimise $(i - j)^2 + g(i, j)^2$" (a distance in
 *   disguise), call for this. Manhattan distance becomes the
 *   $L_oo$ distance $max(|Delta x|, |Delta y|)$ of the points
 *   $(x + y, x - y)$: measure that instead of `dist2`, use `best`
 *   itself as the window and start it at `4e9 + 1`, above any
 *   distance (the type's maximum would overflow `p - d`, and a
 *   tie with the start could return its dummy points). The
 *   closest red-blue pair cannot use this sweep (points of one
 *   color may be close together, so the window is not small); it
 *   needs a Voronoi diagram of one color.
 * Source: kactl ClosestPair (CC0).
 * Status: stress-tested against brute force
 */
#pragma once
#include "point.hpp"

template<class T> array<P<T>, 2> closest_pair(vector<P<T>> v) {
  sort(all(v), [](P<T> a, P<T> b) { return a.y < b.y; });
  pair<T, array<P<T>, 2>> best{numeric_limits<T>::max(), {}};
  set<P<T>> s;
  int j = 0;
  for (P<T> p : v) {
    P<T> d{T(sqrt(best.first)) + is_integral_v<T>, 0};
    while (v[j].y < p.y - d.x) s.erase(v[j++]);
    auto lo = s.lower_bound(p - d), hi = s.upper_bound(p + d);
    for (; lo != hi; ++lo)
      best = min(best, {(*lo - p).dist2(), {*lo, p}});
    s.insert(p);
  }
  return best.second;
}
