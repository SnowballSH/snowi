/**
 * Title: Minkowski sum
 * Description: Returns the Minkowski sum $A + B = {a + b}$ of two
 *   convex polygons, counter-clockwise from its smallest point (by
 *   `x`, then `y`).
 * Precondition: `a`, `b` non-empty, convex, counter-clockwise and
 *   without collinear vertices, as from `convex_hull` (one or two
 *   points are fine). The sum is then the same kind of polygon,
 *   with at most $n + m$ vertices. Exact for `P<ll>` with
 *   $|x|, |y| <= 10^9$; the sum's coordinates reach $2 dot 10^9$.
 * Time: $O(n + m)$.
 * Usage:
 *   vector<P<ll>> s = minkowski(a, b);
 * Tweak: $A - B$: negate the points of `b` first (`p * -1`, a half
 *   turn keeps the order counter-clockwise). Sum of $k$ polygons:
 *   sum them in pairs, or merge all edges by angle at once.
 * Example:
 *   // Do convex polygons A and B + t meet? Iff t is in A - B.
 *   for (auto& p : b) p = p * -1;
 *   vector<P<ll>> d = minkowski(a, b);
 *   bool meet = in_hull(d, t) >= 0; // |coordinates| <= 5e8
 *   // Distance between disjoint A and B: from 0 to A - B.
 *   vector<P<double>> f;
 *   for (auto p : d) f.push_back(P<double>(p.x, p.y));
 *   double dist = 1e18;
 *   rep(i, 0, sz(f)) dist = min(dist,
 *       seg_dist(f[i], f[(i + 1) % sz(f)], P<double>{}));
 * Notes: The extreme point of $A + B$ in a direction is the sum
 *   of the extreme points of $A$ and $B$, so "maximise a linear
 *   function of $a + b$" is one query on the sum. The translations
 *   $t$ for which $A$ and $B + t$ overlap form $A - B$ (collision
 *   checks, free space of a moving robot); the distance between
 *   convex polygons is that from the origin to $A - B$. A sum has
 *   one vertex per distinct edge direction of the summands. Window
 *   problems ("best sum + length times $x$ over all subarrays")
 *   become upper hulls of sums $(i, S_i) + (-j, -S_j)$, built by
 *   divide and conquer on the index with merges like this one. For
 *   `in_hull` on a sum use `P<__int128>`, or input coordinates up
 *   to $5 dot 10^8$. Reverse clockwise input first; the merge
 *   itself starts both polygons at their smallest point.
 * Source: UCup S4-S5 Training Packet 20; cp-algorithms Minkowski
 *   sum (idea).
 * Status: stress-tested against brute force
 */
#pragma once
#include "point.hpp"

template<class T>
vector<P<T>> minkowski(vector<P<T>> a, vector<P<T>> b) {
  rotate(a.begin(), min_element(all(a)), a.end());
  rotate(b.begin(), min_element(all(b)), b.end());
  int n = sz(a), m = sz(b), i = 0, j = 0;
  a.push_back(a[0]), b.push_back(b[0]);
  vector<P<T>> res;
  while (i < n || j < m) {
    res.push_back(a[i] + b[j]);
    T c = 0; // > 0: a's edge turns less, take it first
    if (i < n && j < m)
      c = (a[i + 1] - a[i]).cross(b[j + 1] - b[j]);
    i += i < n && c >= 0, j += j < m && c <= 0;
  }
  return res;
}
