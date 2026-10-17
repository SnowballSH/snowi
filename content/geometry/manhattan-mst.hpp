/**
 * Title: Manhattan MST
 * Description: Returns at most $4n$ edges `{w, i, j}`, with `w` the
 *   distance $|x_i - x_j| + |y_i - y_j|$, that contain a minimum
 *   spanning tree of the complete graph on the points: run Kruskal
 *   on them.
 * Precondition: `P<ll>`; equal points are allowed (distance 0).
 * Time: $O(n log n)$.
 * Usage:
 *   auto e = manhattan_mst(pts);
 *   sort(all(e)); // then Kruskal with a DSU
 * Tweak: Chebyshev distance $max(|Delta x|, |Delta y|)$: use the
 *   points $(x + y, x - y)$ and halve the weights.
 * Notes: Each point is joined to its nearest point in four of the
 *   eight octants, one of each opposite pair: if an MST edge leaves
 *   $p$ into octant $O$, it enters $q$ from the opposite octant, and
 *   among the points of one octant of $q$ the nearest one is at
 *   least as close to every other point of that octant, so only
 *   nearest neighbours are needed. Each sweep keeps the points that
 *   have not found their partner yet, keyed by $-y$.
 * Source: kactl ManhattanMST (CC0).
 * Status: stress-tested against brute force
 */
#pragma once
#include "point.hpp"

vector<array<ll, 3>> manhattan_mst(vector<P<ll>> ps) {
  vi id(sz(ps));
  iota(all(id), 0);
  vector<array<ll, 3>> edges;
  rep(k, 0, 4) {
    sort(all(id), [&](int i, int j) {
      return ps[i].x + ps[i].y < ps[j].x + ps[j].y;
    });
    map<ll, int> sweep; // -y -> point without a partner yet
    for (int i : id) {
      for (auto it = sweep.lower_bound(-ps[i].y);
           it != sweep.end(); sweep.erase(it++)) {
        int j = it->second;
        P<ll> d = ps[i] - ps[j];
        if (d.y > d.x) break;
        edges.push_back({d.y + d.x, i, j});
      }
      sweep[-ps[i].y] = i;
    }
    for (P<ll>& p : ps)
      if (k & 1) p.x = -p.x;
      else swap(p.x, p.y);
  }
  return edges;
}
