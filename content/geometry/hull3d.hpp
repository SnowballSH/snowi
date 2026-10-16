/**
 * Title: 3D convex hull
 * Description: Triangles `{a, b, c}` (indices) covering the surface
 *   of the convex hull, counter-clockwise seen from outside: the
 *   normal `(p[b] - p[a]).cross(p[c] - p[a])` points out. A face
 *   with more than three points comes as several triangles.
 * Precondition: Not all points coplanar; duplicates are fine. Exact
 *   for `P3<ll>` with $|x|, |y|, |z| <= 5 dot 10^5$.
 * Time: $O(n^2)$ time and $n^2$ bytes.
 * Usage:
 *   vector<array<int, 3>> fs = hull3d(pts); // vector<P3<ll>>
 *   ll vol6 = 0; // 6 times the volume
 *   for (auto [a, b, c] : fs)
 *     vol6 += pts[a].dot(pts[b].cross(pts[c]));
 * Tweak: Delaunay triangulation of points in the plane, not all on
 *   one line or circle: `hull3d` of $(x, y, x^2 + y^2)$, keeping the
 *   faces whose normal has $z < 0$ (clockwise in the plane); exact
 *   for $|x|, |y| <= 2 dot 10^4$.
 *   Its edges contain the Euclidean MST and every nearest neighbour;
 *   Voronoi vertices are the circumcenters.
 * Notes: Incremental: a new point removes the faces it sees
 *   strictly and joins the horizon edges to itself. A new face keeps
 *   the direction of the horizon edge from the removed face, so it is
 *   oriented correctly even when it is coplanar with its neighbour;
 *   hence exact arithmetic handles coplanar points. Cocircular points
 *   in the Delaunay lift give an arbitrary triangulation of their
 *   polygon. Surface area: half the sum of the normals' lengths.
 * Source: Incremental construction as in kactl 3dHull (CC0),
 *   orientation by edge direction.
 * Status: stress-tested against brute force
 */
#pragma once

template<class T> struct P3 {
  T x, y, z;
  P3 operator-(P3 p) const {
    return {x - p.x, y - p.y, z - p.z};
  }
  P3 cross(P3 p) const {
    return {y * p.z - z * p.y, z * p.x - x * p.z,
            x * p.y - y * p.x};
  }
  T dot(P3 p) const { return x * p.x + y * p.y + z * p.z; }
};
template<class T>
vector<array<int, 3>> hull3d(const vector<P3<T>>& p) {
  int n = sz(p), a = 0, b = 1, c, d;
  auto nz = [](P3<T> v) { return v.x || v.y || v.z; };
  while (!nz(p[b] - p[a])) b++;
  for (c = b + 1; !nz((p[b] - p[a]).cross(p[c] - p[a])); c++);
  P3<T> q = (p[b] - p[a]).cross(p[c] - p[a]);
  for (d = c + 1; q.dot(p[d] - p[a]) == 0; d++);
  if (q.dot(p[d] - p[a]) > 0) swap(b, c);
  struct F {
    array<int, 3> v;
    P3<T> q;
  };
  vector<F> fs;
  // E[x][y]: the edge x -> y belongs to a live face
  vector<vector<char>> E(n, vector<char>(n));
  auto add = [&](int x, int y, int z) {
    fs.push_back({{x, y, z}, (p[y] - p[x]).cross(p[z] - p[x])});
    E[x][y] = E[y][z] = E[z][x] = 1;
  };
  add(a, b, c), add(a, d, b), add(b, d, c), add(c, d, a);
  rep(i, 0, n) {
    if (i == a || i == b || i == c || i == d) continue;
    vector<F> keep, vis; // faces that see p[i] strictly
    for (F& f : fs)
      (f.q.dot(p[i] - p[f.v[0]]) > 0 ? vis : keep).push_back(f);
    if (vis.empty()) continue;
    for (F& f : vis)
      rep(j, 0, 3) E[f.v[j]][f.v[(j + 1) % 3]] = 0;
    fs = keep;
    for (F& f : keep) // horizon edges: the face across is gone
      rep(j, 0, 3) {
        int x = f.v[j], y = f.v[(j + 1) % 3];
        if (!E[y][x]) add(y, x, i);
      }
  }
  vector<array<int, 3>> res;
  for (F& f : fs) res.push_back(f.v);
  return res;
}
