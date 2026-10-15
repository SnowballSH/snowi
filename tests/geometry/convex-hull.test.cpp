#include "../testlib.hpp"
#include "../../content/geometry/convex-hull.hpp"

typedef __int128 LL;
typedef P<ll> Pl;

Pl rnd_pt(ll c) { return {rnd(-c, c), rnd(-c, c)}; }
LL cr(Pl o, Pl a, Pl b) {
  return (LL)(a.x - o.x) * (b.y - o.y) - (LL)(a.y - o.y) * (b.x - o.x);
}
bool brute_on(Pl a, Pl b, Pl p) {
  return cr(a, b, p) == 0 && min(a.x, b.x) <= p.x && p.x <= max(a.x, b.x) &&
         min(a.y, b.y) <= p.y && p.y <= max(a.y, b.y);
}
bool in_triangle(Pl a, Pl b, Pl c, Pl p) { // closed, non-degenerate
  if (cr(a, b, c) == 0) return false;      // segments cover that case
  LL x = cr(a, b, p), y = cr(b, c, p), z = cr(c, a, p);
  return (x >= 0 && y >= 0 && z >= 0) || (x <= 0 && y <= 0 && z <= 0);
}
// Hull vertices: points not in a closed triangle or segment of the
// others.
set<Pl> brute_vertices(vector<Pl> pts) {
  sort(all(pts)), pts.erase(unique(all(pts)), pts.end());
  int n = sz(pts);
  set<Pl> res;
  rep(i, 0, n) {
    bool covered = false;
    rep(a, 0, n)
      rep(b, a + 1, n) {
        if (a == i || b == i) continue;
        covered |= brute_on(pts[a], pts[b], pts[i]);
        rep(c, b + 1, n)
          if (c != i) covered |= in_triangle(pts[a], pts[b], pts[c], pts[i]);
      }
    if (!covered) res.insert(pts[i]);
  }
  return res;
}
int brute_in_hull(const vector<Pl>& h, Pl p) {
  int n = sz(h);
  if (n <= 2) return n && brute_on(h[0], h[n - 1], p) ? 0 : -1;
  int res = 1;
  rep(i, 0, n) {
    LL c = cr(h[i], h[(i + 1) % n], p);
    if (c < 0) return -1;
    if (c == 0) res = 0;
  }
  return res;
}
void check_hull(const vector<Pl>& pts) {
  vector<Pl> h = convex_hull(pts);
  set<Pl> want = brute_vertices(pts);
  CHECK_EQ(sz(h), sz(want));
  CHECK((set<Pl>(all(h)) == want));
  CHECK(h[0] == *want.begin());
  int n = sz(h);
  if (n >= 3) // strictly convex, counter-clockwise, winding once
    rep(i, 0, n) {
      CHECK(cr(h[i], h[(i + 1) % n], h[(i + 2) % n]) > 0);
      if (i >= 1 && i + 1 < n) CHECK(cr(h[0], h[i], h[i + 1]) > 0);
    }
  ll best = 0;
  for (Pl a : pts)
    for (Pl b : pts) best = max(best, (a - b).dist2());
  auto [a, b] = hull_diameter(h);
  CHECK_EQ((a - b).dist2(), best);
  CHECK(want.count(a) && want.count(b));
}

int main() {
  rep(it, 0, iters(1500)) {
    int n = (int)rnd(1, 10);
    ll C = it % 3 ? rnd(1, 4) : 1'000'000'000;
    vector<Pl> pts(n);
    for (auto& p : pts) p = rnd_pt(C);
    if (it % 4 == 0) // collinear runs
      for (auto& p : pts) p = pts[0] + (Pl{1, 2} * rnd(-3, 3));
    if (it % 5 == 0) pts.push_back(pts[0]);
    check_hull(pts);
    vector<Pl> h = convex_hull(pts);
    if (C < 10)
      rep(x, -C - 1, C + 2)
        rep(y, -C - 1, C + 2)
          CHECK_EQ(in_hull(h, Pl{x, y}), brute_in_hull(h, Pl{x, y}));
    rep(q, 0, 20) {
      Pl p = rnd_pt(C);
      if (q % 3 == 0) p = h[rnd(0, sz(h) - 1)];
      if (q % 3 == 1) {
        Pl a = h[rnd(0, sz(h) - 1)], b = h[rnd(0, sz(h) - 1)];
        p = a + (b - a) * rnd(0, 1) / 2; // edge midpoints, chords
      }
      CHECK_EQ(in_hull(h, p), brute_in_hull(h, p));
    }
  }
  // Degenerate and hand-made cases.
  CHECK(convex_hull(vector<Pl>{}).empty());
  CHECK((convex_hull(vector<Pl>{{2, 2}, {2, 2}}) == vector<Pl>{{2, 2}}));
  vector<Pl> sq{{0, 0}, {2, 0}, {2, 2}, {0, 2}, {1, 1}, {1, 0}, {2, 1}};
  CHECK((convex_hull(sq) == vector<Pl>{{0, 0}, {2, 0}, {2, 2}, {0, 2}}));
  vector<Pl> h = convex_hull(sq);
  CHECK_EQ(in_hull(h, Pl{1, 1}), 1);
  CHECK_EQ(in_hull(h, Pl{0, 0}), 0);
  CHECK_EQ(in_hull(h, Pl{0, 1}), 0);
  CHECK_EQ(in_hull(h, Pl{3, 1}), -1);
  CHECK_EQ(in_hull(h, Pl{0, 3}), -1);
  CHECK_EQ(in_hull(vector<Pl>{}, Pl{0, 0}), -1);
  // Sorted by x only (the sort by x, then y, is still needed).
  vector<Pl> xs{{0, 2}, {0, 3}, {0, 0}, {1, 0}, {2, 3}, {3, 0}};
  CHECK((convex_hull(xs) == vector<Pl>{{0, 0}, {3, 0}, {2, 3}, {0, 3}}));
  vector<Pl> seg = convex_hull(vector<Pl>{{0, 0}, {3, 3}, {1, 1}});
  CHECK((seg == vector<Pl>{{0, 0}, {3, 3}}));
  CHECK_EQ(in_hull(seg, Pl{2, 2}), 0);
  CHECK_EQ(in_hull(seg, Pl{4, 4}), -1);

  // Large: 1e5 points; extremes of the coordinate range.
  ll C = 1'000'000'000;
  int n = 100'000;
  vector<Pl> pts(n);
  for (auto& p : pts) p = rnd_pt(C);
  rep(i, 0, 4) pts[i] = {i % 2 ? C : -C, i / 2 ? C : -C};
  double start = elapsed();
  h = convex_hull(pts);
  for (Pl p : pts) CHECK(in_hull(h, p) >= 0);
  auto [a, b] = hull_diameter(h);
  fprintf(stderr, "  hull + 1e5 in_hull + diameter: %.2fs\n",
          elapsed() - start);
  CHECK((h[0] == Pl{-C, -C}) && sz(h) == 4);
  CHECK_EQ((a - b).dist2(), 8 * C * C);
  // Convex position: points on a circle, almost all on the hull.
  rep(i, 0, n) {
    double t = 2 * M_PI * i / n;
    pts[i] = {llround(C * cos(t)), llround(C * sin(t))};
  }
  h = convex_hull(pts);
  CHECK(sz(h) > n / 2);
  rep(i, 0, sz(h)) CHECK(cr(h[i], h[(i + 1) % sz(h)], h[(i + 2) % sz(h)]) > 0);
  for (Pl p : pts) CHECK(in_hull(h, p) == 0);
  auto far = hull_diameter(h);
  CHECK((far[0] - far[1]).dist2() >= 4 * C * C - 4 * C);
  puts("OK geometry/convex-hull");
}
