#include "../testlib.hpp"
#include "../../content/geometry/polygon.hpp"

// ASan's fake stacks make each by-value call of a small struct
// about 1000x slower on this toolchain; geometry code is all such
// calls and never returns pointers to locals.
extern "C" const char* __asan_default_options() {
  return "detect_stack_use_after_return=0";
}

typedef __int128 LL;
typedef P<ll> Pl;
typedef P<double> Pd;

Pl rnd_pt(ll c) { return {rnd(-c, c), rnd(-c, c)}; }
Pd to_d(Pl p) { return Pd(p.x, p.y); }
vector<Pd> to_d(const vector<Pl>& v) {
  vector<Pd> r;
  for (Pl p : v) r.push_back(to_d(p));
  return r;
}
bool near(double a, double b, double tol = 1e-7) {
  return fabs(a - b) <= tol * max(1.0, fabs(b));
}

// p on segment [a, b], exact, from the line equation.
bool exact_on(Pl a, Pl b, Pl p) {
  LL A = b.y - a.y, B = a.x - b.x;
  return A * p.x + B * p.y == A * a.x + B * a.y && min(a.x, b.x) <= p.x &&
         p.x <= max(a.x, b.x) && min(a.y, b.y) <= p.y && p.y <= max(a.y, b.y);
}
// 1 / 0 / -1 by winding number (angle sum), boundary exact.
int brute_in(const vector<Pl>& v, Pl p) {
  int n = sz(v);
  rep(i, 0, n)
    if (exact_on(v[i], v[(i + 1) % n], p)) return 0;
  double w = 0;
  rep(i, 0, n) {
    Pl a = v[i] - p, b = v[(i + 1) % n] - p;
    double d = atan2(b.y, b.x) - atan2(a.y, a.x);
    while (d > M_PI) d -= 2 * M_PI;
    while (d < -M_PI) d += 2 * M_PI;
    w += d;
  }
  return fabs(w) > M_PI ? 1 : -1;
}

// Is the closed polygon simple? Neighbours share only their vertex.
bool is_simple(const vector<Pl>& v) {
  int n = sz(v);
  rep(i, 0, n)
    rep(j, i + 1, n) {
      auto r = seg_inter(v[i], v[(i + 1) % n], v[j], v[(j + 1) % n]);
      bool adjacent = j == i + 1 || (i == 0 && j == n - 1);
      if (sz(r) > (adjacent ? 1 : 0)) return false;
    }
  return true;
}
// Make counter-clockwise: the turn at the smallest vertex is left.
vector<Pl> make_ccw(vector<Pl> v) {
  int n = sz(v), k = int(min_element(all(v)) - v.begin());
  Pl a = v[(k + n - 1) % n], b = v[k], c = v[(k + 1) % n];
  if ((b - a).cross(c - b) < 0) reverse(all(v));
  return v;
}
// Star-shaped around an interior point (may be non-convex).
vector<Pl> star_polygon(int n, ll C) {
  while (true) {
    Pl c = rnd_pt(C / 2);
    vector<Pl> v;
    rep(i, 0, n) {
      Pl p = rnd_pt(C);
      if (p != c) v.push_back(p);
    }
    sort(all(v), [&](Pl a, Pl b) { return angle_less(a - c, b - c); });
    bool ok = sz(v) >= 3;
    rep(i, 0, sz(v)) ok &= (v[i] - c).cross(v[(i + 1) % sz(v)] - c) > 0;
    if (ok) return v;
  }
}
// x-monotone from a random split into two chains; rejected unless
// simple. Covers zigzags that are not star-shaped.
vector<Pl> monotone_polygon(int n, ll C) {
  while (true) {
    vector<Pl> pts(n);
    for (auto& p : pts) p = rnd_pt(C);
    sort(all(pts)), pts.erase(unique(all(pts)), pts.end());
    if (sz(pts) < 3) continue;
    vector<Pl> lo{pts[0]}, hi;
    rep(i, 1, sz(pts) - 1) (rnd(0, 1) ? lo : hi).push_back(pts[i]);
    lo.push_back(pts.back());
    lo.insert(lo.end(), hi.rbegin(), hi.rend());
    if (is_simple(lo)) return make_ccw(lo);
  }
}

// Area and moments of polygon (intersected with the left side of
// s -> e if cut) by exact Gauss quadrature over vertical slabs.
array<double, 3> slab_brute(const vector<Pd>& v, bool cut, Pd s, Pd e) {
  int n = sz(v);
  vector<double> xs;
  for (Pd p : v) xs.push_back(p.x);
  if (cut) {
    xs.push_back(s.x);
    rep(i, 0, n) {
      Pd a = v[i], b = v[(i + 1) % n];
      double sa = (e - s).cross(a - s), sb = (e - s).cross(b - s);
      if ((sa < 0 && sb > 0) || (sa > 0 && sb < 0))
        xs.push_back(a.x + (b.x - a.x) * sa / (sa - sb));
    }
  }
  sort(all(xs));
  double area = 0, mx = 0, my = 0;
  rep(k, 0, sz(xs) - 1) {
    double w = xs[k + 1] - xs[k], m = (xs[k] + xs[k + 1]) / 2;
    if (w < 1e-12) continue;
    for (double g : {-1.0, 1.0}) {
      double x = m + g * w / (2 * sqrt(3.0)), len = 0, mom = 0;
      if (cut && s.x == e.x && (e.y > s.y ? x > s.x : x < s.x)) continue;
      vector<double> ys;
      rep(i, 0, n) {
        Pd a = v[i], b = v[(i + 1) % n];
        if ((a.x < x) != (b.x < x))
          ys.push_back(a.y + (b.y - a.y) * (x - a.x) / (b.x - a.x));
      }
      sort(all(ys));
      for (int i = 0; i + 1 < sz(ys); i += 2) {
        double lo = ys[i], hi = ys[i + 1];
        if (cut && s.x != e.x) {
          double yl = s.y + (e.y - s.y) * (x - s.x) / (e.x - s.x);
          if (e.x > s.x) lo = max(lo, yl);
          else hi = min(hi, yl);
        }
        if (hi > lo) len += hi - lo, mom += (hi * hi - lo * lo) / 2;
      }
      area += w / 2 * len, mx += w / 2 * x * len, my += w / 2 * mom;
    }
  }
  return {area, mx, my};
}

void check_polygon(vector<Pl> v) {
  int n = sz(v);
  ll C = 0;
  for (Pl p : v) C = max({C, llabs(p.x), llabs(p.y)});
  // Pick's theorem from brute classification of lattice points.
  ll inside = 0, boundary = 0;
  rep(x, -C - 1, C + 2)
    rep(y, -C - 1, C + 2) {
      Pl p{x, y};
      int got = in_polygon(v, p), want = brute_in(v, p);
      CHECK_EQ(got, want);
      inside += want == 1, boundary += want == 0;
    }
  CHECK_EQ(area2(v), 2 * inside + boundary - 2);
  vector<Pl> rv(v.rbegin(), v.rend());
  CHECK_EQ(area2(rv), -area2(v));
  vector<Pd> d = to_d(v);
  CHECK_EQ(area2(d), (double)area2(v));
  auto [area, mx, my] = slab_brute(d, false, {}, {});
  CHECK(near(area, area2(v) / 2.0, 1e-9));
  Pd g = centroid(d), rg = centroid(to_d(rv));
  CHECK(near(g.x, mx / area) && near(g.y, my / area));
  CHECK(near(rg.x, g.x) && near(rg.y, g.y));
  rep(k, 0, 4) {
    Pl s = rnd_pt(C + 1), e = rnd_pt(C + 1);
    if (k == 1) s = v[rnd(0, n - 1)];
    if (k == 2) s = v[rnd(0, n - 1)], e = v[rnd(0, n - 1)];
    if (k == 3) s = v[0], e = v[1]; // along an edge
    if (s == e) continue;
    for (bool rev : {false, true}) {
      vector<Pd> poly = rev ? to_d(rv) : d;
      vector<Pd> cut = polygon_cut(poly, to_d(s), to_d(e));
      auto [ca, cx, cy] = slab_brute(d, true, to_d(s), to_d(e));
      CHECK(near(area2(cut) / 2 * (rev ? -1 : 1), ca, 1e-9));
      for (Pd q : cut)
        CHECK(side_of(to_d(s), to_d(e), q) >= 0 ||
              fabs(line_dist(to_d(s), to_d(e), q)) < 1e-9);
    }
  }
}

// The same polygon moved by up to 1e9, as P<double> (what centroid
// and polygon_cut take): precision must follow its size, not its
// position. Cut points there are rounded to about 6e-8, which a
// sliver's centroid amplifies; so the cut's centroid is compared
// with that of the same rounded points moved back to the origin.
void check_far(const vector<Pl>& v) {
  Pl B = rnd_pt(1'000'000'000);
  vector<Pd> d = to_d(v), far;
  for (Pl p : v) far.push_back(to_d(p + B));
  CHECK_EQ(area2(far), (double)area2(v));
  Pd g = centroid(d), gf = centroid(far) - to_d(B);
  CHECK(near(gf.x, g.x, 1e-6) && near(gf.y, g.y, 1e-6));
  Pl s = rnd_pt(8), e = rnd_pt(8);
  if (s == e) return;
  double ca = slab_brute(d, true, to_d(s), to_d(e))[0];
  vector<Pd> cut = polygon_cut(far, to_d(s + B), to_d(e + B)), back;
  CHECK(near(area2(cut) / 2, ca, 1e-5));
  for (Pd p : cut) back.push_back(p - to_d(B));
  CHECK(near(area2(cut), area2(back), 1e-9));
  if (ca < 1e-3) return;
  Pd gc = centroid(cut) - to_d(B), gb = centroid(back);
  CHECK(near(gc.x, gb.x, 1e-6) && near(gc.y, gb.y, 1e-6));
}

int main() {
  rep(it, 0, iters(1500)) {
    int n = (int)rnd(3, 10);
    ll C = rnd(2, 7);
    vector<Pl> v = it % 2 ? star_polygon(n, C) : monotone_polygon(n, C);
    check_polygon(v);
    rep(k, 0, 5) check_far(v);
  }
  // Reproducer: area2 128 instead of 36, centroid off by 2e8.
  double B = 987654321;
  vector<Pd> t3{{B, B + 7}, {B + 3, B + 1}, {B + 5, B + 9}};
  CHECK_EQ(area2(t3), 36.0);
  Pd g3 = centroid(t3);
  CHECK(fabs(g3.x - (B + 8 / 3.0)) < 1e-6 &&
        fabs(g3.y - (B + 17 / 3.0)) < 1e-6);
  // The Example: cut a rectangle at 1e9 in half, then area2 (exact
  // here: every cut point is a half-integer).
  rep(it, 0, iters(5000)) {
    Pd o = to_d(rnd_pt(1'000'000'000)), m = o + Pd{0.5 * rnd(1, 199), 0};
    double w = rnd(1, 100), h = rnd(1, 100);
    vector<Pd> rect{o, o + Pd{w, 0}, o + Pd{w, h}, o + Pd{0, h}};
    vector<Pd> left = polygon_cut(rect, m, m + Pd{0, 1});
    double lw = min(m.x - o.x, w);
    CHECK_EQ(area2(left), 2 * lw * h);
    Pd g = centroid(left) - o;
    CHECK(near(g.x, lw / 2, 1e-6) && near(g.y, h / 2, 1e-6));
  }
  // Hand-made: comb (not star-shaped from any lattice point inside
  // all teeth), collinear vertices, horizontal edges at query y.
  vector<Pl> comb{{0, 0}, {6, 0}, {6, 4}, {5, 4}, {5, 1}, {4, 1},
                  {4, 4}, {3, 4}, {3, 1}, {2, 1}, {2, 4}, {0, 4}};
  vector<Pl> line_pts{{-3, -3}, {0, -3}, {3, -3}, {3, 0},
                      {3, 3},   {0, 3},  {-3, 3}, {-3, 0}};
  vector<Pl> spiral{{0, 0}, {7, 0}, {7, 7}, {0, 7}, {0, 2}, {5, 2},
                    {5, 5}, {2, 5}, {2, 4}, {4, 4}, {4, 3}, {1, 3},
                    {1, 6}, {6, 6}, {6, 1}, {0, 1}};
  for (auto& v : {comb, line_pts, spiral}) {
    CHECK(is_simple(v));
    check_polygon(v);
    rep(k, 0, 20) check_far(v);
  }
  CHECK_EQ(area2(comb), 2 * 18LL); // 6 x 4 minus two 1 x 3 notches
  Pd half_sq = centroid(polygon_cut(to_d(comb), Pd{1, 0}, Pd{1, 1}));
  CHECK(near(half_sq.x, 0.5) && near(half_sq.y, 2)); // x <= 1 part
  CHECK(polygon_cut(to_d(comb), Pd{0, 0}, Pd{0, 1}).size() == 2);

  // Empty and one-vertex polygons.
  CHECK_EQ(area2(vector<Pl>{}), 0LL);
  CHECK_EQ(in_polygon(vector<Pl>{}, Pl{0, 0}), -1);
  CHECK_EQ(in_polygon(vector<Pl>{{1, 1}}, Pl{1, 1}), 0);
  CHECK_EQ(in_polygon(vector<Pl>{{1, 1}}, Pl{0, 1}), -1);
  CHECK(polygon_cut(vector<Pd>{}, Pd{0, 0}, Pd{1, 0}).empty());

  // Extremes: coordinates up to 1e9 stay exact in ll.
  ll C = 1'000'000'000;
  vector<Pl> sq{{-C, -C}, {C, -C}, {C, C}, {-C, C}};
  CHECK_EQ(area2(sq), 8'000'000'000'000'000'000LL);
  rep(it, 0, iters(20000)) {
    vector<Pl> tri{rnd_pt(C), rnd_pt(C), rnd_pt(C)};
    LL want = 0;
    rep(i, 0, 3) {
      Pl a = tri[i], b = tri[(i + 1) % 3];
      want += (LL)a.x * b.y - (LL)a.y * b.x;
    }
    CHECK(area2(tri) == want);
    if (want < 0) swap(tri[0], tri[1]);
    Pl p = it % 3 ? rnd_pt(C) : tri[0] + (tri[1] - tri[0]) / 2;
    if (it % 7 == 0) p = (tri[0] + tri[1] + tri[2]) / 3;
    int expect = 1;
    rep(i, 0, 3) {
      Pl a = tri[i], b = tri[(i + 1) % 3];
      LL c = (LL)(b.x - a.x) * (p.y - a.y) - (LL)(b.y - a.y) * (p.x - a.x);
      if (c < 0) expect = -1;
      else if (c == 0 && expect == 1) expect = exact_on(a, b, p) ? 0 : -1;
    }
    if (want != 0) CHECK_EQ(in_polygon(tri, p), expect);
  }

  // Large: 1e5-gon (star-shaped), area vs __int128, a few queries.
  int n = 100'000;
  vector<Pl> big;
  rep(i, 0, n) {
    double t = 2 * M_PI * i / n, r = rnd(5e8, 1e9);
    big.push_back({llround(r * cos(t)), llround(r * sin(t))});
  }
  double start = elapsed();
  LL want = 0;
  rep(i, 0, n) {
    Pl a = big[i], b = big[(i + 1) % n];
    want += (LL)a.x * b.y - (LL)a.y * b.x;
  }
  CHECK(area2(big) == want);
  vector<Pd> bd = to_d(big);
  Pd g = centroid(bd);
  CHECK(fabs(g.x) < 5e7 && fabs(g.y) < 5e7);
  rep(q, 0, 20) {
    Pl p = q % 2 ? rnd_pt(C) : big[rnd(0, n - 1)];
    CHECK_EQ(in_polygon(big, p), brute_in(big, p));
  }
  auto cut = polygon_cut(bd, Pd{0, -1}, Pd{0, 1});
  CHECK(cut.size() > (size_t)n / 3);
  fprintf(stderr, "  1e5-gon area, centroid, 20 queries, cut: %.2fs\n",
          elapsed() - start);
  puts("OK geometry/polygon");
}
