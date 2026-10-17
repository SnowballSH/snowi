#include "../testlib.hpp"
#include "../../content/geometry/circles.hpp"

typedef P<double> Pd;

Pd rnd_pt(double c) { return {rnd_real(-c, c), rnd_real(-c, c)}; }
Pd polar(double r, double t) { return {r * cos(t), r * sin(t)}; }
bool near(double a, double b, double tol = 1e-9) {
  return fabs(a - b) <= tol * max(1.0, fabs(b));
}
bool near(Pd a, Pd b, double tol = 1e-9) {
  return (a - b).dist() <= tol * max(1.0, b.dist());
}
// got and want hold the same points (up to order).
bool same_set(vector<Pd> got, vector<Pd> want, double tol = 1e-7) {
  if (sz(got) != sz(want)) return false;
  for (Pd w : want) {
    bool found = false;
    for (Pd g : got) found |= near(g, w, tol);
    if (!found) return false;
  }
  return true;
}
// Minimum of a convex function on [lo, hi].
template<class F> double ternary(double lo, double hi, F f) {
  rep(it, 0, 200) {
    double m1 = lo + (hi - lo) / 3, m2 = hi - (hi - lo) / 3;
    if (f(m1) < f(m2)) hi = m2;
    else lo = m1;
  }
  return f((lo + hi) / 2);
}
// Points at angle t0 +- acos(x) on circle (c, r), or none / one.
vector<Pd> by_angle(Pd c, double r, double t0, double x) {
  if (x > 1 || x < -1) return {};
  double al = acos(x);
  return {c + polar(r, t0 + al), c + polar(r, t0 - al)};
}

void check_line(Pd c, double r, Pd a, Pd b) {
  vector<Pd> v = circle_line(c, r, a, b);
  double dmin = ternary(-1e4, 1e4,
                        [&](double t) { return (a + (b - a) * t - c).dist(); });
  if (dmin > r + 1e-6) CHECK(v.empty());
  if (dmin < r - 1e-6) CHECK_EQ(sz(v), 2);
  for (Pd q : v) {
    CHECK(near((q - c).dist(), r, 1e-6));
    CHECK(fabs((b - a).cross(q - a)) / (b - a).dist() < 1e-6);
  }
  if (sz(v) == 2) CHECK((v[0] - a).dot(b - a) < (v[1] - a).dot(b - a));
}

void check_inter(Pd a, double r1, Pd b, double r2) {
  vector<Pd> v = circle_inter(a, r1, b, r2);
  double d = (b - a).dist();
  if (d > r1 + r2 + 1e-6 || d < fabs(r1 - r2) - 1e-6) CHECK(v.empty());
  if (d < r1 + r2 - 1e-6 && d > fabs(r1 - r2) + 1e-6) {
    // Law of cosines on circle a: independent of the formula used.
    double x = (r1 * r1 + d * d - r2 * r2) / (2 * r1 * d);
    CHECK(same_set(v, by_angle(a, r1, (b - a).angle(), x)));
    CHECK((b - a).cross(v[0] - a) > 0);
  }
  for (Pd q : v) {
    CHECK(near((q - a).dist(), r1, 1e-6));
    CHECK(near((q - b).dist(), r2, 1e-6));
  }
}

// Tangents: outer (r2 > 0), inner (r2 < 0) or from a point (0).
void check_tangents(Pd c1, double r1, Pd c2, double r2) {
  auto t = tangents(c1, r1, c2, r2);
  double d = (c2 - c1).dist(), dr = r1 - r2;
  if (d < fabs(dr) - 1e-6) CHECK(t.empty());
  if (d > fabs(dr) + 1e-6) {
    // The unit normal n of the line has n . (c2 - c1) = r1 - r2.
    vector<Pd> want = by_angle(c1, r1, (c2 - c1).angle(), dr / d);
    vector<Pd> got;
    for (auto [p, q] : t) got.push_back(p);
    CHECK(same_set(got, want));
  }
  for (auto [p, q] : t) {
    CHECK(near((p - c1).dist(), r1, 1e-6));
    CHECK(near((q - c2).dist(), fabs(r2), 1e-6));
    if ((q - p).dist() < 1e-6) continue; // touching circles
    Pd dir = (q - p).unit();
    CHECK(fabs(dir.dot(p - c1)) < 1e-6 && fabs(dir.dot(q - c2)) < 1e-6);
    double s1 = dir.cross(c1 - p), s2 = dir.cross(c2 - p);
    if (r2 != 0) CHECK(r2 > 0 ? s1 * s2 > 0 : s1 * s2 < 0);
  }
}

// Brute force minimum enclosing circle over pairs and triples.
pair<Pd, double> brute_mec(const vector<Pd>& ps) {
  pair<Pd, double> best{ps[0], 0};
  bool found = false;
  auto consider = [&](Pd o, double r) {
    for (Pd p : ps)
      if ((p - o).dist() > r + 1e-9) return;
    if (!found || r < best.second) best = {o, r}, found = true;
  };
  int n = sz(ps);
  consider(ps[0], 0);
  rep(i, 0, n)
    rep(j, i + 1, n) {
      Pd o = (ps[i] + ps[j]) / 2;
      consider(o, (ps[i] - o).dist());
      rep(k, j + 1, n) {
        // Perpendicular bisectors: 2 (b - a) . o = |b|^2 - |a|^2.
        Pd a = ps[i], b = ps[j], c = ps[k];
        double a1 = 2 * (b.x - a.x), b1 = 2 * (b.y - a.y);
        double c1 = b.dist2() - a.dist2();
        double a2 = 2 * (c.x - a.x), b2 = 2 * (c.y - a.y);
        double c2 = c.dist2() - a.dist2(), det = a1 * b2 - a2 * b1;
        if (fabs(det) < 1e-12) continue;
        Pd q{(c1 * b2 - c2 * b1) / det, (a1 * c2 - a2 * c1) / det};
        consider(q, (q - a).dist());
      }
    }
  return best;
}

// Area of disc (c, r) inside a simple polygon by Green's theorem, a
// different method from circle_poly: the polygon's edges cut at the
// circle (pieces inside the disc), then the circle's arcs between
// crossings (arcs inside the polygon). Generic input only.
bool in_poly(const vector<Pd>& v, Pd p) {
  bool in = false;
  rep(i, 0, sz(v)) {
    Pd a = v[i], b = v[(i + 1) % sz(v)];
    if ((a.y > p.y) != (b.y > p.y) &&
        p.x < a.x + (b.x - a.x) * (p.y - a.y) / (b.y - a.y))
      in = !in;
  }
  return in;
}
double green_area(Pd c, double r, vector<Pd> v) {
  double sign = 1, area2 = 0;
  rep(i, 0, sz(v)) area2 += v[i].cross(v[(i + 1) % sz(v)]);
  if (area2 < 0) reverse(all(v)), sign = -1;
  for (Pd& p : v) p = p - c;
  double sum = 0;
  vector<double> angles;
  rep(i, 0, sz(v)) {
    Pd a = v[i], d = v[(i + 1) % sz(v)] - a;
    vector<double> ts{0, 1};
    // |a + d t| = r: d.d t^2 + 2 a.d t + a.a - r^2 = 0.
    double A = d.dist2(), B = 2 * a.dot(d), C = a.dist2() - r * r;
    double disc = B * B - 4 * A * C;
    if (disc > 0)
      for (double sg : {-1.0, 1.0}) {
        double t = (-B + sg * sqrt(disc)) / (2 * A);
        if (t > 0 && t < 1) {
          ts.push_back(t);
          Pd q = a + d * t;
          angles.push_back(atan2(q.y, q.x));
        }
      }
    sort(all(ts));
    rep(j, 0, sz(ts) - 1) {
      Pd p = a + d * ts[j], q = a + d * ts[j + 1];
      if (((p + q) / 2).dist() < r) sum += p.cross(q) / 2;
    }
  }
  if (angles.empty())
    return sign * (sum + (in_poly(v, Pd{r, 0}) ? M_PI * r * r : 0));
  sort(all(angles));
  rep(i, 0, sz(angles)) {
    double lo = angles[i];
    double hi = i + 1 < sz(angles) ? angles[i + 1] : angles[0] + 2 * M_PI;
    if (in_poly(v, polar(r, (lo + hi) / 2))) sum += r * r * (hi - lo) / 2;
  }
  return sign * sum;
}
// Clip v to the left of every edge of the convex polygon w (both
// orientations of v keep their signed area).
vector<Pd> clip(vector<Pd> v, const vector<Pd>& w) {
  rep(k, 0, sz(w)) {
    Pd s = w[k], e = w[(k + 1) % sz(w)];
    vector<Pd> out;
    rep(i, 0, sz(v)) {
      Pd cur = v[i], prev = v[(i + sz(v) - 1) % sz(v)];
      double a = (e - s).cross(cur - s), b = (e - s).cross(prev - s);
      if ((a < 0) != (b < 0)) out.push_back(cur + (prev - cur) * (a / (a - b)));
      if (a >= 0) out.push_back(cur);
    }
    v = out;
  }
  return v;
}
double signed_area(const vector<Pd>& v) {
  double a = 0;
  rep(i, 0, sz(v)) a += v[i].cross(v[(i + 1) % sz(v)]);
  return a / 2;
}
// Bounds from regular m-gons inside and around the circle.
pair<double, double> ngon_bounds(Pd c, double r, const vector<Pd>& v) {
  int m = 256;
  vector<Pd> in, out;
  rep(k, 0, m) {
    in.push_back(c + polar(r, 2 * M_PI * k / m));
    out.push_back(c + polar(r / cos(M_PI / m), 2 * M_PI * k / m));
  }
  double x = signed_area(clip(v, in)), y = signed_area(clip(v, out));
  return {min(x, y), max(x, y)};
}
// Star-shaped (so simple) polygon around o, counter-clockwise: all
// angular gaps below a half turn.
vector<Pd> rnd_star(Pd o, int n, double rmax) {
  vector<double> t(n);
  for (bool ok = false; !ok;) {
    for (double& x : t) x = rnd_real(0, 2 * M_PI);
    sort(all(t));
    ok = t[0] + 2 * M_PI - t[n - 1] < 3;
    rep(i, 1, n) ok &= t[i] - t[i - 1] < 3;
  }
  vector<Pd> v;
  for (double x : t) v.push_back(o + polar(rnd_real(0.05, rmax), x));
  return v;
}

int main() {
  rep(it, 0, iters(3000)) {
    Pd c = rnd_pt(10), a = rnd_pt(10), b = rnd_pt(10);
    double r = rnd_real(0.1, 10), t = rnd_real(0, 2 * M_PI);
    check_line(c, r, a, b);
    // Tangent line through T, perpendicular to the radius.
    Pd T = c + polar(r, t), dir = polar(1, t).perp();
    Pd ta = T - dir * rnd_real(0.5, 5), tb = T + dir * rnd_real(0.5, 5);
    vector<Pd> v = circle_line(c, r, ta, tb);
    CHECK(sz(v) == 1 && near(v[0], T, 1e-7));

    double r1 = rnd_real(0.1, 10), r2 = rnd_real(0.1, 10);
    check_inter(a, r1, b, r2);
    check_inter(a, r1, a + polar(r1 + r2 + rnd_real(-1, 1), t), r2);
    Pd ext = a + polar(r1 + r2, t), in = a + polar(fabs(r1 - r2), t);
    v = circle_inter(a, r1, ext, r2);
    CHECK(sz(v) == 1 && near(v[0], a + polar(r1, t), 1e-7));
    v = circle_inter(a, r1, in, r2);
    CHECK(sz(v) == 1);
    CHECK(near(v[0], a + polar(r1 < r2 ? -r1 : r1, t), 1e-7));
    CHECK(circle_inter(a, r1, a, r2).empty());

    check_tangents(a, r1, b, r2);
    check_tangents(a, r1, b, -r2);
    check_tangents(a, r1, b, 0.0);
    check_tangents(a, r1, a + polar(fabs(r1 - r2) + rnd_real(-1, 1), t), r2);
    CHECK_EQ(sz(tangents(a, r1, ext, r2)), 2);
    CHECK_EQ(sz(tangents(a, r1, ext, -r2)), 1);
    CHECK_EQ(sz(tangents(a, r1, in, r2)), 1);
    CHECK_EQ(sz(tangents(a, r1, a + polar(r1, t), 0.0)), 1);
    CHECK(tangents(a, r1, a, r2).empty());

    Pd p = rnd_pt(10), q = rnd_pt(10), s = rnd_pt(10);
    if (fabs((q - p).cross(s - p)) > 1e-2) {
      auto [o, rad] = circumcircle(p, q, s);
      CHECK(near((p - o).dist(), rad, 1e-7));
      CHECK(near((q - o).dist(), rad, 1e-7) && near((s - o).dist(), rad, 1e-7));
    }
  }
  // Exact integer configurations.
  Pd o{0, 0};
  CHECK((circle_line(o, 5.0, Pd{3, -10}, Pd{3, 10}) ==
         vector<Pd>{{3, -4}, {3, 4}}));
  CHECK((circle_line(o, 5.0, Pd{5, -1}, Pd{5, 1}) == vector<Pd>{{5, 0}}));
  CHECK(circle_line(o, 5.0, Pd{6, -1}, Pd{6, 1}).empty());
  CHECK((circle_inter(o, 5.0, Pd{8, 6}, 5.0) == vector<Pd>{{4, 3}}));
  CHECK(same_set(circle_inter(o, 5.0, Pd{6, 0}, 5.0), {{3, -4}, {3, 4}}));
  CHECK(same_set({tangents(o, 5.0, Pd{0, 10}, 5.0)[0].first,
                  tangents(o, 5.0, Pd{0, 10}, 5.0)[1].first},
                 {{5, 0}, {-5, 0}}));
  auto [cc, cr] = circumcircle(Pd{0, 0}, Pd{4, 0}, Pd{0, 3});
  CHECK(near(cc, Pd{2, 1.5}) && near(cr, 2.5));

  rep(it, 0, iters(1500)) {
    int n = (int)rnd(1, 8);
    vector<Pd> ps(n);
    for (auto& p : ps) p = it % 2 ? rnd_pt(10) : Pd(rnd(-3, 3), rnd(-3, 3));
    if (it % 5 == 0) ps.push_back(ps[0]);
    auto [ctr, rad] = mec(ps);
    auto [bc, br] = brute_mec(ps);
    CHECK(near(rad, br, 1e-7));
    for (Pd p : ps) CHECK((p - ctr).dist() <= rad * (1 + 1e-7) + 1e-12);
    CHECK(near(ctr, bc, 1e-6));
  }

  // Far from the origin the precision must follow the spread of the
  // points, not their position. Reproducer: NaN at 1e9 (a circle
  // through a collinear triple after a rounded boundary test).
  vector<Pd> repro{{0, 2}, {1, 0}, {2, 1}, {0, 2}, {0, 2}};
  rep(it, 0, iters(4000)) {
    bool fixed = it < iters(1000);
    Pd B = fixed ? Pd{1e9, 1e9} : Pd(rnd(-1e9, 1e9), rnd(-1e9, 1e9));
    double spread = vector<double>{1, 3, 10, 100}[it % 4];
    vector<Pd> ps = repro;
    if (!fixed) {
      ps.resize(rnd(1, 8));
      for (Pd& p : ps)
        p = it % 3 ? Pd(rnd(0, (ll)spread), rnd(0, (ll)spread))
                   : Pd{rnd_real(0, spread), rnd_real(0, spread)};
      if (it % 2) ps.push_back(ps[rnd(0, sz(ps) - 1)]);
    }
    for (Pd& p : ps) p = p + B;
    shuffle(all(ps), test_rng); // mec's own seed is the clock
    vector<Pd> rel;             // exact differences (Sterbenz)
    for (Pd p : ps) rel.push_back(p - B);
    auto [ctr, rad] = mec(ps);
    auto [bc, br] = brute_mec(rel);
    CHECK(near(rad, br, 1e-7));
    CHECK(near(ctr - B, bc, 1e-6));
    // The returned center is rounded to about 1e-7 at 1e9.
    for (Pd p : ps) CHECK((p - ctr).dist() <= rad * (1 + 1e-7) + 1e-6);
  }

  // circle_poly: generic polygons against Green's theorem.
  rep(it, 0, iters(3000)) {
    int n = (int)rnd(3, 9);
    Pd ctr0 = rnd_pt(5), c = rnd_pt(5);
    vector<Pd> v = rnd_star(ctr0, n, 6);
    if (it % 2) reverse(all(v));
    double r = vector<double>{0.3, 1, 3, 8, 20}[it % 5] * rnd_real(0.5, 1);
    double got = circle_poly(c, r, v), want = green_area(c, r, v);
    CHECK(near(got, want, 1e-7));
    reverse(all(v));
    CHECK(near(circle_poly(c, r, v), -want, 1e-7));
  }
  // Degenerate integer cases (vertices on the circle, tangent edges,
  // the center on a vertex or an edge): inside regular-polygon bounds
  // and close to the area for a slightly larger radius.
  int degenerate = 0;
  rep(it, 0, iters(600)) {
    vector<Pd> v;
    if (it % 4 == 0) v = {{0, 0}, {6, 0}, {6, 6}, {0, 6}};
    else
      rep(i, 0, rnd(3, 6)) v.push_back(Pd(rnd(-5, 5), rnd(-5, 5)));
    if (it % 8 == 4) reverse(all(v));
    int n = sz(v);
    // Keep simple polygons: no two non-adjacent edges touch.
    bool simple = fabs(signed_area(v)) > 1e-9;
    rep(i, 0, n) simple &= v[i] != v[(i + 1) % n];
    rep(i, 0, n)
      rep(j, i + 2, n) {
        if (i == 0 && j == n - 1) continue;
        Pd a = v[i], b = v[(i + 1) % n], p = v[j], q = v[(j + 1) % n];
        double d1 = (b - a).cross(p - a), d2 = (b - a).cross(q - a);
        double d3 = (q - p).cross(a - p), d4 = (q - p).cross(b - p);
        if (d1 * d2 <= 0 && d3 * d4 <= 0) simple = false;
      }
    if (!simple) continue;
    vector<Pd> cs{v[0], (v[0] + v[1]) / 2, {3, 3}};
    Pd c = it % 2 ? Pd(rnd(-3, 3), rnd(-3, 3)) : cs[it / 2 % 3];
    double r = (double)rnd(1, 5);
    auto [lo, hi] = ngon_bounds(c, r, v);
    double got = circle_poly(c, r, v);
    CHECK(lo - 1e-9 <= got && got <= hi + 1e-9);
    // A slightly larger radius makes the case generic.
    CHECK(near(got, green_area(c, r + 1e-8, v), 1e-6));
    degenerate++;
  }
  CHECK(degenerate > iters(600) / 4);
  {
    vector<Pd> sq{{0, 0}, {2, 0}, {2, 2}, {0, 2}};
    CHECK(near(circle_poly(Pd{0, 0}, 1.0, sq), M_PI / 4, 1e-12));
    CHECK(near(circle_poly(Pd{1, 1}, 1.0, sq), M_PI, 1e-12));
    CHECK(near(circle_poly(Pd{1, 1}, 5.0, sq), 4, 1e-12));
    CHECK(near(circle_poly(Pd{9, 9}, 1.0, sq), 0, 1e-12));
    CHECK(near(circle_poly(Pd{1, 0}, 1.0, sq), M_PI / 2, 1e-12));
  }

  int n = 100'000;
  vector<Pd> ps(n);
  rep(i, 0, n)
    ps[i] = polar(i % 10 ? rnd_real(0, 1e3) : 1e3, rnd_real(0, 2 * M_PI));
  double start = elapsed(), sum = 0;
  auto [ctr, rad] = mec(ps);
  rep(i, 0, n) {
    sum += sz(circle_line(ps[i], 1.0, ps[(i + 1) % n], o));
    sum += sz(circle_inter(ps[i], 3.0, ps[(i + 1) % n], 4.0));
    sum += sz(tangents(ps[i], 3.0, ps[(i + 1) % n], -4.0));
  }
  sum += circle_poly(Pd{1, 2}, 500.0, ps);
  fprintf(stderr, "  mec of 1e5 + 3e5 circle ops + circle_poly: %.2fs (%g)\n",
          elapsed() - start, sum);
  CHECK(rad <= 1e3 * (1 + 1e-9) && rad > 999);
  for (Pd p : ps) CHECK((p - ctr).dist() <= rad * (1 + 1e-7));
  puts("OK geometry/circles");
}
