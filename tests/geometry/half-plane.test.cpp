#include "../testlib.hpp"
#include "../../content/geometry/half-plane.hpp"

typedef __int128 LL;
typedef long double ld;
typedef P<ll> Pl;
typedef P<double> Pd;
typedef array<Pl, 2> Hl;

Pl rnd_pt(ll c) { return {rnd(-c, c), rnd(-c, c)}; }
Pd to_d(Pl p) { return Pd(p.x, p.y); }
LL cr(Pl a, Pl b) { return (LL)a.x * b.y - (LL)a.y * b.x; }

// Exact brute force. A vertex of the region is the meeting point of
// two boundary lines that satisfies every half-plane (closed); the
// region is the convex hull of those points.
struct Q { // the point (x / d, y / d), d > 0
  LL x, y, d;
  ld fx() const { return (ld)x / d; }
  ld fy() const { return (ld)y / d; }
};
bool meet(Hl a, Hl b, Q& q) {
  Pl d1 = a[1] - a[0], d2 = b[1] - b[0];
  LL den = cr(d1, d2), t = cr(b[0] - a[0], d2);
  if (den == 0) return false;
  q = {a[0].x * den + d1.x * t, a[0].y * den + d1.y * t, den};
  if (den < 0) q = {-q.x, -q.y, -den};
  // Reduce, so that points on the box at 1e9 keep small terms.
  LL g = gcd(gcd(q.x < 0 ? -q.x : q.x, q.y < 0 ? -q.y : q.y), q.d);
  q = {q.x / g, q.y / g, q.d / g};
  return true;
}
int side(Hl h, Q q) { // sign of the point against the line
  Pl d = h[1] - h[0];
  LL v = d.x * (q.y - h[0].y * q.d) - d.y * (q.x - h[0].x * q.d);
  return (v > 0) - (v < 0);
}
bool less_q(Q a, Q b) {
  LL u = a.x * b.d - b.x * a.d;
  return u < 0 || (u == 0 && a.y * b.d < b.y * a.d);
}
bool same_q(Q a, Q b) { return !less_q(a, b) && !less_q(b, a); }
int orient(Q o, Q a, Q b) {
  LL ax = a.x * o.d - o.x * a.d, ay = a.y * o.d - o.y * a.d;
  LL bx = b.x * o.d - o.x * b.d, by = b.y * o.d - o.y * b.d;
  LL v = ax * by - ay * bx;
  return (v > 0) - (v < 0);
}
vector<Q> exact_region(vector<Hl> hs, ll big) {
  Pl c[] = {{-big, -big}, {big, -big}, {big, big}, {-big, big}};
  rep(i, 0, 4) hs.push_back({c[i], c[(i + 1) % 4]});
  vector<Q> pts;
  rep(i, 0, sz(hs))
    rep(j, 0, i) {
      Q q;
      if (!meet(hs[i], hs[j], q)) continue;
      bool ok = true;
      for (Hl h : hs) ok &= side(h, q) >= 0;
      if (ok) pts.push_back(q);
    }
  sort(all(pts), less_q);
  pts.erase(unique(all(pts), same_q), pts.end());
  if (sz(pts) < 3) return {};
  vector<Q> h(2 * sz(pts));
  int k = 0;
  rep(it, 0, 2) { // monotone chain without collinear points
    for (int s = k; Q p : pts) {
      while (k >= s + 2 && orient(h[k - 2], h[k - 1], p) <= 0) k--;
      h[k++] = p;
    }
    k--, reverse(all(pts));
  }
  h.resize(k);
  return sz(h) >= 3 ? h : vector<Q>{}; // collinear: no interior
}

vector<array<Pd, 2>> to_d(const vector<Hl>& hs) {
  vector<array<Pd, 2>> r;
  for (auto [s, e] : hs) r.push_back({to_d(s), to_d(e)});
  return r;
}
// Same vertices, both counter-clockwise, up to the starting vertex.
// Returns the number of vertices.
int same_region(const vector<Pd>& got, const vector<Q>& want) {
  CHECK_EQ(sz(got), sz(want));
  int n = sz(got), k = 0;
  if (!n) return 0;
  auto d0 = [&](Pd p) { return hypot(p.x - want[0].fx(), p.y - want[0].fy()); };
  rep(i, 0, n)
    if (d0(got[i]) < d0(got[k])) k = i;
  rep(i, 0, n) {
    Pd g = got[(k + i) % n];
    CHECK(fabsl(g.x - want[i].fx()) < 1e-6 && fabsl(g.y - want[i].fy()) < 1e-6);
  }
  return n;
}
// With own_box the box is passed as half-planes, inside the default
// one.
int check_region(vector<Hl> hs, ll big, bool own_box = false) {
  vector<Q> want = exact_region(hs, big);
  Pl c[] = {{-big, -big}, {big, -big}, {big, big}, {-big, big}};
  if (own_box)
    rep(i, 0, 4) hs.push_back({c[i], c[(i + 1) % 4]});
  return same_region(
      own_box ? half_plane(to_d(hs)) : half_plane(to_d(hs), (double)big), want);
}

vector<Hl> random_planes(int n, ll C) {
  Pl o = rnd_pt(C); // most half-planes contain o: non-empty results
  vector<Hl> hs;
  while (sz(hs) < n) {
    Pl s = rnd_pt(C), e = rnd_pt(C);
    if (s == e) continue;
    if (rnd(0, 4) && s.cross(e, o) < 0) swap(s, e);
    hs.push_back({s, e});
  }
  return hs;
}

int main() {
  int nonempty = 0, iterations = iters(4000);
  rep(it, 0, iterations) {
    int n = (int)rnd(0, 9);
    ll C = rnd(1, 5), big = rnd(C, 3 * C);
    vector<Hl> hs = random_planes(n, C);
    if (it % 5 == 1 && n) // repeated and reversed half-planes
      rep(i, 0, 3) {
        Hl h = hs[rnd(0, sz(hs) - 1)];
        if (rnd(0, 1)) swap(h[0], h[1]);
        if (rnd(0, 1)) h[1] = h[0] + (h[1] - h[0]) * 2;
        hs.push_back(h);
      }
    if (it % 5 == 2) { // all lines through one point: a cone or a point
      Pl o = rnd_pt(C);
      for (auto& h : hs) h[1] = h[0] == o ? h[0] + Pl{1, 1} : o;
    }
    if (it % 5 == 3) // few directions: many parallel lines, strips
      for (auto& h : hs)
        h[1] = h[0] + Pl{rnd(0, 1), rnd(0, 1) ? 1 : -1} * rnd(-2, 2);
    hs.erase(remove_if(all(hs), [](Hl h) { return h[0] == h[1]; }), hs.end());
    nonempty += check_region(hs, big, it % 3 == 0) > 0;
  }
  CHECK(nonempty > iterations / 3 && nonempty < iterations * 9 / 10);

  // The default box at 1e9, where doubles resolve only about 1e-7.
  // Lines y = x and y = -x pass through its corners; their points,
  // integer or decimal, must not split a corner into two vertices.
  rep(it, 0, iters(2000)) {
    vector<Hl> hs = random_planes((int)rnd(0, 4), 5);
    vector<array<Pd, 2>> hd = to_d(hs);
    rep(k, 0, rnd(1, 3)) {
      ll a = rnd(-5000, 5000), b = rnd(-5000, 5000), sg = rnd(0, 1) ? 1 : -1;
      if (a == b) continue;
      ll dir = a < b ? 1 : -1;
      hs.push_back({Pl{0, 0}, Pl{dir, sg * dir}}); // the same line
      double sc = it % 2 ? 1 : 1000;
      hd.push_back({Pd{a / sc, sg * (a / sc)}, Pd{b / sc, sg * (b / sc)}});
    }
    same_region(half_plane(hd), exact_region(hs, 1'000'000'000));
  }

  // Hand-made cases.
  Pl o{0, 0}, X{1, 0}, Y{0, 1};
  vector<Pd> sq = half_plane(to_d(vector<Hl>{}), 1.0);
  CHECK((sq == vector<Pd>{{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}));
  check_region({{o, X}, {Y, Y - X}}, 3);                 // strip 0 <= y <= 1
  check_region({{o, X}, {o, o - X}}, 3);                 // line y = 0: empty
  check_region({{Y, Y + X}, {o, o - X}}, 3);             // y >= 1, y <= 0
  check_region({{o, X}, {X, X + Y}, {o + Y, o}}, 3);     // box-bounded strip
  check_region({{o, X}, {X, Y}, {Y, o}, {o, X + Y}}, 3); // 2 regions meet
  // A point: x >= 0, x <= 0, y >= 0, y <= 0.
  check_region({{o, o - Y}, {o, Y}, {o, X}, {o, o - X}}, 3);
  // Unbounded with the default box: the half-plane y >= 0.
  vector<Pd> up = half_plane(vector<array<Pd, 2>>{{Pd{0, 0}, Pd{1, 0}}});
  CHECK((up == vector<Pd>{{-1e9, 0}, {1e9, 0}, {1e9, 1e9}, {-1e9, 1e9}}));
  // x + y >= 0 cuts the default box through two corners: a triangle.
  CHECK_EQ(sz(half_plane(vector<array<Pd, 2>>{{Pd{-3, 3}, Pd{-1, 1}}})), 3);
  CHECK_EQ(sz(half_plane(
               vector<array<Pd, 2>>{{Pd{3.48, 3.48}, Pd{-3.141, -3.141}}})),
           3);
  // Kernel of an L-shaped polygon (counter-clockwise).
  vector<Pd> L{{0, 0}, {4, 0}, {4, 1}, {1, 1}, {1, 4}, {0, 4}};
  vector<array<Pd, 2>> hs;
  rep(i, 0, 6) hs.push_back({L[i], L[(i + 1) % 6]});
  vector<Pd> ker = half_plane(hs);
  CHECK((ker == vector<Pd>{{0, 0}, {1, 0}, {1, 1}, {0, 1}}));

  // A direction with y = -0.0 (angle -pi in atan2) and a parallel
  // one with y = 0.0: -1 <= y <= 0 inside the box.
  vector<array<Pd, 2>> z{
      {Pd{1, 0.0}, Pd{0, -0.0}}, {Pd{1, 1}, Pd{0, 1}}, {Pd{0, -1}, Pd{1, -1}}};
  vector<Pd> rect = half_plane(z, 3.0);
  CHECK((rect == vector<Pd>{{-3, -1}, {3, -1}, {3, 0}, {-3, 0}}));

  // Large: 1e5 tangents of a circle at jittered angles (gaps far
  // above eps, so all of them are edges), area from the angles; then
  // random half-planes that contain the disc.
  int n = 100'000;
  double R = 1000;
  vector<double> ang(n);
  rep(i, 0, n) ang[i] = 2 * M_PI * (i + rnd_real(0.1, 0.9)) / n;
  hs.clear();
  for (double a : ang) {
    Pd t{R * cos(a), R * sin(a)};
    hs.push_back({t, t + Pd{-sin(a), cos(a)}});
  }
  shuffle(all(hs), test_rng);
  double start = elapsed();
  vector<Pd> poly = half_plane(hs);
  double area = 0, want = 0;
  rep(i, 0, sz(poly)) area += poly[i].cross(poly[(i + 1) % sz(poly)]) / 2;
  rep(i, 0, n)
    want +=
        R * R * tan(fmod(ang[(i + 1) % n] - ang[i] + 2 * M_PI, 2 * M_PI) / 2);
  CHECK_EQ(sz(poly), n);
  CHECK(fabs(area - want) < 1e-6 * want);
  hs.clear();
  rep(i, 0, n) {
    double a = rnd_real(0, 2 * M_PI), r = rnd_real(R, 3 * R);
    Pd t{r * cos(a), r * sin(a)};
    hs.push_back({t, t + Pd{-sin(a), cos(a)}});
  }
  poly = half_plane(hs);
  fprintf(stderr, "  2 x half-plane intersection of 1e5: %.2fs\n",
          elapsed() - start);
  area = 0;
  rep(i, 0, sz(poly)) {
    Pd p = poly[i], q = poly[(i + 1) % sz(poly)];
    area += p.cross(q) / 2;
    CHECK(p.cross(q) > 0 && (q - p).cross(poly[(i + 2) % sz(poly)] - q) > 0);
  }
  CHECK(area > M_PI * R * R);
  for (Pd p : poly)
    for (auto [s, e] : hs) CHECK(line_dist(s, e, p) > -1e-6);
  puts("OK geometry/half-plane");
}
