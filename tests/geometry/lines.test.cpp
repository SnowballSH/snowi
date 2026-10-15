#include "../testlib.hpp"
#include "../../content/geometry/lines.hpp"

// ASan's fake stacks make each by-value call of a small struct
// about 1000x slower on this toolchain; geometry code is all such
// calls and never returns pointers to locals.
extern "C" const char* __asan_default_options() {
  return "detect_stack_use_after_return=0";
}

typedef __int128 LL;
typedef P<ll> Pl;
typedef P<double> Pd;
typedef P<long double> PL;

Pd to_d(Pl p) { return Pd(p.x, p.y); }
PL to_ld(Pl p) { return PL(p.x, p.y); }
Pl rnd_pt(ll c) { return {rnd(-c, c), rnd(-c, c)}; }
bool near(double a, double b, double tol = 1e-9) {
  return fabs(a - b) <= tol * max(1.0, fabs(b));
}

// Orientation from angles, without cross products.
int brute_side(Pl s, Pl e, Pl p) {
  if (p == s) return 0;
  double t = sin(atan2(p.y - s.y, p.x - s.x) - atan2(e.y - s.y, e.x - s.x));
  return fabs(t) < 1e-9 ? 0 : t > 0 ? 1 : -1;
}
// p on [s, e] iff |sp| + |pe| = |se| (small coordinates).
bool brute_on(Pl s, Pl e, Pl p) {
  return fabs((p - s).dist() + (e - p).dist() - (e - s).dist()) < 1e-9;
}

// Exact rational point xn / den, yn / den with den > 0.
struct Q {
  LL xn, yn, den;
};
// Line through a and b as A x + B y = C.
array<LL, 3> eq(Pl a, Pl b) {
  LL A = b.y - a.y, B = a.x - b.x;
  return {A, B, A * a.x + B * a.y};
}
// 1 and the point, 0 parallel, -1 same line (Cramer's rule).
pair<int, Q> brute_line(Pl a, Pl b, Pl c, Pl d) {
  auto [A1, B1, C1] = eq(a, b);
  auto [A2, B2, C2] = eq(c, d);
  LL det = A1 * B2 - A2 * B1;
  if (det == 0) {
    bool same = A1 * C2 == A2 * C1 && B1 * C2 == B2 * C1;
    return {same ? -1 : 0, {}};
  }
  Q q{C1 * B2 - C2 * B1, A1 * C2 - A2 * C1, det};
  if (det < 0) q = {-q.xn, -q.yn, -det};
  return {1, q};
}
bool in_box(Q q, Pl a, Pl b) {
  return min(a.x, b.x) * q.den <= q.xn && q.xn <= max(a.x, b.x) * q.den &&
         min(a.y, b.y) * q.den <= q.yn && q.yn <= max(a.y, b.y) * q.den;
}
bool exact_on(Pl s, Pl e, Pl p) {
  if (s == e) return p == s;
  auto [A, B, C] = eq(s, e);
  return A * p.x + B * p.y == C && in_box(Q{p.x, p.y, 1}, s, e);
}
// Expected seg_inter: rational points, sorted when exact.
vector<Q> brute_seg(Pl a, Pl b, Pl c, Pl d) {
  if (a != b && c != d) {
    auto [kind, q] = brute_line(a, b, c, d);
    if (kind == 1) {
      if (in_box(q, a, b) && in_box(q, c, d)) return {q};
      return {};
    }
  }
  set<Pl> s;
  for (Pl p : {a, b, c, d})
    if (exact_on(a, b, p) && exact_on(c, d, p)) s.insert(p);
  vector<Q> res;
  for (Pl p : s) res.push_back({p.x, p.y, 1});
  return res;
}
template<class T> bool same_pt(P<T> p, Q q) {
  long double x = (long double)q.xn / q.den;
  long double y = (long double)q.yn / q.den;
  auto close = [](long double u, long double v) {
    return fabsl(u - v) <= 1e-12L * max(1.0L, fabsl(v));
  };
  return close(p.x, x) && close(p.y, y);
}
template<class T> void check_seg(Pl a, Pl b, Pl c, Pl d) {
  auto conv = [](Pl p) { return P<T>(p.x, p.y); };
  vector<P<T>> got = seg_inter(conv(a), conv(b), conv(c), conv(d));
  vector<Q> want = brute_seg(a, b, c, d);
  CHECK_EQ(sz(got), sz(want));
  rep(i, 0, sz(got)) CHECK(same_pt(got[i], want[i]));
}

// Minimise a convex function of t on [lo, hi].
template<class F> double ternary(double lo, double hi, F f) {
  rep(it, 0, 200) {
    double m1 = lo + (hi - lo) / 3, m2 = hi - (hi - lo) / 3;
    if (f(m1) < f(m2)) hi = m2;
    else lo = m1;
  }
  return f((lo + hi) / 2);
}

void check_all(Pl a, Pl b, Pl c, Pl d, Pl p) {
  if (a != b) CHECK_EQ(side_of(a, b, p), brute_side(a, b, p));
  CHECK_EQ(on_segment(a, b, p), brute_on(a, b, p));
  CHECK_EQ(on_segment(a, b, p), exact_on(a, b, p));
  check_seg<double>(a, b, c, d);
  check_seg<ll>(a, b, a, b);
  Pd A = to_d(a), B = to_d(b), Pp = to_d(p);
  auto dist_at = [&](double t) { return (A + (B - A) * t - Pp).dist(); };
  CHECK(near(seg_dist(a, b, p), ternary(0, 1, dist_at), 1e-7));
  CHECK(near(seg_dist(A, B, Pp), seg_dist(a, b, p)));
  if (a == b) return;
  double ld = line_dist(a, b, p);
  CHECK(near(fabs(ld), ternary(-200, 200, dist_at), 1e-7));
  CHECK_EQ(ld > 1e-9 ? 1 : ld < -1e-9 ? -1 : 0, brute_side(a, b, p));
  Pd q = line_proj(A, B, Pp), r = line_proj(A, B, Pp, true);
  CHECK(fabs((B - A).cross(q - A)) < 1e-9);
  CHECK(fabs((B - A).dot(Pp - q)) < 1e-9);
  CHECK(near((Pp - q).dist(), fabs(ld)));
  CHECK(near(r.x, 2 * q.x - Pp.x) && near(r.y, 2 * q.y - Pp.y));
  CHECK(near((r - A).dist(), (Pp - A).dist()));
  if (c == d) return;
  auto [kind, x] = line_inter(A, B, to_d(c), to_d(d));
  auto [bk, bq] = brute_line(a, b, c, d);
  CHECK_EQ(kind, bk);
  if (kind == 1) CHECK(same_pt(x, bq));
}

int main() {
  rep(it, 0, iters(30000)) {
    ll C = it % 3 ? 4 : 30; // small grids force degenerate cases
    Pl a = rnd_pt(C), b = rnd_pt(C), c = rnd_pt(C), d = rnd_pt(C),
       p = rnd_pt(C);
    if (it % 5 == 0) c = a + (b - a) * rnd(-2, 2), d = a + (b - a) * rnd(-2, 2);
    if (it % 7 == 0) p = a + (b - a) * rnd(-1, 2);
    check_all(a, b, c, d, p);
  }
  // Hand-made cases: cross, T-junction, touching ends, overlap,
  // disjoint collinear, identical, parallel, point segments.
  Pl o{0, 0}, x2{2, 0}, x4{4, 0}, y2{0, 2}, m{1, 1}, n{2, 2};
  CHECK_EQ(sz(seg_inter(to_d(o), to_d(n), to_d(x2), to_d(y2))), 1);
  CHECK(
      same_pt(seg_inter(to_d(o), to_d(n), to_d(x2), to_d(y2))[0], Q{1, 1, 1}));
  CHECK(
      (seg_inter(to_d(o), to_d(x4), to_d(x2), to_d(n)) == vector<Pd>{{2, 0}}));
  CHECK(
      (seg_inter(to_d(o), to_d(x2), to_d(x2), to_d(x4)) == vector<Pd>{{2, 0}}));
  CHECK((seg_inter(to_d(o), to_d(x4), to_d(Pl{3, 0}), to_d(Pl{1, 0})) ==
         vector<Pd>{{1, 0}, {3, 0}}));
  CHECK(seg_inter(to_d(o), to_d(Pl{1, 0}), to_d(x2), to_d(x4)).empty());
  CHECK((seg_inter(to_d(o), to_d(x4), to_d(x4), to_d(o)) ==
         vector<Pd>{{0, 0}, {4, 0}}));
  CHECK(seg_inter(to_d(o), to_d(x2), to_d(y2), to_d(n)).empty());
  CHECK((seg_inter(to_d(m), to_d(m), to_d(o), to_d(n)) == vector<Pd>{{1, 1}}));
  CHECK(seg_inter(to_d(m), to_d(m), to_d(o), to_d(x2)).empty());
  CHECK_EQ(line_inter(to_d(o), to_d(x2), to_d(y2), to_d(n)).first, 0);
  CHECK_EQ(line_inter(to_d(o), to_d(x2), to_d(x4), to_d(o)).first, -1);
  CHECK(on_segment(o, o, o) && !on_segment(o, o, m));
  CHECK_EQ(seg_dist(o, o, Pl{3, 4}), 5.0);
  CHECK_EQ(line_dist(o, x2, Pl{7, -3}), -3.0);

  // Extremes: |coordinates| <= 1e9. The cases taken stay exact in ll
  // and in long double; the brute force runs in __int128.
  ll C = 1'000'000'000;
  rep(it, 0, iters(20000)) {
    Pl a = rnd_pt(C), b = rnd_pt(C), c = rnd_pt(C), d = rnd_pt(C);
    if (it % 2) { // collinear: overlapping, touching or disjoint
      Pl v = rnd_pt(1000);
      a = rnd_pt(C - (ll)2e8), b = a + v * rnd(-1e5, 1e5);
      c = a + v * rnd(-1e5, 1e5), d = a + v * rnd(-1e5, 1e5);
    }
    if (it % 3 == 0) c = {a.x, b.y}, d = {b.x, a.y}; // cross
    if (it % 5 == 0) c = (a + b) / 2, d = c + (b - a).perp() / 3;
    check_seg<long double>(a, b, c, d);
    Pl p = it % 2 ? a + (b - a) / 2 : c;
    CHECK_EQ(on_segment(a, b, p), exact_on(a, b, p));
    LL cr = (LL)(b.x - a.x) * (p.y - a.y) - (LL)(b.y - a.y) * (p.x - a.x);
    CHECK_EQ(side_of(a, b, p), cr > 0 ? 1 : cr < 0 ? -1 : 0);
  }
  // line_inter in P<double> at 1e9 with the point far from s1 compared
  // with |e1 - s1|: the weights of s1 and e1 are huge and opposite.
  // Reproducer was off by about 3600; the rest is rounding of cross
  // products near 4e18 (measured worst 4e-5 over 1e6 such cases).
  auto far_error = [&](Pl s1, Pl e1, Pl s2, Pl e2) {
    auto [bk, bq] = brute_line(s1, e1, s2, e2);
    auto [kind, x] = line_inter(to_d(s1), to_d(e1), to_d(s2), to_d(e2));
    CHECK_EQ(kind, bk);
    long double ex = (long double)bq.xn / bq.den;
    long double ey = (long double)bq.yn / bq.den;
    return bk == 1 ? max(fabsl(x.x - ex), fabsl(x.y - ey)) : 0.0L;
  };
  CHECK(far_error({999999655, 999999162}, {999999658, 999999165},
                  {-810232897, -796596951}, {731982243, 720915984}) < 1e-4);
  rep(it, 0, iters(30000)) {
    Pl s1 = rnd_pt(C), e1 = s1 + rnd_pt(10), s2 = rnd_pt(C), e2 = rnd_pt(C);
    if (it % 2) e1 = rnd_pt(C);
    if (s1 == e1 || s2 == e2) continue;
    auto [bk, bq] = brute_line(s1, e1, s2, e2);
    if (bk == 1 && !in_box(bq, Pl{-C, -C}, Pl{C, C})) continue;
    CHECK(far_error(s1, e1, s2, e2) < 1e-12 * C);
  }
  CHECK(on_segment(Pl{-C, -C}, Pl{C, C}, Pl{C - 1, C - 1}));
  CHECK(!on_segment(Pl{-C, -C}, Pl{C, C - 1}, Pl{C - 1, C - 2}));
  CHECK_EQ(side_of(Pl{-C, -C}, Pl{C, C - 1}, Pl{C - 1, C - 2}), -1);

  double start = elapsed(), sum = 0;
  rep(it, 0, 100'000) {
    Pd a{rnd_real(-1, 1), rnd_real(-1, 1)}, b{rnd_real(-1, 1), 0},
        c{0, rnd_real(-1, 1)}, d{rnd_real(-1, 1), rnd_real(-1, 1)};
    sum += sz(seg_inter(a, b, c, d)) + seg_dist(a, b, c);
    sum += line_proj(a, b, d).x + line_inter(a, b, c, d).second.y;
  }
  fprintf(stderr, "  4e5 line ops: %.2fs (%g)\n", elapsed() - start, sum);
  puts("OK geometry/lines");
}
