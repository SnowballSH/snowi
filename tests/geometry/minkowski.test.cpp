#include "../testlib.hpp"
#include "../../content/geometry/minkowski.hpp"
#include "../../content/geometry/convex-hull.hpp"

// ASan's fake stacks make each by-value call of a small struct
// about 1000x slower on this toolchain; geometry code is all such
// calls and never returns pointers to locals.
extern "C" const char* __asan_default_options() {
  return "detect_stack_use_after_return=0";
}

typedef __int128 LL;
typedef P<ll> Pl;

Pl rnd_pt(ll c) { return {rnd(-c, c), rnd(-c, c)}; }
LL cr(Pl o, Pl a, Pl b) {
  return (LL)(a.x - o.x) * (b.y - o.y) - (LL)(a.y - o.y) * (b.x - o.x);
}
LL dot(Pl a, Pl d) { return (LL)a.x * d.x + (LL)a.y * d.y; }
bool on_seg(Pl a, Pl b, Pl p) { // closed segment, exact
  return cr(a, b, p) == 0 && min(a.x, b.x) <= p.x && p.x <= max(a.x, b.x) &&
         min(a.y, b.y) <= p.y && p.y <= max(a.y, b.y);
}

// The sum is hull(S) for S = {a + b}: its vertices lie in S, it is
// strictly convex and counter-clockwise from its smallest point, and
// contains all of S. One or two points if S is a point or collinear.
void check_sum(const vector<Pl>& a, const vector<Pl>& b) {
  vector<Pl> s = minkowski(a, b), all_sums;
  for (Pl p : a)
    for (Pl q : b) all_sums.push_back(p + q);
  set<Pl> S(all(all_sums));
  int n = sz(s);
  CHECK(n >= 1 && n <= sz(a) + sz(b));
  for (Pl p : s) CHECK(S.count(p));
  CHECK(s[0] == *S.begin());
  if (n == 1) CHECK_EQ(sz(S), 1);
  if (n == 2) {
    CHECK(s[0] != s[1]);
    for (Pl p : S) CHECK(on_seg(s[0], s[1], p));
  }
  if (n >= 3) {
    rep(i, 0, n) CHECK(cr(s[i], s[(i + 1) % n], s[(i + 2) % n]) > 0);
    rep(i, 1, n - 1) CHECK(cr(s[0], s[i], s[i + 1]) > 0); // winds once
    for (Pl p : S)
      rep(i, 0, n) CHECK(cr(s[i], s[(i + 1) % n], p) >= 0);
  }
}
// Support function: the extreme value in direction d adds up.
void check_support(const vector<Pl>& a, const vector<Pl>& b,
                   const vector<Pl>& s, Pl d) {
  auto best = [&](const vector<Pl>& v) {
    LL r = dot(v[0], d);
    for (Pl p : v) r = max(r, dot(p, d));
    return r;
  };
  CHECK(best(s) == best(a) + best(b));
}

vector<Pl> rnd_convex(int n, ll c) {
  vector<Pl> pts(n);
  for (auto& p : pts) p = rnd_pt(c);
  return convex_hull(pts);
}

int main() {
  rep(it, 0, iters(3000)) {
    ll C = it % 4 == 0 ? 999'999'997 : rnd(1, 6); // room for +-3
    vector<Pl> a = rnd_convex((int)rnd(1, 12), C);
    vector<Pl> b = rnd_convex((int)rnd(1, 12), C);
    if (it % 5 == 1) { // a parallel copy: every edge direction shared
      b = a;
      for (auto& p : b) p = p * (C < 10 ? rnd(1, 2) : 1) + rnd_pt(3);
      b = convex_hull(b);
    }
    rotate(a.begin(), a.begin() + rnd(0, sz(a) - 1), a.end());
    check_sum(a, b);
    vector<Pl> s = minkowski(a, b);
    rep(k, 0, 10) check_support(a, b, s, rnd_pt(1'000'000'000));
  }
  // Hand-made: the triangle + unit square of Packet 20 (area 5).
  vector<Pl> tri{{0, 0}, {2, 0}, {0, 1}}, sq{{0, 0}, {1, 0}, {1, 1}, {0, 1}};
  vector<Pl> s = minkowski(tri, sq);
  CHECK((s == vector<Pl>{{0, 0}, {3, 0}, {3, 1}, {1, 2}, {0, 2}}));
  CHECK((minkowski(vector<Pl>{{1, 2}}, vector<Pl>{{3, 4}}) ==
         vector<Pl>{{4, 6}}));
  CHECK((minkowski(vector<Pl>{{1, 1}}, tri) ==
         vector<Pl>{{1, 1}, {3, 1}, {1, 2}}));
  // Parallel segments, and crossing segments (a parallelogram).
  vector<Pl> s1{{0, 0}, {2, 1}}, s2{{5, 3}, {1, 1}};
  CHECK((minkowski(s1, s2) == vector<Pl>{{1, 1}, {7, 4}}));
  vector<Pl> s3{{0, 3}, {0, 0}};
  CHECK((minkowski(s1, s3) == vector<Pl>{{0, 0}, {2, 1}, {2, 4}, {0, 3}}));
  ll C = 1'000'000'000;
  vector<Pl> big{{-C, -C}, {C, -C}, {C, C}, {-C, C}};
  s = minkowski(big, big);
  CHECK((s == vector<Pl>{{-2 * C, -2 * C},
                         {2 * C, -2 * C},
                         {2 * C, 2 * C},
                         {-2 * C, 2 * C}}));

  // The Example: collision of A and B + t, distance between them.
  rep(it, 0, iters(300)) {
    ll D = it % 2 ? 500'000'000 : 10;
    vector<Pl> A = rnd_convex((int)rnd(1, 8), D),
               B = rnd_convex((int)rnd(1, 8), D), b = B;
    Pl t = rnd_pt(D);
    for (auto& p : b) p = p * -1;
    vector<Pl> d = minkowski(A, b);
    bool meet = in_hull(d, t) >= 0;
    vector<Pl> Bt;
    for (Pl p : B) Bt.push_back(p + t);
    bool want = false; // brute: a vertex inside the other, or edges cross
    for (Pl p : A) want |= in_hull(Bt, p) >= 0;
    for (Pl p : Bt) want |= in_hull(A, p) >= 0;
    rep(i, 0, sz(A))
      rep(j, 0, sz(Bt)) {
        Pl p = A[i], q = A[(i + 1) % sz(A)], u = Bt[j],
           v = Bt[(j + 1) % sz(Bt)];
        typedef P<long double> Q; // exact for these integers
        want |= !seg_inter(Q(p.x, p.y), Q(q.x, q.y), Q(u.x, u.y), Q(v.x, v.y))
                     .empty();
      }
    CHECK_EQ(meet, want);
    if (meet || it % 2) continue;
    // Distance between disjoint A and B + t: from t to A - B.
    vector<P<double>> f;
    for (auto p : d) f.push_back(P<double>(p.x - t.x, p.y - t.y));
    double dist = 1e18, brute = 1e18;
    rep(i, 0, sz(f))
      dist = min(dist, seg_dist(f[i], f[(i + 1) % sz(f)], P<double>{}));
    auto pd = [](Pl p) { return P<double>(p.x, p.y); };
    for (auto [X, Y] : {pair{A, Bt}, pair{Bt, A}})
      rep(i, 0, sz(X))
        for (Pl p : Y)
          brute = min(brute, seg_dist(pd(X[i]), pd(X[(i + 1) % sz(X)]), pd(p)));
    CHECK(fabs(dist - brute) < 1e-9 * max(1.0, brute));
  }

  // in_hull on a sum of polygons with |coordinates| <= 1e9 needs
  // P<__int128> (the note in the header).
  rep(it, 0, iters(300)) {
    vector<Pl> A = rnd_convex((int)rnd(1, 8), C),
               B = rnd_convex((int)rnd(1, 8), C), d = minkowski(A, B);
    vector<P<LL>> wide;
    for (Pl p : d) wide.push_back(P<LL>(p.x, p.y));
    rep(q, 0, 20) {
      Pl t = q % 2 ? d[rnd(0, sz(d) - 1)] : rnd_pt(2 * C);
      int n = sz(d), want = 1;
      if (n <= 2) want = on_seg(d[0], d[n - 1], t) ? 0 : -1;
      else
        rep(i, 0, n) {
          LL c = cr(d[i], d[(i + 1) % n], t);
          if (c < 0) want = -1;
          if (c == 0 && want == 1) want = 0;
        }
      CHECK_EQ(in_hull(wide, P<LL>(t.x, t.y)), want);
    }
  }

  // Large: two polygons with 5e4 vertices near circles of radius
  // 1e9 (sharing many edge directions), plus a support check.
  int n = 50'000;
  vector<Pl> pa(n), pb(n);
  rep(i, 0, n) {
    double u = 2 * M_PI * i / n, v = 2 * M_PI * (i + 0.5 * (i % 2)) / n;
    pa[i] = {llround(C * cos(u)), llround(C * sin(u))};
    pb[i] = {llround(C / 2 * cos(v)), llround(C / 2 * sin(v))};
  }
  pa = convex_hull(pa), pb = convex_hull(pb);
  double start = elapsed();
  s = minkowski(pa, pb);
  fprintf(stderr, "  sum of %d + %d vertices: %.2fs\n", sz(pa), sz(pb),
          elapsed() - start);
  set<pair<ll, ll>> dirs;
  for (auto* v : {&pa, &pb})
    rep(i, 0, sz(*v)) {
      Pl e = (*v)[(i + 1) % sz(*v)] - (*v)[i];
      ll g = gcd(e.x, e.y);
      dirs.insert({e.x / g, e.y / g});
    }
  CHECK_EQ(sz(s), sz(dirs));
  rep(i, 0, sz(s)) CHECK(cr(s[i], s[(i + 1) % sz(s)], s[(i + 2) % sz(s)]) > 0);
  rep(k, 0, 200) check_support(pa, pb, s, rnd_pt(C));
  puts("OK geometry/minkowski");
}
