#include "../testlib.hpp"
#include "../../content/geometry/point.hpp"

typedef __int128 LL;
typedef P<ll> Pl;
typedef P<double> Pd;

// Polar angle in [0, 2 pi) of the reduced direction, so equal
// directions get bitwise-equal keys. Exact enough for small input.
double angle_key(Pl v) {
  ll g = gcd(llabs(v.x), llabs(v.y));
  double a = atan2((double)(v.y / g), (double)(v.x / g));
  return a < 0 ? a + 2 * M_PI : a;
}

// Orientation of o -> a -> b from angles, without cross products.
int brute_turn(Pl o, Pl a, Pl b) {
  if (a == o || b == o) return 0;
  double s = sin(atan2(b.y - o.y, b.x - o.x) - atan2(a.y - o.y, a.x - o.x));
  return fabs(s) < 1e-9 ? 0 : s > 0 ? 1 : -1;
}

// Exact angle order by quadrant, then __int128 cross product.
int quadrant(Pl v) {
  if (v.x > 0 && v.y >= 0) return 0;
  if (v.x <= 0 && v.y > 0) return 1;
  if (v.x < 0 && v.y <= 0) return 2;
  return 3;
}
bool exact_less(Pl a, Pl b) {
  if (quadrant(a) != quadrant(b)) return quadrant(a) < quadrant(b);
  return (LL)a.x * b.y - (LL)a.y * b.x > 0;
}

Pl rnd_pt(ll c) { return {rnd(-c, c), rnd(-c, c)}; }
Pl rnd_nonzero(ll c) {
  Pl v;
  while (v == Pl{}) v = rnd_pt(c);
  return v;
}
bool near(double a, double b, double tol = 1e-9) {
  return fabs(a - b) <= tol * max(1.0, fabs(b));
}

int main() {
  rep(it, 0, iters(20000)) {
    Pl o = rnd_pt(20), a = rnd_pt(20), b = rnd_pt(20);
    ll k = rnd(-5, 5);
    CHECK((a + b == Pl{a.x + b.x, a.y + b.y}));
    CHECK((a - b == Pl{a.x - b.x, a.y - b.y}));
    CHECK((a * k == Pl{a.x * k, a.y * k}));
    CHECK_EQ(a == b, a.x == b.x && a.y == b.y);
    CHECK_EQ(a < b, a.x < b.x || (a.x == b.x && a.y < b.y));
    CHECK_EQ(a.dot(b), a.x * b.x + a.y * b.y);
    CHECK_EQ(a.dist2(), a.x * a.x + a.y * a.y);
    CHECK(near(a.dist(), hypot(a.x, a.y)));
    CHECK_EQ(sgn(o.cross(a, b)), brute_turn(o, a, b));
    CHECK_EQ(sgn(a.cross(b)), brute_turn(Pl{}, a, b));
    double sine =
        sin(atan2(b.y - o.y, b.x - o.x) - atan2(a.y - o.y, a.x - o.x));
    CHECK(near(fabs((double)o.cross(a, b)),
               fabs((a - o).dist() * (b - o).dist() * sine)));
    CHECK((a.perp() == Pl{-a.y, a.x}));
    CHECK_EQ(sgn(k), k > 0 ? 1 : k < 0 ? -1 : 0);
    if (a == Pl{} || b == Pl{}) continue;
    CHECK_EQ(angle_less(a, b), angle_key(a) < angle_key(b));
    CHECK_EQ(a.half(), angle_key(a) >= M_PI);
    CHECK(near(a.angle(), atan2(a.y, a.x)));
  }

  rep(it, 0, iters(5000)) {
    Pd a{rnd_real(-1e3, 1e3), rnd_real(-1e3, 1e3)};
    double t = rnd_real(-10, 10);
    Pd r = a.rotate(t);
    CHECK(near(r.dist(), a.dist()));
    CHECK(near(r.x, a.x * cos(t) - a.y * sin(t), 1e-6));
    CHECK(near(a.cross(r), a.dist2() * sin(t), 1e-6));
    CHECK(near(a.dot(r), a.dist2() * cos(t), 1e-6));
    Pd q = a.rotate(M_PI / 2);
    CHECK(near(q.x, a.perp().x, 1e-9) && near(q.y, a.perp().y, 1e-9));
    Pd u = a.unit();
    CHECK(near(u.dist(), 1) && u.dot(a) > 0);
    CHECK(fabs(u.cross(a)) < 1e-9 * a.dist());
    double ang = a.angle();
    CHECK(-M_PI < ang && ang <= M_PI);
    CHECK(near(a.dist() * cos(ang), a.x, 1e-9));
    CHECK(near((a / 4).x, a.x / 4) && near((a * 0.5).y, a.y / 2));
    P<long double> al{(long double)a.x, (long double)a.y};
    static_assert(is_same_v<decltype(al.dist()), long double>);
    CHECK(fabsl(al.rotate(t).dist() - al.dist()) < 1e-12L * al.dist());
  }

  // Extremes: |coordinates| <= 1e9 keeps cross products in ll.
  ll C = 1'000'000'000;
  rep(it, 0, iters(20000)) {
    Pl o = rnd_pt(C), a = rnd_pt(C), b = rnd_pt(C);
    if (it % 4 == 0) o = {-C, -C}, a = {C, -C}, b = {C, C};
    LL want = (LL)(a.x - o.x) * (b.y - o.y) - (LL)(a.y - o.y) * (b.x - o.x);
    CHECK(o.cross(a, b) == want);
    LL d2 = (LL)(a.x - b.x) * (a.x - b.x) + (LL)(a.y - b.y) * (a.y - b.y);
    CHECK((a - b).dist2() == d2);
    Pl u = a - o, w = b - o;
    if (u == Pl{} || w == Pl{}) continue;
    CHECK_EQ(angle_less(u, w), exact_less(u, w));
  }
  // Nearly parallel: cross = -1, beyond double atan2 resolution.
  Pl u{C, C - 1}, w{C - 1, C - 2};
  CHECK_EQ(u.cross(w), -1LL);
  CHECK(angle_less(w, u) && !angle_less(u, w));
  CHECK(!angle_less(Pl{2, 0}, Pl{5, 0}) && !angle_less(Pl{5, 0}, Pl{2, 0}));
  CHECK(angle_less(Pl{1, 0}, Pl{-1, 0}));
  CHECK(angle_less(Pl{-1, 1}, Pl{-1, 0}));
  CHECK(angle_less(Pl{-1, 0}, Pl{0, -1}));
  CHECK(angle_less(Pl{-1, 0}, Pl{1, -1}));
  // Negative zero y: atan2 gives -pi, the angle order treats it as pi.
  Pd nz = Pd{1, 0} * -1.0;
  CHECK(signbit(nz.y) && nz.angle() == -M_PI && nz.half());
  CHECK(!angle_less(nz, Pd{-1, 0}) && !angle_less(Pd{-1, 0}, nz));
  CHECK(angle_less(Pd{-1, 1}, nz) && angle_less(nz, Pd{-1, -1}));
  Pd pz = Pd{-1, 0} * -1.0; // (1, -0.0) is angle 0, not 2 pi
  CHECK(signbit(pz.y) && !pz.half() && angle_less(pz, Pd{1, 1}));

  int n = 100'000;
  vector<Pl> v(n);
  for (auto& p : v) p = rnd_nonzero(C);
  rep(i, 0, 100) v[i] = v[n - 1 - i] * 2; // repeated directions
  double start = elapsed();
  sort(all(v), angle_less<ll>);
  rep(i, 0, n - 1) CHECK(!exact_less(v[i + 1], v[i]));
  fprintf(stderr, "  sort 1e5 by angle: %.2fs\n", elapsed() - start);
  puts("OK geometry/point");
}
