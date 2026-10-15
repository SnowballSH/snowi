#include "../testlib.hpp"
#include "../../content/misc/ternary-search.hpp"

// Strictly decreasing up to a random minimum, then increasing;
// with `flat`, the part after the minimum may have flat steps.
vector<ll> valley(int n, bool flat) {
  int m = (int)rnd(0, n - 1);
  vector<ll> f(n);
  f[m] = rnd(-1e12, 1e12);
  for (int i = m - 1; i >= 0; i--) f[i] = f[i + 1] + rnd(1, 5);
  rep(i, m + 1, n) f[i] = f[i - 1] + rnd(!flat, 5);
  return f;
}

// Convex: sorted differences, so flat runs only at the bottom.
vector<ll> convex(int n) {
  vector<ll> d = rnd_vec(n - 1, -6, 6), f(n, rnd(-99, 99));
  sort(all(d));
  rep(i, 1, n) f[i] = f[i - 1] + d[i - 1];
  return f;
}

void check_int(const vector<ll>& f) {
  int n = sz(f), calls = 0;
  ll lo = rnd(-1e9, 1e9);
  ll got = ternary(lo, lo + n - 1, [&](ll x) {
    CHECK(lo <= x && x < lo + n);
    CHECK(++calls <= 2 * __lg(2 * n - 1)); // 2 ceil(log2 n)
    return f[x - lo];
  });
  CHECK_EQ(got - lo, (ll)(min_element(all(f)) - f.begin()));
}

// golden must find a minimiser up to what doubles resolve: x is
// within a few ulps of the minimisers [cl, cr], or f(x) is within
// 1e-9 of the minimum (a smooth minimum is flat).
template<class F>
double check_real(double a, double b, F f, double cl, double cr) {
  int calls = 0;
  double x = golden(a, b, [&](double t) {
    calls++;
    return f(t);
  });
  CHECK_EQ(calls, 102);
  double d = 1e-13 * max({1.0, abs(a), abs(b)}), fmin = f(cl);
  CHECK(a - d <= x && x <= b + d);
  CHECK((cl - d <= x && x <= cr + d) ||
        f(x) <= fmin + 1e-9 * max(1.0, abs(fmin)));
  return x;
}

void random_real() {
  double a = rnd_real(-1e3, 1e3), b = a + rnd_real(0, 1e3);
  if (rnd(0, 3) == 0) a *= 1e6, b *= 1e6;
  double c = rnd_real(a, b), s = rnd_real(1e-3, 1e3);
  double k = rnd_real(-1e6, 1e6), tol = 1e-9 * max(1.0, b - a);
  switch (rnd(0, 5)) {
  case 0:
    check_real(a, b, [&](double x) { return s * (x - c) * (x - c) + k; }, c, c);
    break;
  case 1: {
    double x = check_real(a, b, [&](double t) { return s * abs(t - c); }, c, c);
    CHECK(abs(x - c) <= tol);
    break;
  }
  case 2: // asymmetric: linear on the left, cubic on the right
    check_real(
        a, b,
        [&](double x) { return x < c ? s * (c - x) + k : pow(x - c, 3) + k; },
        c, c);
    break;
  case 3: { // flat bottom of half-width w
    double w = rnd_real(0, (b - a) / 2);
    double x = check_real(
        a, b, [&](double t) { return max(0.0, abs(t - c) - w); }, max(a, c - w),
        min(b, c + w));
    CHECK(abs(x - c) <= w + tol);
    break;
  }
  case 4: // monotone: minimum at an endpoint
    CHECK(check_real(a, b, [&](double x) { return s * x; }, a, a) - a <= tol);
    CHECK(b - check_real(a, b, [&](double x) { return -s * x; }, b, b) <= tol);
    break;
  case 5: // smooth, not a polynomial
    b = a + rnd_real(0, 40), c = rnd_real(a, b);
    check_real(a, b, [&](double x) { return exp(x - c) + exp(c - x); }, c, c);
  }
}

int main() {
  rep(it, 0, iters(3000)) {
    int n = (int)rnd(1, 40);
    check_int(valley(n, false));
    check_int(valley(n, true));
    check_int(convex(n));
    random_real();
  }
  // Edges: one point, a constant function, a huge range.
  CHECK_EQ(ternary(5, 5, [](ll) { return 0; }), 5LL);
  CHECK_EQ(ternary(-3, 9, [](ll) { return 7; }), -3LL);
  rep(it, 0, 100) {
    ll c = rnd(-1e18, 1e18);
    CHECK_EQ(ternary(-1e18, 1e18, [&](ll x) { return abs(x - c); }), c);
    CHECK_EQ(ternary(-1e18, 1e18, [&](ll x) { return pair(abs(x - c), 0); }),
             c);
  }
  CHECK_EQ(golden(2.5, 2.5, [](double x) { return x * x; }), 2.5);
  // golden, like ternary, accepts an f that returns a pair.
  rep(it, 0, 100) {
    double c = rnd_real(-1e3, 1e3);
    double x = golden(-1e3, 1e3, [&](double t) { return pair(abs(t - c), 0); });
    CHECK(abs(x - c) < 1e-9);
  }
  // Speed: 1e5 integer searches over 1e9 values, 1e4 real ones.
  ll sum = 0;
  rep(it, 0, 100'000) {
    ll c = rnd(0, 1e9);
    sum += ternary(0, 1e9, [&](ll x) { return (x - c) * (x - c); }) - c;
  }
  CHECK_EQ(sum, 0LL);
  rep(it, 0, 10'000) {
    double c = rnd_real(-1e3, 1e3);
    double x = golden(-1e3, 1e3, [&](double t) { return abs(t - c); });
    CHECK(abs(x - c) < 1e-9);
  }
  puts("OK misc/ternary-search");
}
