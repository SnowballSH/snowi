#include "../testlib.hpp"
#include "../../content/linalg/gauss-real.hpp"

using Mat = vector<vector<ll>>;
using MatD = vector<vector<double>>;

ll perm_det(const Mat& a) {
  int n = sz(a);
  vi p(n);
  iota(all(p), 0);
  ll res = 0;
  do {
    ll term = 1;
    int inversions = 0;
    rep(i, 0, n) {
      term *= a[i][p[i]];
      rep(j, i + 1, n) inversions += p[i] > p[j];
    }
    res += inversions % 2 ? -term : term;
  } while (next_permutation(all(p)));
  return res;
}

Mat sub(const Mat& a, int rows, int cols) {
  Mat s;
  rep(i, 0, sz(a))
    if (rows >> i & 1) {
      s.emplace_back();
      rep(j, 0, sz(a[i]))
        if (cols >> j & 1) s.back().push_back(a[i][j]);
    }
  return s;
}

// Largest k with a nonzero k x k minor (exact).
int minor_rank(const Mat& a, int m) {
  int n = sz(a), best = 0;
  rep(rows, 0, 1 << n)
    rep(cols, 0, 1 << m) {
      int k = __builtin_popcount(rows);
      if (k != __builtin_popcount(cols) || k <= best) continue;
      if (perm_det(sub(a, rows, cols))) best = k;
    }
  return best;
}

// Exact unique solution by Cramer's rule on m independent rows.
vector<double> cramer(const Mat& a, const vector<ll>& b, int m) {
  int n = sz(a);
  rep(rows, 0, 1 << n) {
    if (__builtin_popcount(rows) != m) continue;
    Mat s = sub(a, rows, (1 << m) - 1);
    ll d = perm_det(s);
    if (!d) continue;
    vector<ll> bs;
    rep(i, 0, n)
      if (rows >> i & 1) bs.push_back(b[i]);
    vector<double> x(m);
    rep(j, 0, m) {
      Mat t = s;
      rep(i, 0, m) t[i][j] = bs[i];
      x[j] = (double)perm_det(t) / d;
    }
    return x;
  }
  return {};
}

MatD to_double(const Mat& a) {
  MatD d;
  for (auto& row : a) d.emplace_back(all(row));
  return d;
}

int main() {
  rep(it, 0, iters(4000)) {
    int n = (int)rnd(0, 4), m = (int)rnd(0, 4), hi = (int)rnd(1, 4);
    Mat a(n, vector<ll>(m));
    for (auto& row : a)
      for (ll& v : row) v = rnd(-hi, hi);
    if (n > 1 && rnd(0, 2) == 0) { // dependent row
      int i = (int)rnd(0, n - 1), j = (int)rnd(0, n - 1);
      ll f = rnd(-2, 2);
      rep(c, 0, m) a[i][c] = f * a[j][c];
    }
    MatD ad = to_double(a);
    auto [piv, d] = gauss_real(ad, m);
    CHECK_EQ(sz(piv), minor_rank(a, m));
    if (n == m) {
      ll want = perm_det(a);
      CHECK(abs(d - want) < 1e-9 * max(1.0, (double)abs(want)));
      CHECK(abs(det_real(to_double(a)) - want) < 1e-6);
      MatD tiny = to_double(a); // entries ~1e-10 need a smaller eps
      for (auto& row : tiny)
        for (double& v : row) v *= 1e-10;
      double scale = pow(1e-10, n), dt = det_real(tiny, 1e-20);
      CHECK(abs(dt - want * scale) <
            1e-9 * max(1.0, (double)abs(want)) * scale);
    }
    vector<ll> b(n), x0(m);
    for (ll& v : x0) v = rnd(-3, 3);
    rep(i, 0, n) {
      rep(j, 0, m) b[i] += a[i][j] * x0[j];
      if (rnd(0, 2) == 0) b[i] += rnd(-2, 2);
    }
    Mat ab = a;
    rep(i, 0, n) ab[i].push_back(b[i]);
    int ra = minor_rank(a, m), rab = minor_rank(ab, m + 1);
    vector<double> x(m, 7);
    int kind = solve_real(to_double(a), {all(b)}, x);
    CHECK_EQ(kind, ra < rab ? 0 : ra < m ? 2 : 1);
    if (kind) {
      CHECK_EQ(sz(x), m);
      rep(i, 0, n) {
        double s = 0;
        rep(j, 0, m) s += a[i][j] * x[j];
        CHECK(abs(s - b[i]) < 1e-9);
      }
    }
    if (kind == 1) {
      vector<double> want = cramer(a, b, m);
      rep(j, 0, m) CHECK(abs(x[j] - want[j]) < 1e-9);
    }
  }
  vector<double> x(2, 7);
  CHECK_EQ(solve_real({}, {}, x), 2);
  CHECK(x == vector<double>(2, 0));
  x.clear();
  CHECK_EQ(solve_real({{}, {}}, {0, 1}, x), 0);
  CHECK_EQ(det_real({}), 1.0);
  CHECK_EQ(det_real({{0, 1}, {1, 0}}), -1.0);
  CHECK_EQ(det_real({{1e-10, 0}, {0, 1}}), 0.0); // below eps
  CHECK_EQ(det_real({{1e-10, 0}, {0, 1}}, 1e-12), 1e-10);
  { // The Usage lines in order: only the last one, gauss_real, is
    // in place, so the others all see the original matrix.
    MatD u = {{2, 1}, {1, 3}};
    vector<double> ub = {3, 4}, ux(2);
    CHECK_EQ(solve_real(u, ub, ux), 1);
    CHECK(abs(ux[0] - 1) < 1e-12 && abs(ux[1] - 1) < 1e-12);
    CHECK(abs(det_real(u) - 5) < 1e-12);
    vi piv = gauss_real(u, 2).first;
    CHECK_EQ(sz(piv), 2);
  }
  // Real entries: residuals and permutation determinants.
  rep(it, 0, iters(2000)) {
    int n = (int)rnd(1, 6);
    MatD a(n, vector<double>(n));
    for (auto& row : a)
      for (double& v : row) v = rnd_real(-1e3, 1e3);
    vector<double> b(n);
    for (double& v : b) v = rnd_real(-1e3, 1e3);
    vector<double> xs(n);
    CHECK_EQ(solve_real(a, b, xs), 1);
    rep(i, 0, n) {
      double s = 0;
      rep(j, 0, n) s += a[i][j] * xs[j];
      CHECK(abs(s - b[i]) < 1e-6);
    }
  }

  int n = 300;
  MatD a(n, vector<double>(n));
  for (auto& row : a)
    for (double& v : row) v = rnd_real(-1, 1);
  vector<double> x0(n), b(n), xs(n);
  for (double& v : x0) v = rnd_real(-10, 10);
  rep(i, 0, n)
    rep(j, 0, n) b[i] += a[i][j] * x0[j];
  double t0 = elapsed();
  CHECK_EQ(solve_real(a, b, xs), 1);
  fprintf(stderr, "  300 x 300 solve: %.2fs\n", elapsed() - t0);
  rep(i, 0, n) CHECK(abs(xs[i] - x0[i]) < 1e-6);
  puts("OK linalg/gauss-real");
}
