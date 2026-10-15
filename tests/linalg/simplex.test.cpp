#include "../testlib.hpp"
#include "../../content/linalg/simplex.hpp"

using Mat = vector<vector<ll>>;
using i128 = __int128;
const double INF = numeric_limits<double>::infinity();

ll det(const Mat& a) {
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

// Exact optimum over the vertices of {A x <= b, 0 <= x <= box}:
// -INF if empty. Vertices solve n tight constraints (Cramer).
double vertex_opt(const Mat& A, const vector<ll>& b, const vector<ll>& c,
                  ll box) {
  int n = sz(c);
  Mat rows = A;
  vector<ll> rhs = b;
  rep(j, 0, n) {
    vector<ll> lo(n), hi(n);
    lo[j] = -1, hi[j] = 1;
    rows.push_back(lo), rhs.push_back(0);
    rows.push_back(hi), rhs.push_back(box);
  }
  int k = sz(rows);
  bool found = false;
  i128 bp = 0, bq = 1; // best objective bp / bq, bq > 0
  rep(mask, 0, 1 << k) {
    if (__builtin_popcount(mask) != n) continue;
    Mat s;
    vector<ll> sb;
    rep(i, 0, k)
      if (mask >> i & 1) s.push_back(rows[i]), sb.push_back(rhs[i]);
    ll d = det(s);
    if (!d) continue;
    vector<ll> num(n); // x_j = num[j] / d
    rep(j, 0, n) {
      Mat t = s;
      rep(i, 0, n) t[i][j] = sb[i];
      num[j] = det(t);
    }
    bool ok = true;
    rep(i, 0, k) {
      i128 lhs = 0;
      rep(j, 0, n) lhs += (i128)rows[i][j] * num[j];
      if (d > 0 ? lhs > (i128)rhs[i] * d : lhs < (i128)rhs[i] * d) ok = false;
    }
    if (!ok) continue;
    i128 p = 0, q = d;
    rep(j, 0, n) p += (i128)c[j] * num[j];
    if (q < 0) p = -p, q = -q;
    if (!found || p * bq > bp * q) bp = p, bq = q, found = true;
  }
  return found ? (double)bp / (double)bq : -INF;
}

// Feasible, optimal and dual-certified (strong duality). Errors
// are relative to the size of the numbers compared.
void check_solution(LPSolver& lp, const vector<vector<double>>& A,
                    const vector<double>& b, const vector<double>& c,
                    const vector<double>& x, double val, double tol) {
  int m = sz(b), n = sz(c);
  auto close = [&](double u, double v) {
    return abs(u - v) <= tol * (1 + abs(u) + abs(v));
  };
  CHECK_EQ(sz(x), n);
  double cx = 0;
  rep(j, 0, n) {
    CHECK(x[j] >= -tol);
    cx += c[j] * x[j];
  }
  CHECK(close(cx, val));
  rep(i, 0, m) {
    double s = 0;
    rep(j, 0, n) s += A[i][j] * x[j];
    CHECK(s <= b[i] || close(s, b[i]));
  }
  vector<double> y(m);
  rep(j, 0, n + 1)
    if (lp.N[j] >= n) y[lp.N[j] - n] = lp.D[m][j];
  double by = 0;
  rep(i, 0, m) {
    CHECK(y[i] >= -tol);
    by += b[i] * y[i];
  }
  CHECK(close(by, val));
  rep(j, 0, n) {
    double s = 0;
    rep(i, 0, m) s += A[i][j] * y[i];
    CHECK(s >= c[j] || close(s, c[j]));
  }
}

double run(const vector<vector<double>>& A, const vector<double>& b,
           const vector<double>& c, double tol) {
  LPSolver lp(A, b, c);
  vector<double> x;
  double val = lp.solve(x);
  if (abs(val) != INF) check_solution(lp, A, b, c, x, val, tol);
  return val;
}

int main() {
  rep(it, 0, iters(3000)) {
    int n = (int)rnd(1, 3), m = (int)rnd(0, 4);
    Mat A(m, vector<ll>(n));
    vector<ll> b(m), c(n);
    for (auto& row : A)
      for (ll& v : row) v = rnd(0, 2) ? rnd(-4, 4) : 0;
    for (ll& v : b) v = rnd(0, 2) ? rnd(-6, 6) : 0;
    for (ll& v : c) v = rnd(-4, 4);
    if (m > 1 && rnd(0, 3) == 0) A[1] = A[0], b[1] = b[0];
    double small = vertex_opt(A, b, c, 1000);
    double big = vertex_opt(A, b, c, 2000);
    double want = small == -INF ? -INF : big > small + 1e-9 ? INF : small;
    vector<vector<double>> Ad;
    for (auto& row : A) Ad.emplace_back(all(row));
    double got = run(Ad, {all(b)}, {all(c)}, 1e-6);
    if (abs(want) == INF) CHECK_EQ(got, want);
    else CHECK(abs(got - want) < 1e-6);
  }
  // Beale's example: Dantzig's rule alone can cycle on it.
  vector<vector<double>> A = {
      {0.25, -60, -0.04, 9}, {0.5, -90, -0.02, 3}, {0, 0, 1, 0}};
  CHECK(abs(run(A, {0, 0, 1}, {0.75, -150, 0.02, -6}, 1e-9) - 0.05) < 1e-9);
  CHECK_EQ(run({}, {}, {1}, 1e-9), INF); // no constraints
  CHECK_EQ(run({}, {}, {-1, 0}, 1e-9), 0.0);
  CHECK_EQ(run({{1}, {-1}}, {1, -2}, {1}, 1e-9), -INF);
  CHECK_EQ(run({{0}}, {-1}, {1}, 1e-9), -INF); // 0 <= -1
  CHECK(abs(run({{1, 1}, {-1, -1}}, {2, -2}, {1, 2}, 1e-9) - 4) < 1e-9);
  // Klee-Minty cube: Dantzig's rule visits 2^n - 1 vertices.
  int km = 12;
  vector<vector<double>> K(km, vector<double>(km));
  vector<double> kb(km), kc(km);
  rep(i, 0, km) {
    rep(j, 0, i) K[i][j] = pow(2, i - j + 1);
    K[i][i] = 1, kb[i] = pow(5, i + 1), kc[i] = pow(2, km - 1 - i);
  }
  CHECK(abs(run(K, kb, kc, 1e-9) / pow(5, km) - 1) < 1e-9);

  double t0 = elapsed();
  rep(round, 0, 2) { // random 100 x 100, some b_i < 0 (phase 1)
    int n = 100, m = 100;
    vector<vector<double>> R(m, vector<double>(n));
    vector<double> x0(n), rb(m), rc(n);
    for (auto& row : R)
      for (double& v : row) v = rnd_real(-1, 1);
    for (double& v : x0) v = rnd_real(0, 1);
    rep(i, 0, m) {
      rep(j, 0, n) rb[i] += R[i][j] * x0[j];
      rb[i] += round ? rnd_real(0, 1) : 0.5;
    }
    for (double& v : rc) v = rnd_real(-1, 1);
    rep(j, 0, n) { // keep it bounded: x_j <= 2
      R.emplace_back(n), R.back()[j] = 1;
      rb.push_back(2);
    }
    double val = run(R, rb, rc, 1e-6);
    CHECK(abs(val) != INF);
  }
  fprintf(stderr, "  two random 200 x 100 LPs: %.2fs\n", elapsed() - t0);
  puts("OK linalg/simplex");
}
