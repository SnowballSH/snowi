#include "../testlib.hpp"
#include "../../content/linalg/gauss-mod.hpp"

using Mat = vector<vector<ll>>;

ll perm_det(const Mat& a) {
  int n = sz(a);
  vi p(n);
  iota(all(p), 0);
  ll res = 0;
  do {
    ll term = 1;
    int inversions = 0;
    rep(i, 0, n) {
      term = term * a[i][p[i]] % MOD;
      rep(j, i + 1, n) inversions += p[i] > p[j];
    }
    res = (res + (inversions % 2 ? MOD - term : term)) % MOD;
  } while (next_permutation(all(p)));
  return res;
}

// Largest k with a nonzero k x k minor.
int minor_rank(const Mat& a, int m) {
  int n = sz(a), best = 0;
  rep(rows, 0, 1 << n)
    rep(cols, 0, 1 << m) {
      int k = __builtin_popcount(rows);
      if (k != __builtin_popcount(cols) || k <= best) continue;
      Mat sub;
      rep(i, 0, n)
        if (rows >> i & 1) {
          sub.emplace_back();
          rep(j, 0, m)
            if (cols >> j & 1) sub.back().push_back(a[i][j]);
        }
      if (perm_det(sub)) best = k;
    }
  return best;
}

Mat random_matrix(int n, int m) {
  Mat a(n, vector<ll>(m));
  int kind = (int)rnd(0, 3);
  if (kind == 3) { // rank <= k: product of n x k and k x m
    int k = (int)rnd(0, min(n, m));
    Mat u(n, vector<ll>(k)), v(k, vector<ll>(m));
    for (auto& row : u)
      for (ll& x : row) x = rnd(0, MOD - 1);
    for (auto& row : v)
      for (ll& x : row) x = rnd(0, MOD - 1);
    rep(i, 0, n)
      rep(j, 0, m)
        rep(t, 0, k) a[i][j] = (a[i][j] + u[i][t] * v[t][j]) % MOD;
    return a;
  }
  for (auto& row : a)
    for (ll& x : row)
      x = kind == 0   ? rnd(0, 2)
          : kind == 1 ? rnd(0, MOD - 1)
                      : (rnd(0, 1) ? MOD - 1 : rnd(0, 1));
  if (n > 1 && rnd(0, 3) == 0) a[rnd(0, n - 1)] = a[rnd(0, n - 1)];
  return a;
}

vector<ll> mul(const Mat& a, const vector<ll>& x) {
  vector<ll> y(sz(a));
  rep(i, 0, sz(a))
    rep(j, 0, sz(x)) y[i] = (y[i] + a[i][j] * x[j]) % MOD;
  return y;
}

void check_rref(const Mat& a, int m, const vi& piv) {
  int r = sz(piv);
  rep(i, 0, r) {
    CHECK(i == 0 || piv[i - 1] < piv[i]);
    rep(j, 0, piv[i]) CHECK_EQ(a[i][j], 0LL);
    rep(k, 0, sz(a)) CHECK_EQ(a[k][piv[i]], (ll)(k == i));
  }
  rep(i, r, sz(a))
    rep(j, 0, m) CHECK_EQ(a[i][j], 0LL);
}

int main() {
  rep(it, 0, iters(3000)) {
    int n = (int)rnd(0, 4), m = (int)rnd(0, 4);
    const Mat a = random_matrix(n, m); // only gauss may change it
    Mat g = a;
    auto [piv, d] = gauss(g, m);
    CHECK_EQ(sz(piv), minor_rank(a, m));
    check_rref(g, m, piv);
    if (n == m) {
      CHECK_EQ(d, perm_det(a));
      CHECK_EQ(det(a), d);
      Mat inv = mat_inv(a);
      CHECK_EQ(sz(inv), d ? n : 0);
      if (d)
        rep(j, 0, n) {
          vector<ll> col(n), e(n);
          rep(i, 0, n) col[i] = inv[i][j];
          e[j] = 1;
          CHECK(mul(a, col) == e);
        }
    }
    vector<ll> b(n), x0(m);
    for (ll& v : x0) v = rnd(0, MOD - 1);
    if (rnd(0, 1)) b = mul(a, x0);
    else
      for (ll& v : b) v = rnd(0, 2);
    Mat ab = a;
    rep(i, 0, n) ab[i].push_back(b[i]);
    int ra = minor_rank(a, m), rab = minor_rank(ab, m + 1);
    vector<ll> x(m, 7);
    int r = solve_linear(a, b, x);
    CHECK_EQ(r, ra == rab ? ra : -1);
    if (r >= 0) {
      CHECK_EQ(sz(x), m);
      CHECK(mul(a, x) == b);
    }
  }
  vector<ll> x(3, 7); // no equations: x = 0 is a solution
  CHECK_EQ(solve_linear({}, {}, x), 0);
  CHECK(x == vector<ll>(3, 0));
  CHECK(mat_inv({}).empty());
  CHECK_EQ(det({}), 1LL);
  const Mat s = {{1, 1}, {1, 1}}; // singular: {}, s untouched
  CHECK(mat_inv(s).empty());
  { // The Usage lines in order: only the last one, gauss, is in
    // place, so the others all see the original matrix.
    Mat u = {{2, 1}, {1, 3}};
    vector<ll> b = {3, 4}, ux(2);
    CHECK_EQ(solve_linear(u, b, ux), 2);
    CHECK(mul({{2, 1}, {1, 3}}, ux) == b);
    CHECK_EQ(det(u), 5LL);
    Mat inv = mat_inv(u);
    CHECK(mul(u, {inv[0][0], inv[1][0]}) == vector<ll>({1, 0}));
    CHECK(mul(u, {inv[0][1], inv[1][1]}) == vector<ll>({0, 1}));
    vi piv = gauss(u, 2).first;
    CHECK(sz(piv) == 2 && u == Mat({{1, 0}, {0, 1}}));
  }

  int n = 300;
  Mat a(n, vector<ll>(n));
  for (auto& row : a)
    for (ll& v : row) v = rnd(0, MOD - 1);
  vector<ll> x0(n), xs(n);
  for (ll& v : x0) v = rnd(0, MOD - 1);
  double t0 = elapsed();
  CHECK_EQ(solve_linear(a, mul(a, x0), xs), n);
  CHECK(xs == x0);
  Mat inv = mat_inv(a);
  CHECK_EQ(sz(inv), n);
  rep(j, 0, 5) {
    vector<ll> col(n), ej(n);
    rep(i, 0, n) col[i] = inv[i][j];
    ej[j] = 1;
    CHECK(mul(a, col) == ej);
  }
  fprintf(stderr, "  300 x 300 solve + inverse: %.2fs\n", elapsed() - t0);
  puts("OK linalg/gauss-mod");
}
