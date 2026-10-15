#include "../testlib.hpp"
#include "../../content/linalg/matrix.hpp"

using Mat = vector<vector<ll>>;

Mat naive_mul(const Mat& a, const Mat& b) {
  int n = sz(a);
  Mat c(n, vector<ll>(n));
  rep(i, 0, n)
    rep(j, 0, n) {
      __int128 s = 0;
      rep(k, 0, n) s += (__int128)a[i][k] * b[k][j];
      c[i][j] = (ll)(s % MOD);
    }
  return c;
}

Mat naive_pow(const Mat& a, ll e) {
  int n = sz(a);
  Mat r(n, vector<ll>(n));
  rep(i, 0, n) r[i][i] = 1;
  while (e--) r = naive_mul(r, a);
  return r;
}

Matrix from(const Mat& a) {
  Matrix m(sz(a));
  m.a = a;
  return m;
}

// Walks of exactly len steps from s to t, by enumeration.
ll count_walks(const Mat& adj, int s, int t, int len) {
  if (len == 0) return s == t;
  ll total = 0;
  rep(v, 0, sz(adj))
    if (adj[s][v]) total += adj[s][v] * count_walks(adj, v, t, len - 1);
  return total;
}

int main() {
  rep(it, 0, iters(1500)) {
    int n = (int)rnd(0, 6);
    ll hi = rnd(0, 1) ? 3 : MOD - 1;
    Mat a(n, vector<ll>(n));
    for (auto& row : a)
      for (ll& v : row) v = rnd(0, 2) ? rnd(0, hi) : 0;
    ll e = rnd(0, 25);
    CHECK(from(a).pow(e).a == naive_pow(a, e));
    Mat b(n, vector<ll>(n));
    for (auto& row : b)
      for (ll& v : row) v = rnd(0, MOD - 1);
    CHECK((from(a) * from(b)).a == naive_mul(a, b));
    ll e1 = rnd(0, (ll)1e18), e2 = rnd(0, (ll)1e18);
    CHECK(from(a).pow(e1 + e2).a ==
          naive_mul(from(a).pow(e1).a, from(a).pow(e2).a));
  }
  rep(it, 0, iters(300)) { // walk counting on small graphs
    int n = (int)rnd(1, 5), len = (int)rnd(0, 7);
    Mat adj(n, vector<ll>(n));
    for (auto& row : adj)
      for (ll& v : row) v = rnd(0, 2);
    Matrix p = from(adj).pow(len);
    rep(s, 0, n)
      rep(t, 0, n) CHECK_EQ(p.a[s][t], count_walks(adj, s, t, len));
  }
  rep(it, 0, iters(300)) { // huge e: permutation, shear, diagonal
    int n = (int)rnd(1, 8);
    vi perm(n);
    iota(all(perm), 0);
    shuffle(all(perm), test_rng);
    Matrix pm(n);
    rep(i, 0, n) pm.a[i][perm[i]] = 1;
    ll e = rnd(0, (ll)1e18);
    vi q(n); // perm applied e times, e reduced by the order
    ll order = 1;
    rep(i, 0, n) {
      int len = 1;
      for (int j = perm[i]; j != i; j = perm[j]) len++;
      order = lcm(order, (ll)len);
    }
    rep(i, 0, n) {
      q[i] = i;
      rep(s, 0, e % order) q[i] = perm[q[i]];
    }
    Matrix pe = pm.pow(e);
    rep(i, 0, n)
      rep(j, 0, n) CHECK_EQ(pe.a[i][j], (ll)(q[i] == j));
    Matrix sh(2), dg(n);
    sh.a = {{1, 1}, {0, 1}};
    CHECK(sh.pow(e).a == Mat({{1, e % MOD}, {0, 1}}));
    rep(i, 0, n) dg.a[i][i] = rnd(0, MOD - 1);
    Matrix de = dg.pow(e);
    rep(i, 0, n)
      rep(j, 0, n) CHECK_EQ(de.a[i][j], i == j ? modpow(dg.a[i][i], e) : 0);
  }
  CHECK(Matrix(0).pow(5).a.empty());
  CHECK(Matrix(3).pow(0).a == Matrix(3, true).a);
  CHECK(Matrix(3).pow(1).a == Mat(3, vector<ll>(3)));

  int n = 64;
  Matrix m(n);
  for (auto& row : m.a)
    for (ll& v : row) v = rnd(0, MOD - 1);
  double t0 = elapsed();
  Matrix big = m.pow((ll)1e18);
  fprintf(stderr, "  64 x 64 to the power 1e18: %.2fs\n", elapsed() - t0);
  CHECK((big * m).a == (m * big).a);
  puts("OK linalg/matrix");
}
