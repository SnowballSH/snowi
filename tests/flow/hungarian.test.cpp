#include "../testlib.hpp"
#include "../../content/flow/hungarian.hpp"

using Mat = vector<vector<ll>>;

// Minimum over all injections rows -> columns.
ll brute(const Mat& a, int n, int m) {
  vi cols(m);
  iota(all(cols), 0);
  ll best = LLONG_MAX;
  do {
    ll c = 0;
    rep(i, 0, n) c += a[i][cols[i]];
    best = min(best, c);
  } while (next_permutation(all(cols)));
  return n ? best : 0;
}

void check_match(const Mat& a, int n, int m, ll cost, const vi& match) {
  CHECK_EQ(sz(match), n);
  vector<bool> used(m);
  ll c = 0;
  rep(i, 0, n) {
    CHECK(0 <= match[i] && match[i] < m && !used[match[i]]);
    used[match[i]] = true;
    c += a[i][match[i]];
  }
  CHECK_EQ(c, cost);
}

// Random duals u, v (v <= 0, and 0 on free columns) and a planted
// injection whose entries are tight: it is the unique optimum.
pair<Mat, vi> planted(int n, int m, ll& opt) {
  vi perm(m);
  iota(all(perm), 0);
  shuffle(all(perm), test_rng);
  perm.resize(n);
  vector<ll> u(n), v(m, 0);
  opt = 0;
  rep(i, 0, n) u[i] = rnd(-1e6, 1e6), opt += u[i];
  for (int j : perm) v[j] = rnd(-1e6, 0), opt += v[j];
  Mat a(n, vector<ll>(m));
  rep(i, 0, n)
    rep(j, 0, m) a[i][j] = u[i] + v[j] + rnd(1, 1e6);
  rep(i, 0, n) a[i][perm[i]] = u[i] + v[perm[i]];
  return {a, perm};
}

int main() {
  rep(it, 0, iters(3000)) {
    int m = (int)rnd(0, 6), n = (int)rnd(0, m);
    int kind = (int)rnd(0, 2);
    ll A = kind == 0 ? 5 : kind == 1 ? (ll)1e9 : LLONG_MAX / (n + 4);
    Mat a(n, vector<ll>(m));
    for (auto& row : a)
      for (ll& x : row)
        x = kind == 2 && rnd(0, 1) ? A * (rnd(0, 1) ? 1 : -1) : rnd(-A, A);
    auto [cost, match] = hungarian(a);
    CHECK_EQ(cost, brute(a, n, m));
    check_match(a, n, m, cost, match);
  }
  rep(it, 0, iters(200)) {
    int m = (int)rnd(1, 40), n = (int)rnd(1, m);
    ll opt;
    auto [a, perm] = planted(n, m, opt);
    auto [cost, match] = hungarian(a);
    CHECK_EQ(cost, opt);
    CHECK(match == perm);
  }
  elapsed();
  {
    ll opt;
    auto [a, perm] = planted(200, 300, opt);
    auto [cost, match] = hungarian(a);
    CHECK_EQ(cost, opt);
    CHECK(match == perm);
  }
  { // a[i][j] = -i j: identity is optimal (rearrangement).
    int n = 300;
    Mat a(n, vector<ll>(n));
    ll opt = 0;
    rep(i, 0, n)
      rep(j, 0, n) a[i][j] = -(ll)i * j;
    rep(i, 0, n) opt -= (ll)i * i;
    auto [cost, match] = hungarian(a);
    CHECK_EQ(cost, opt);
    check_match(a, n, n, cost, match);
  }
  fprintf(stderr, "  large cases: %.2fs\n", elapsed());
  puts("OK flow/hungarian");
}
