#include "../testlib.hpp"
#include "../../content/dp/knuth.hpp"

// O(n^3) interval DP, dp[i][i + 1] = 0.
template<class W> vector<vector<ll>> brute(int n, W w) {
  vector<vector<ll>> dp(n + 1, vector<ll>(n + 1));
  rep(len, 2, n + 1)
    rep(l, 0, n - len + 1) {
      int r = l + len;
      ll best = LLONG_MAX;
      rep(k, l + 1, r) best = min(best, dp[l][k] + dp[k][r]);
      dp[l][r] = best + w(l, r);
    }
  return dp;
}

// Random w with the quadrangle inequality and monotone under
// inclusion: f(sum) with f convex increasing, plus pair weights.
struct RandomW {
  vector<ll> pre;
  vector<vector<ll>> pairs; // pairs[l][r]: weight of pairs in [l, r)
  ll lin, quad;
  RandomW(int n, ll c) : pre(n + 1), pairs(n + 1, vector<ll>(n + 1)) {
    lin = rnd(0, 3), quad = rnd(0, 2);
    vector<vector<ll>> p(n, vector<ll>(n));
    rep(x, 0, n)
      rep(y, x + 1, n) p[x][y] = rnd(0, 2) ? 0 : rnd(0, c);
    rep(l, 0, n + 1)
      rep(r, l + 1, n + 1) {
        pairs[l][r] = pairs[l][r - 1];
        rep(x, l, r - 1) pairs[l][r] += p[x][r - 1];
      }
    rep(i, 0, n) pre[i + 1] = pre[i] + rnd(0, c);
  }
  ll operator()(int l, int r) const {
    ll s = pre[r] - pre[l];
    return lin * s + quad * s * s + pairs[l][r];
  }
};

void stress(int rounds, int maxn, ll c) {
  rep(it, 0, iters(rounds)) {
    int n = (int)rnd(0, maxn);
    RandomW w(n, c);
    CHECK(knuth(n, w) == brute(n, w));
  }
}

// For large n: f(sum) plus pairs x < y inside with weight
// p[x] * q[y], in O(1) from prefix sums.
struct ProductW {
  vector<ll> pre, P, Q, QP;
  ProductW(int n) : pre(n + 1), P(n + 1), Q(n + 1), QP(n + 1) {
    rep(i, 0, n) {
      ll p = rnd(0, 100), q = rnd(0, 100);
      pre[i + 1] = pre[i] + rnd(0, 1000);
      P[i + 1] = P[i] + p, Q[i + 1] = Q[i] + q;
      QP[i + 1] = QP[i] + q * P[i];
    }
  }
  ll operator()(int l, int r) const {
    ll s = pre[r] - pre[l];
    return s * s + (QP[r] - QP[l]) - P[l] * (Q[r] - Q[l]);
  }
};

// Cost of merging len unit piles: a balanced tree.
ll units(int len) {
  int q = __lg(len);
  return (ll)len * q + 2 * (len - (1 << q));
}

int main() {
  stress(5000, 7, 3); // tiny values: many ties
  stress(2000, 14, 20);
  stress(200, 40, 1'000'000);
  CHECK(knuth(0, [](int, int) { return 1LL; }) == vector<vector<ll>>{{0}});
  ll one = knuth(1, [](int, int) { return 1LL; })[0][1];
  ll two = knuth(2, [](int, int) { return 7LL; })[0][2];
  CHECK_EQ(one, 0LL);
  CHECK_EQ(two, 7LL);
  // Optimal BST via n + 1 boundaries (the Example), against a
  // direct O(n^3) DP over roots.
  rep(it, 0, iters(500)) {
    int n = (int)rnd(0, 9);
    vector<ll> f = rnd_vec(n, 0, 10), pre(n + 1);
    rep(i, 0, n) pre[i + 1] = pre[i] + f[i];
    vector<vector<ll>> e(n + 1, vector<ll>(n + 1)); // keys [l, r)
    rep(len, 1, n + 1)
      rep(l, 0, n - len + 1) {
        int r = l + len;
        e[l][r] = LLONG_MAX;
        rep(k, l, r) e[l][r] = min(e[l][r], e[l][k] + e[k + 1][r]);
        e[l][r] += pre[r] - pre[l];
      }
    auto dp = knuth(n + 1, [&](int l, int r) { return pre[r - 1] - pre[l]; });
    CHECK_EQ(dp[0][n + 1], e[0][n]);
  }
  // Large: unit piles, every interval known in closed form.
  elapsed();
  int n = 2000;
  auto dp = knuth(n, [](int l, int r) { return (ll)(r - l); });
  rep(l, 0, n)
    rep(r, l + 1, n + 1) CHECK_EQ(dp[l][r], units(r - l));
  ProductW w(n);
  dp = knuth(n, w);
  vector<vector<ll>> small = brute(150, w);
  rep(l, 0, 151)
    rep(r, l + 1, 151) CHECK_EQ(dp[l][r], small[l][r]);
  fprintf(stderr, "  large inputs: %.2fs\n", elapsed());
  puts("OK dp/knuth");
}
