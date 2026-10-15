#include "../testlib.hpp"
#include "../../content/dp/dc-opt.hpp"

const ll BIG = LLONG_MAX / 4;

// new[i] = min over j < i of old[j] + cost(j, i), capped at BIG.
template<class F> vector<ll> brute(const vector<ll>& old, F cost) {
  vector<ll> res(sz(old), BIG);
  rep(i, 1, sz(old))
    rep(j, 0, i) res[i] = min(res[i], old[j] + cost(j, i));
  return res;
}

// A random cost with the quadrangle inequality: f(S_i - S_j) with
// f convex, plus weights of pairs inside [j, i), plus terms in i
// alone and in j alone (which keep the inequality, may be < 0).
struct RandomCost {
  vector<ll> pre, gi, hj;
  vector<vector<ll>> pairs; // pairs[j][i]: weight of pairs in [j, i)
  ll quad;
  RandomCost(int n, ll c) : pre(n + 1), gi(n + 1), hj(n + 1) {
    quad = rnd(0, 2);
    vector<vector<ll>> w(n, vector<ll>(n));
    rep(x, 0, n)
      rep(y, x + 1, n) w[x][y] = rnd(0, 1) * rnd(0, c);
    pairs.assign(n + 1, vector<ll>(n + 1));
    rep(j, 0, n + 1)
      rep(i, j + 1, n + 1) {
        pairs[j][i] = pairs[j][i - 1];
        rep(x, j, i - 1) pairs[j][i] += w[x][i - 1];
      }
    rep(i, 0, n) pre[i + 1] = pre[i] + rnd(0, c);
    rep(i, 0, n + 1) gi[i] = rnd(-c, c), hj[i] = rnd(-c, c);
  }
  ll operator()(int j, int i) const {
    ll s = pre[i] - pre[j];
    return quad * s * s + pairs[j][i] + gi[i] + hj[j];
  }
};

void stress(int rounds, int maxn, ll c) {
  rep(it, 0, iters(rounds)) {
    int n = (int)rnd(0, maxn);
    RandomCost cost(n, c);
    vector<ll> dp(n + 1, BIG);
    if (rnd(0, 3)) dp[0] = 0;
    else
      rep(i, 0, n + 1)
        if (rnd(0, 1)) dp[i] = rnd(-c, c);
    int layers = (int)rnd(1, n + 2);
    rep(k, 0, layers) {
      vector<ll> want = brute(dp, cost);
      dp = dc_layer(dp, cost);
      CHECK(dp == want);
    }
  }
}

// Regression for the Precondition: no quadrangle inequality,
// only the smallest optimal j non-decreasing in i. Larger j often
// tie with it, so taking the largest optimal j (`<=`) must fail.
void stress_smallest_opt(int rounds, int maxn) {
  rep(it, 0, iters(rounds)) {
    int n = (int)rnd(1, maxn), o = 0;
    vector<ll> old = rnd_vec(n + 1, -5, 5);
    vector<vector<ll>> w(n + 1, vector<ll>(n + 1));
    rep(i, 1, n + 1) {
      o = (int)rnd(o, i - 1); // the smallest optimal j for i
      ll best = rnd(-5, 5);
      rep(j, 0, i) {
        ll v = best + (j < o ? rnd(1, 3) : j > o ? rnd(0, 2) : 0);
        w[j][i] = v - old[j];
      }
    }
    auto cost = [&](int j, int i) { return w[j][i]; };
    CHECK(dc_layer(old, cost) == brute(old, cost));
  }
}

int main() {
  stress(4000, 6, 3);   // tiny values: many ties
  stress(3000, 14, 10); // longer arrays
  stress(500, 30, 1'000'000);
  stress_smallest_opt(3000, 12);
  // One element: only new[0], which has no j.
  CHECK(dc_layer(vector<ll>{5}, [](int, int) { return 0LL; }) ==
        vector<ll>{BIG});
  // Unreachable states stay capped over many layers.
  vector<ll> dp(4, BIG);
  dp[0] = 0;
  auto huge = [](int j, int i) { return (ll)(i - j) << 60; };
  rep(k, 0, 5) dp = dc_layer(dp, huge);
  CHECK(dp == vector<ll>(4, BIG));
  // Large: squared group sums. Check sampled i against a direct
  // scan and that cost is called O(n log n) times.
  elapsed();
  int n = 100'000;
  vector<ll> pre(n + 1), old(n + 1);
  rep(i, 0, n) pre[i + 1] = pre[i] + rnd(0, 1000);
  rep(i, 0, n + 1) old[i] = rnd(0, 1'000'000'000'000);
  ll calls = 0;
  auto sq = [&](int j, int i) {
    calls++;
    return (pre[i] - pre[j]) * (pre[i] - pre[j]);
  };
  vector<ll> res = dc_layer(old, sq);
  ll dc_calls = calls;
  CHECK(dc_calls <= (ll)(n + 1) * 2 * 18);
  rep(t, 0, 100) {
    int i = (int)rnd(1, n);
    ll want = BIG;
    rep(j, 0, i) want = min(want, old[j] + sq(j, i));
    CHECK_EQ(res[i], want);
  }
  fprintf(stderr, "  large inputs: %.2fs, %lld cost calls\n", elapsed(),
          dc_calls);
  puts("OK dp/dc-opt");
}
