#include "../testlib.hpp"
#include "../../content/dp/wqs.hpp"

typedef pair<ll, ll> P; // (value, count): min prefers fewer items
const ll BIG = LLONG_MAX / 4;

// Problem 1 (minimum, convex): split a (values >= 0) into exactly k
// non-empty groups, minimising the sum of squared group sums.
vector<ll> brute_split(const vector<ll>& a) {
  int n = sz(a);
  vector<ll> F(n + 1, BIG);
  rep(mask, 0, 1 << max(n - 1, 0)) { // bit i: cut after a[i]
    ll total = 0, s = 0;
    rep(i, 0, n) {
      s += a[i];
      if (i == n - 1 || mask >> i & 1) total += s * s, s = 0;
    }
    int k = __builtin_popcount(mask) + 1;
    if (n) F[k] = min(F[k], total);
  }
  if (!n) F[0] = 0;
  return F;
}
P solve_split(const vector<ll>& a, ll lam) {
  int n = sz(a);
  vector<ll> pre(n + 1);
  rep(i, 0, n) pre[i + 1] = pre[i] + a[i];
  vector<P> dp(n + 1, {BIG, 0});
  dp[0] = {0, 0};
  rep(i, 1, n + 1)
    rep(j, 0, i) {
      ll s = pre[i] - pre[j];
      dp[i] = min(dp[i], {dp[j].first + s * s + lam, dp[j].second + 1});
    }
  return dp[n];
}

// Problem 2 (maximum, concave): exactly k pairwise non-adjacent
// elements with the largest sum. Minimise -sum + lam * count.
vector<ll> brute_nonadj(const vector<ll>& a) {
  int n = sz(a);
  vector<ll> F(n + 1, -BIG);
  rep(mask, 0, 1 << n)
    if (!(mask & mask >> 1)) {
      ll s = 0;
      rep(i, 0, n)
        if (mask >> i & 1) s += a[i];
      int k = __builtin_popcount(mask);
      F[k] = max(F[k], s);
    }
  return F;
}
// T = __int128 when lam * count can pass 9e18 (the Example).
template<class T = ll> pair<T, ll> solve_nonadj(const vector<ll>& a, ll lam) {
  pair<T, ll> skip{0, 0}, take{BIG, 0};
  for (ll x : a) {
    pair<T, ll> t{skip.first - x + lam, skip.second + 1};
    skip = min(skip, take), take = t;
  }
  return min(skip, take);
}

// Problem 3 (maximum, concave): exactly k disjoint non-empty
// subarrays (adjacent allowed) with the largest total.
vector<ll> brute_subarrays(const vector<ll>& a) {
  int n = sz(a);
  vector<ll> F(n + 1, -BIG);
  int total = 1;
  rep(i, 0, n) total *= 3;
  rep(code, 0, total) { // per element: 0 out, 1 opens, 2 extends
    ll s = 0, k = 0, c = code, prev = 0;
    bool ok = true;
    rep(i, 0, n) {
      int t = c % 3;
      c /= 3;
      if (t == 2 && prev == 0) ok = false;
      if (t) s += a[i];
      if (t == 1) k++;
      prev = t;
    }
    if (ok) F[k] = max(F[k], s);
  }
  return F;
}
P solve_subarrays(const vector<ll>& a, ll lam) {
  P out{0, 0}, in{BIG, 0}; // in: the last element is chosen
  for (ll x : a) {
    P best = min(out, in);
    P open{best.first - x + lam, best.second + 1};
    P extend{in.first - x, in.second};
    out = best, in = min(open, extend);
  }
  return min(out, in);
}

ll abs_sum(const vector<ll>& a) {
  ll s = 0;
  for (ll x : a) s += abs(x);
  return s;
}

// Moderate n: exactly k non-adjacent against an O(n k) DP.
template<class T> void moderate(int n, ll maxv) {
  vector<ll> a = rnd_vec(n, 0, maxv);
  vector<vector<ll>> best(2, vector<ll>(n / 2 + 2, -BIG));
  best[0][0] = 0; // best[taken last][count]
  for (ll x : a) {
    vector<vector<ll>> nxt(2, vector<ll>(n / 2 + 2, -BIG));
    rep(c, 0, n / 2 + 1) {
      ll prev = max(best[0][c], best[1][c]);
      nxt[0][c] = prev;
      if (best[0][c] > -BIG) nxt[1][c + 1] = best[0][c] + x;
    }
    best = nxt;
  }
  auto solve = [&](ll lam) { return solve_nonadj<T>(a, lam); };
  ll C = abs_sum(a);
  vector<int> ks = {0, 1, n / 2, (n + 1) / 2};
  rep(t, 0, 20) ks.push_back((int)rnd(0, (n + 1) / 2));
  for (int k : ks) CHECK_EQ(-wqs(k, -C, C, solve), max(best[0][k], best[1][k]));
}

int main() {
  rep(it, 0, iters(1500)) {
    int n = (int)rnd(0, 10);
    vector<ll> a = rnd_vec(n, 0, it % 2 ? 1 : 5);
    vector<ll> F = brute_split(a);
    ll C = (abs_sum(a) + 1) * (abs_sum(a) + 1);
    auto solve = [&](ll lam) { return solve_split(a, lam); };
    rep(k, n ? 1 : 0, n + 1) CHECK_EQ(wqs(k, -C, C, solve), F[k]);
  }
  rep(it, 0, iters(1500)) {
    int n = (int)rnd(0, 12);
    vector<ll> a = rnd_vec(n, it % 2 ? 0 : -3, 3);
    vector<ll> F = brute_nonadj(a);
    ll C = abs_sum(a);
    auto solve = [&](ll lam) { return solve_nonadj(a, lam); };
    rep(k, 0, (n + 1) / 2 + 1) {
      CHECK_EQ(-wqs(k, -C, C, solve), F[k]);
      // At most k items: lo = 0 (the Tweak).
      ll want = *max_element(begin(F), begin(F) + k + 1);
      CHECK_EQ(-wqs(k, 0, C, solve), want);
    }
  }
  rep(it, 0, iters(400)) {
    int n = (int)rnd(0, 8);
    vector<ll> a = rnd_vec(n, -3, it % 2 ? 0 : 3);
    vector<ll> F = brute_subarrays(a);
    ll C = abs_sum(a) + 1;
    auto solve = [&](ll lam) { return solve_subarrays(a, lam); };
    rep(k, 0, n + 1) CHECK_EQ(-wqs(k, -C, C, solve), F[k]);
  }
  // Moderate n against an O(n k) DP; values up to 1e15 make
  // slopes near 1e18, so lam * count needs __int128.
  elapsed();
  moderate<ll>(3000, 1'000'000);
  moderate<__int128>(3000, 1'000'000'000'000'000);
  // Regression: a = 0, M, 0, ..., M, 0 at contest size. Slopes
  // reach m * M = 1e14 (k = m + 1 forces every 0), and lam * k
  // overflowed ll in wqs itself.
  int m = 100'000;
  ll M = 1'000'000'000;
  vector<ll> alt(2 * m + 1);
  rep(i, 0, m) alt[2 * i + 1] = M;
  ll C = abs_sum(alt);
  auto solve_alt = [&](ll lam) { return solve_nonadj<__int128>(alt, lam); };
  CHECK_EQ(-wqs(m, -C, C, solve_alt), m * M);
  CHECK_EQ(-wqs(m + 1, -C, C, solve_alt), 0LL);
  CHECK_EQ(-wqs(1, -C, C, solve_alt), M);
  // A large straight F (all values equal: every count ties at
  // lam = value).
  int n = 100'000;
  vector<ll> same(n, 7);
  auto solve_same = [&](ll lam) { return solve_nonadj(same, lam); };
  for (int k : {0, 1, 2, n / 4, n / 2})
    CHECK_EQ(-wqs(k, -7 * n, 7 * n, solve_same), 7LL * k);
  fprintf(stderr, "  large inputs: %.2fs\n", elapsed());
  puts("OK dp/wqs");
}
