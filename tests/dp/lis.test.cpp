#include "../testlib.hpp"
#include "../../content/dp/lis.hpp"

// Length of the longest strictly increasing subsequence, O(n^2).
template<class T> int brute(const vector<T>& a) {
  int best = 0;
  vi len(sz(a), 1);
  rep(i, 0, sz(a)) {
    rep(j, 0, i)
      if (a[j] < a[i]) len[i] = max(len[i], len[j] + 1);
    best = max(best, len[i]);
  }
  return best;
}

// idx must be increasing positions with strictly increasing values.
template<class T> void check_valid(const vector<T>& a, const vi& idx) {
  rep(t, 0, sz(idx)) {
    CHECK(0 <= idx[t] && idx[t] < sz(a));
    if (t) CHECK(idx[t - 1] < idx[t] && a[idx[t - 1]] < a[idx[t]]);
  }
}

// Independent O(n log n) length: Fenwick prefix max over value ranks.
int fenwick_length(const vector<ll>& a) {
  vector<ll> vals = a;
  sort(all(vals));
  vals.erase(unique(all(vals)), vals.end());
  vi t(sz(vals) + 1);
  int best = 0;
  for (ll x : a) {
    int r = lower_bound(all(vals), x) - begin(vals), cur = 0;
    for (int i = r; i > 0; i &= i - 1) cur = max(cur, t[i]);
    cur++;
    for (int i = r + 1; i <= sz(vals); i += i & -i) t[i] = max(t[i], cur);
    best = max(best, cur);
  }
  return best;
}

int main() {
  rep(it, 0, iters(10000)) {
    int n = (int)rnd(0, 12);
    ll c = it % 3 == 0 ? 2 : it % 3 == 1 ? 8 : 1'000'000'000'000;
    vector<ll> a = rnd_vec(n, -c, c);
    vi idx = lis(a);
    check_valid(a, idx);
    CHECK_EQ(sz(idx), brute(a));
  }
  // Other value types: int, pairs (lexicographic), doubles.
  rep(it, 0, iters(2000)) {
    int n = (int)rnd(0, 10);
    vector<int> a(n);
    vector<pii> p(n);
    vector<double> d(n);
    rep(i, 0, n) {
      a[i] = (int)rnd(INT_MIN, INT_MAX);
      p[i] = {(int)rnd(0, 2), (int)rnd(0, 2)};
      d[i] = rnd(0, 3) / 2.0;
    }
    check_valid(a, lis(a));
    CHECK_EQ(sz(lis(a)), brute(a));
    check_valid(p, lis(p));
    CHECK_EQ(sz(lis(p)), brute(p));
    check_valid(d, lis(d));
    CHECK_EQ(sz(lis(d)), brute(d));
  }
  // Edge cases.
  CHECK(lis(vector<ll>{}).empty());
  CHECK(lis(vector<ll>{7}) == vi{0});
  CHECK(lis(vector<ll>{5, 5, 5}) == vi{2});
  CHECK(lis(vector<ll>{3, 2, 1}) == vi{2});
  CHECK((lis(vector<ll>{LLONG_MIN, LLONG_MAX}) == vi{0, 1}));
  CHECK((lis(vector<ll>{3, 1, 4, 1, 5, 9, 2, 6}) == vi{1, 2, 4, 7}));
  // Large.
  elapsed();
  int n = 100'000;
  vector<ll> a = rnd_vec(n, 0, 1'000'000);
  vi idx = lis(a);
  check_valid(a, idx);
  CHECK_EQ(sz(idx), fenwick_length(a));
  vector<ll> inc(n);
  iota(all(inc), -50'000);
  CHECK_EQ(sz(lis(inc)), n);
  fprintf(stderr, "  large inputs: %.2fs\n", elapsed());
  puts("OK dp/lis");
}
