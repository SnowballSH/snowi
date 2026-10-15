#include "../testlib.hpp"
#include "../../content/ds2/wavelet-matrix.hpp"

int brute_less(const vector<ll>& a, int l, int r, ll x) {
  int c = 0;
  rep(i, l, r) c += a[i] < x;
  return c;
}

// Values of a[l, r) in [lo, hi); none when lo >= hi.
int brute_count(const vector<ll>& a, int l, int r, ll lo, ll hi) {
  int c = 0;
  rep(i, l, r) c += lo <= a[i] && a[i] < hi;
  return c;
}

// The Tweak "sum of the k smallest": prefix sums per level h of
// the values in their order after level h's partition.
struct SumSmallest {
  WaveletMatrix wm;
  vector<vector<ll>> ps;
  SumSmallest(vector<ll> a) : wm(a), ps(wm.lg) {
    for (int h = wm.lg; h--;) {
      stable_partition(all(a), [&](ll x) { return ~x >> h & 1; });
      ps[h].assign(sz(a) + 1, 0);
      rep(i, 0, sz(a)) ps[h][i + 1] = ps[h][i] + a[i];
    }
  }
  ll sum(int l, int r, int k) { // 0 <= k <= r - l
    ll res = 0, s = 0;
    for (int h = wm.lg; h--;) {
      int l0 = wm.rank0(h, l), r0 = wm.rank0(h, r);
      if (k < r0 - l0) l = l0, r = r0;
      else {
        k -= r0 - l0, res |= 1LL << h, s += ps[h][r0] - ps[h][l0];
        l += wm.zeros[h] - l0, r += wm.zeros[h] - r0;
      }
    }
    return s + k * res; // the k copies of the k-th value left
  }
};

void stress_sum(int rounds, int maxn, ll hi) {
  rep(it, 0, iters(rounds)) {
    int n = (int)rnd(0, maxn);
    vector<ll> a = rnd_vec(n, 0, hi);
    SumSmallest ss(a);
    rep(q, 0, 30) {
      int l = (int)rnd(0, n), r = (int)rnd(l, n), k = (int)rnd(0, r - l);
      vector<ll> b(a.begin() + l, a.begin() + r);
      sort(all(b));
      CHECK_EQ(ss.sum(l, r, k), accumulate(b.begin(), b.begin() + k, 0LL));
    }
  }
}

void stress(int rounds, int maxn, ll hi) {
  rep(it, 0, iters(rounds)) {
    int n = (int)rnd(0, maxn);
    vector<ll> a = rnd_vec(n, 0, hi);
    WaveletMatrix wm(a);
    rep(q, 0, 30) {
      int l = (int)rnd(0, n), r = (int)rnd(l, n);
      if (l < r) {
        vector<ll> b(a.begin() + l, a.begin() + r);
        sort(all(b));
        int k = (int)rnd(0, r - l - 1);
        CHECK_EQ(wm.kth(l, r, k), b[k]);
      }
      ll x = rnd(-2, hi + 2), y = rnd(-2, hi + 2);
      if (rnd(0, 3) == 0) x = a.empty() ? 0 : a[rnd(0, n - 1)];
      CHECK_EQ(wm.count_less(l, r, x), brute_less(a, l, r, x));
      CHECK_EQ(wm.count(l, r, x, y), brute_count(a, l, r, x, y));
      CHECK_EQ(wm.rank(x, r), (int)count(a.begin(), a.begin() + r, x));
    }
  }
}

int main() {
  { // the Example
    WaveletMatrix wm({5, 1, 4, 1, 3, 9});
    CHECK_EQ(wm.kth(1, 5, 2), 3LL);
    CHECK_EQ(wm.count_less(1, 5, 4), 3);
    CHECK_EQ(wm.count(0, 6, 3, 6), 3);
    CHECK_EQ(wm.rank(1, 4), 2);
    int c = wm.count_less(0, 6, 5); // largest value < 5
    CHECK_EQ(wm.kth(0, 6, c - 1), 4LL);
    // regressions: x + 1 overflowed in rank, lo > hi went negative
    CHECK_EQ(wm.rank(LLONG_MAX, 6), 0);
    CHECK_EQ(wm.count(0, 6, 5, 3), 0);
  }
  { // the Tweak on duplicates: a = [1, 1], sum of the 1 smallest
    SumSmallest ss({1, 1});
    CHECK_EQ(ss.sum(0, 2, 1), 1LL);
  }
  stress(3000, 10, 0);                // all zeros: no levels
  stress(3000, 20, 3);                // many duplicates
  stress(2000, 70, 1000);             // more than one 64-bit word
  stress(1000, 200, (1LL << 62) - 1); // all 62 levels
  stress_sum(2000, 20, 3);
  stress_sum(1000, 100, 1'000'000'000);
  // extremes for x and k
  {
    vector<ll> a{(1LL << 62) - 1, 0, 1LL << 61};
    WaveletMatrix wm(a);
    CHECK_EQ(wm.kth(0, 3, 2), (1LL << 62) - 1);
    CHECK_EQ(wm.count_less(0, 3, LLONG_MAX), 3);
    CHECK_EQ(wm.count_less(0, 3, LLONG_MIN), 0);
    CHECK_EQ(wm.count_less(0, 3, 1LL << 61), 1);
    CHECK_EQ(wm.rank(0, 3), 1);
    CHECK_EQ(wm.rank(LLONG_MAX, 3), 0);
    CHECK_EQ(wm.rank(LLONG_MIN, 3), 0);
    CHECK_EQ(wm.count(0, 3, LLONG_MAX, LLONG_MIN), 0);
  }
  // Large: 1e5 values, 3e5 queries; a sample checked by brute force.
  elapsed();
  int n = 100'000;
  vector<ll> a = rnd_vec(n, 0, 1'000'000'000);
  WaveletMatrix wm(a);
  ll total = 0;
  rep(q, 0, n) {
    int l = (int)rnd(0, n - 1), r = (int)rnd(l + 1, n);
    int k = (int)rnd(0, r - l - 1);
    ll x = rnd(0, 1'000'000'000);
    ll v = wm.kth(l, r, k);
    int c = wm.count_less(l, r, x), d = wm.rank(x, r);
    total += v + c + d;
    if (q % 1000 == 0) {
      vector<ll> b(a.begin() + l, a.begin() + r);
      nth_element(b.begin(), b.begin() + k, b.end());
      CHECK_EQ(v, b[k]);
      CHECK_EQ(c, brute_less(a, l, r, x));
    }
  }
  CHECK(total > 0);
  fprintf(stderr, "  large inputs: %.2fs\n", elapsed());
  puts("OK ds2/wavelet-matrix");
}
