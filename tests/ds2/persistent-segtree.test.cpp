#include "../testlib.hpp"
#include "../../content/ds2/persistent-segtree.hpp"

// Smallest i whose count over [0, i] exceeds k.
int brute_kth(const vector<ll>& d, ll k) {
  ll pre = 0;
  rep(i, 0, sz(d))
    if ((pre += d[i]) > k) return i;
  return -1;
}

// Random adds on random old versions; every version stays intact.
void stress_versions(int rounds, ll c) {
  rep(it, 0, iters(rounds)) {
    int n = (int)rnd(1, 12);
    PersistentSegTree st(n);
    vi root{0};
    vector<vector<ll>> arr{vector<ll>(n)};
    rep(op, 0, 40) {
      int v = (int)rnd(0, sz(root) - 1), kind = (int)rnd(0, 2);
      if (kind == 0) {
        int i = (int)rnd(0, n - 1);
        ll d = rnd(-c, c);
        root.push_back(st.add(root[v], i, d));
        arr.push_back(arr[v]), arr.back()[i] += d;
      } else if (kind == 1) {
        int l = (int)rnd(0, n), r = (int)rnd(l, n);
        CHECK_EQ(st.query(root[v], l, r),
                 accumulate(arr[v].begin() + l, arr[v].begin() + r, 0LL));
      } else { // b = version v plus some positive adds
        int b = root[v];
        vector<ll> d(n);
        rep(j, 0, (int)rnd(1, 4)) {
          int i = (int)rnd(0, n - 1);
          ll x = rnd(1, 3);
          b = st.add(b, i, x), d[i] += x;
        }
        ll k = rnd(0, accumulate(all(d), 0LL) - 1);
        CHECK_EQ(st.kth(root[v], b, k), brute_kth(d, k));
      }
    }
    rep(v, 0, sz(root))
      rep(i, 0, n) CHECK_EQ(st.query(root[v], i, i + 1), arr[v][i]);
    CHECK_EQ(st.query(0, 0, n), 0LL);
  }
}

// The Example: k-th smallest and "count < y" in a[l, r).
struct RangeKth {
  vector<ll> vals;
  PersistentSegTree st;
  vi root{0};
  RangeKth(const vector<ll>& a) : vals(a), st(0) {
    sort(all(vals)), vals.erase(unique(all(vals)), vals.end());
    st = PersistentSegTree(sz(vals));
    for (ll x : a) {
      int c = int(lower_bound(all(vals), x) - vals.begin());
      root.push_back(st.add(root.back(), c, 1));
    }
  }
  ll kth(int l, int r, int k) { return vals[st.kth(root[l], root[r], k)]; }
  ll less(int l, int r, ll y) {
    int c = int(lower_bound(all(vals), y) - vals.begin());
    return st.query(root[r], 0, c) - st.query(root[l], 0, c);
  }
};

void stress_range_kth(int rounds, ll lo, ll hi) {
  rep(it, 0, iters(rounds)) {
    int n = (int)rnd(1, 15);
    vector<ll> a = rnd_vec(n, lo, hi);
    RangeKth rk(a);
    rep(q, 0, 20) {
      int l = (int)rnd(0, n - 1), r = (int)rnd(l + 1, n);
      vector<ll> b(a.begin() + l, a.begin() + r);
      sort(all(b));
      int k = (int)rnd(0, r - l - 1);
      CHECK_EQ(rk.kth(l, r, k), b[k]);
      ll y = rnd(lo - 1, hi + 1);
      ll want = lower_bound(all(b), y) - b.begin();
      CHECK_EQ(rk.less(l, r, y), want);
    }
  }
}

int main() {
  stress_versions(2000, 5);
  stress_versions(500, 100'000'000'000'000'000); // |sums| <= 4e18
  stress_range_kth(2000, 0, 5);                  // duplicates
  stress_range_kth(1000, -1'000'000'000'000'000'000, 1'000'000'000'000'000'000);
  { // the Usage, i = 4, j = 1: kth needs v2 - v1 >= 0
    PersistentSegTree st(6);
    int v1 = st.add(0, 4, 5), v2 = st.add(v1, 1, 3);
    vector<ll> d(6);
    rep(p, 0, 6) d[p] = st.query(v2, p, p + 1) - st.query(v1, p, p + 1);
    CHECK(*min_element(all(d)) >= 0);
    rep(k, 0, 3) CHECK_EQ(st.kth(v1, v2, k), brute_kth(d, k));
    CHECK_EQ(st.query(v1, 0, 6), 5LL);
  }
  { // n = 1: a single leaf
    PersistentSegTree st(1);
    int v = st.add(st.add(0, 0, 2), 0, 3);
    CHECK_EQ(st.query(v, 0, 1), 5LL);
    CHECK_EQ(st.kth(0, v, 4), 0);
  }
  // Large: prefix versions of 1e5 values, 1e5 queries; a sample
  // of them checked against nth_element.
  elapsed();
  int n = 100'000;
  vector<ll> a = rnd_vec(n, 0, 1'000'000'000);
  RangeKth rk(a);
  CHECK(sz(rk.st.t) <= 1 + n * 18); // depth <= 17 for n <= 2^17
  ll total = 0;
  rep(q, 0, n) {
    int l = (int)rnd(0, n - 1), r = (int)rnd(l + 1, n);
    int k = (int)rnd(0, r - l - 1);
    ll got = rk.kth(l, r, k);
    total += got;
    if (q % 500 == 0) {
      vector<ll> b(a.begin() + l, a.begin() + r);
      nth_element(b.begin(), b.begin() + k, b.end());
      CHECK_EQ(got, b[k]);
    }
  }
  CHECK(total > 0);
  fprintf(stderr, "  large inputs: %.2fs\n", elapsed());
  puts("OK ds2/persistent-segtree");
}
