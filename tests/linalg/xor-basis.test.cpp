#include "../testlib.hpp"
#include "../../content/linalg/xor-basis.hpp"

// All subset XORs of a, sorted.
vector<ll> span_of(const vector<ll>& a) {
  set<ll> s{0};
  for (ll x : a) {
    set<ll> t = s;
    for (ll v : s) t.insert(v ^ x);
    s = t;
  }
  return {all(s)};
}

ll random_value(int width, const vector<ll>& prev) {
  if (!prev.empty() && rnd(0, 3) == 0) { // dependent on earlier ones
    ll x = 0;
    for (ll v : prev)
      if (rnd(0, 1)) x ^= v;
    return x;
  }
  return rnd(0, (1LL << width) - 1);
}

int main() {
  rep(it, 0, iters(1500)) {
    int width = rnd(0, 4) ? (int)rnd(1, 6) : 60;
    int n = (int)rnd(0, 10);
    XorBasis xb;
    vector<ll> a;
    rep(i, 0, n) {
      ll x = random_value(width, a);
      vector<ll> before = span_of(a);
      a.push_back(x);
      vector<ll> s = span_of(a);
      CHECK_EQ(xb.insert(x), sz(s) > sz(before));
      CHECK_EQ(1LL << xb.rank, (ll)sz(s));
      rep(k, 0, sz(s)) CHECK_EQ(xb.kth(k), s[k]);
      rep(q, 0, 5) {
        ll y = q == 0 ? 0 : random_value(width, a);
        ll mx = 0, mn = LLONG_MAX;
        for (ll v : s) mx = max(mx, y ^ v), mn = min(mn, y ^ v);
        CHECK_EQ(xb.max_xor(y), mx);
        CHECK_EQ(xb.min_xor(y), mn);
        CHECK_EQ(xb.contains(y), binary_search(all(s), y));
      }
    }
  }
  XorBasis full; // all-ones minus one bit: independent (60 even)
  const ll ones = (1LL << 60) - 1;
  rep(b, 0, 60) CHECK(full.insert(ones ^ 1LL << b));
  CHECK(full.rank == 60 && !full.insert(1LL << 59));
  CHECK_EQ(full.max_xor(), ones);
  CHECK_EQ(full.min_xor(ones), 0LL);
  for (ll k : {0LL, 1LL, 12345LL, ones}) CHECK_EQ(full.kth(k), k);

  rep(it, 0, iters(1500)) {
    int width = rnd(0, 4) ? (int)rnd(1, 5) : 60;
    int n = (int)rnd(1, width == 60 ? 8 : 12);
    vector<ll> a;
    PrefixXorBasis pb;
    rep(r, 0, n) {
      a.push_back(random_value(width, a));
      pb.insert(a[r], r);
      set<ll> s{0}; // span of a[l..r], growing as l decreases
      for (int l = r; l >= 0; l--) {
        set<ll> t = s;
        for (ll v : s) t.insert(v ^ a[l]);
        s = t;
        CHECK_EQ(pb.max_xor(l), *s.rbegin());
        XorBasis sub; // the slots with pos >= l are a basis
        rep(b, 0, 60)
          if (pb.p[b] && pb.pos[b] >= l) CHECK(sub.insert(pb.p[b]));
        CHECK_EQ(1LL << sub.rank, (ll)sz(s));
        rep(k, 0, sz(s)) CHECK_EQ(sub.kth(k), *next(s.begin(), k));
      }
    }
  }

  int n = 100000;
  vector<ll> a(n);
  for (ll& x : a) x = rnd(0, (1LL << 60) - 1);
  double t0 = elapsed();
  XorBasis xb;
  PrefixXorBasis pb;
  ll sum = 0;
  rep(r, 0, n) {
    xb.insert(a[r]), pb.insert(a[r], r);
    sum += pb.max_xor((int)rnd(0, r)) & 1;
  }
  CHECK_EQ(pb.max_xor(0), xb.max_xor());
  CHECK_EQ(xb.max_xor(), (1LL << 60) - 1);
  fprintf(stderr, "  1e5 inserts + range queries: %.2fs (%lld)\n",
          elapsed() - t0, sum);
  puts("OK linalg/xor-basis");
}
