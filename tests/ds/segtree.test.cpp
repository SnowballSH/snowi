#include "../testlib.hpp"

// The snippet fixes T = ll. Included in a template whose parameter
// is named ll, the same code also runs over Seq below, whose op is
// not commutative, so swapped operands are caught.
template<class ll> struct Generic {
#include "../../content/ds/segtree.hpp"
};
using SegTree = Generic<long long>::SegTree;

// Free monoid on digits 0..3: the exact sequence, 2 bits a digit,
// at most 30 digits. Seq(0) is the empty sequence (the snippet's e).
struct Seq {
  ll code = 0;
  int len = 0;
  constexpr Seq(ll = 0) {}
  friend Seq operator+(Seq a, Seq b) {
    Seq s;
    s.code = a.code << 2 * b.len | b.code, s.len = a.len + b.len;
    return s;
  }
  bool operator==(const Seq&) const = default;
};
using SeqTree = Generic<Seq>::SegTree;

Seq leaf(int v) {
  Seq s;
  s.code = v, s.len = 1;
  return s;
}

int runs(Seq s) { // maximal blocks of equal digits
  int r = s.len > 0;
  rep(i, 1, s.len) r += (s.code >> 2 * i & 3) != (s.code >> 2 * (i - 1) & 3);
  return r;
}

Seq fold(const vector<Seq>& a, int l, int r) {
  Seq s;
  rep(i, l, r) s = s + a[i];
  return s;
}

void stress_seq() {
  rep(it, 0, iters(800)) {
    int n = (int)rnd(0, 20);
    vector<Seq> a(n);
    for (auto& x : a) x = leaf((int)rnd(0, 2));
    SeqTree st(a);
    rep(op, 0, 30) {
      if (n && rnd(0, 2) == 0) {
        int i = (int)rnd(0, n - 1);
        st.set(i, a[i] = leaf((int)rnd(0, 2)));
      }
      int l = (int)rnd(0, n), r = (int)rnd(l, n), i = (int)rnd(0, n);
      CHECK(st.query(l, r) == fold(a, l, r));
      CHECK(st.query_all() == fold(a, 0, n));
      if (n) CHECK(st.get(i % n) == a[i % n]);
      // The walks may only test folds of [i, j) and [j, i), in order.
      int k = (int)rnd(0, 5);
      auto right = [&](Seq s) {
        CHECK(i + s.len <= n && s == fold(a, i, i + s.len));
        return runs(s) <= k;
      };
      auto left = [&](Seq s) {
        CHECK(s.len <= i && s == fold(a, i - s.len, i));
        return runs(s) <= k;
      };
      int want = i;
      while (want < n && runs(fold(a, i, want + 1)) <= k) want++;
      CHECK_EQ(st.max_right(i, right), want);
      want = i;
      while (want > 0 && runs(fold(a, want - 1, i)) <= k) want--;
      CHECK_EQ(st.min_left(i, left), want);
    }
  }
}

ll brute_sum(const vector<ll>& a, int l, int r) {
  return accumulate(a.begin() + l, a.begin() + r, 0LL);
}

int brute_max_right(const vector<ll>& a, int l, auto pred) {
  ll acc = 0;
  int r = l;
  while (r < sz(a) && pred(acc + a[r])) acc += a[r++];
  return r;
}

int brute_min_left(const vector<ll>& a, int r, auto pred) {
  ll acc = 0;
  int l = r;
  while (l > 0 && pred(acc + a[l - 1])) acc += a[--l];
  return l;
}

// Non-negative values make prefix predicates monotone; the other
// modes test plain sums with negative and huge values.
ll rnd_value(int mode) {
  if (mode == 0) return rnd(0, 5);
  if (mode == 1) return rnd(-5, 5);
  return rnd(-4e17, 4e17);
}

void check_random(SegTree& st, const vector<ll>& a, int mode) {
  int n = sz(a);
  int l = (int)rnd(0, n), r = (int)rnd(l, n), i = (int)rnd(0, n);
  CHECK_EQ(st.query(l, r), brute_sum(a, l, r));
  CHECK_EQ(st.query_all(), brute_sum(a, 0, n));
  if (n) CHECK_EQ(st.get(i % n), a[i % n]);
  if (mode != 0) return;
  ll x = rnd(-1, 25);
  auto lt = [&](ll s) { return s < x; };
  auto le = [&](ll s) { return s <= x; };
  if (x >= 0) { // pred(e) must hold
    CHECK_EQ(st.max_right(i, le), brute_max_right(a, i, le));
    CHECK_EQ(st.min_left(i, le), brute_min_left(a, i, le));
  }
  if (x > 0) {
    CHECK_EQ(st.max_right(i, lt), brute_max_right(a, i, lt));
    CHECK_EQ(st.min_left(i, lt), brute_min_left(a, i, lt));
  }
}

int main() {
  rep(it, 0, iters(1500)) {
    int n = (int)rnd(0, 20), mode = (int)rnd(0, 3) % 3;
    vector<ll> a(n);
    for (auto& x : a) x = rnd_value(mode);
    SegTree st(a);
    rep(op, 0, 30) {
      if (n && rnd(0, 2) == 0) {
        int i = (int)rnd(0, n - 1);
        a[i] = rnd_value(mode);
        st.set(i, a[i]);
      }
      check_random(st, a, mode);
    }
  }
  rep(n, 0, 6) {
    SegTree st(n);
    vector<ll> a(n);
    rep(it, 0, 20) check_random(st, a, 0);
  }
  stress_seq();
  // Large: k-th one of a 0/1 array under toggles, against counts
  // per block of 256.
  int n = 100'000, B = 256;
  vector<ll> a(n);
  vi block(n / B + 1);
  rep(i, 0, n)
    if (rnd(0, 1)) a[i] = 1, block[i / B]++;
  auto prefix = [&](int r) { // ones in [0, r)
    int c = 0, i = 0;
    for (; i + B <= r; i += B) c += block[i / B];
    for (; i < r; i++) c += (int)a[i];
    return c;
  };
  auto kth = [&](int k) { // position of the k-th one, or n
    int i = 0;
    for (; i < n && block[i / B] < k; i += B) k -= block[i / B];
    for (; i < n; i++)
      if ((k -= (int)a[i]) == 0) return i;
    return n;
  };
  SegTree st(a);
  rep(it, 0, n) {
    int i = (int)rnd(0, n - 1);
    block[i / B] += a[i] ? -1 : 1;
    a[i] ^= 1, st.set(i, a[i]);
    int k = (int)rnd(1, n / 2 + 10);
    CHECK_EQ(st.max_right(0, [&](ll s) { return s < k; }), kth(k));
    int r = (int)rnd(0, n), cnt = prefix(r);
    int want = cnt < k ? 0 : kth(cnt - k + 1) + 1;
    CHECK_EQ(st.min_left(r, [&](ll s) { return s < k; }), want);
  }
  CHECK_EQ(st.query_all(), (ll)prefix(n));
  puts("OK ds/segtree");
}
