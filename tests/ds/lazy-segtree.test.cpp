#include "../testlib.hpp"

// The snippet fixes S = F = ll. Included in a template whose parameter
// is named ll, the same code also runs over W<P> below, which carries a
// node and a tag of a policy P with non-commutative op or composition.
template<class ll> struct Generic {
#include "../../content/ds/lazy-segtree.hpp"
};
using LazySegTree = Generic<long long>::LazySegTree;

// The snippet's min(a, b), s + f and f + g become P::op, P::mapping and
// P::composition on the matching part; the other part is ignored. Its
// e = LLONG_MAX and id = 0 both become W(): P::e and P::id.
template<class P> struct W {
  P::S s = P::e;
  P::F f = P::id;
  constexpr W(ll = 0) {}
  friend W min(W a, W b) {
    W w;
    w.s = P::op(a.s, b.s);
    return w;
  }
  friend W operator+(W x, W y) {
    W w;
    w.s = P::mapping(y.f, x.s), w.f = P::composition(x.f, y.f);
    return w;
  }
};

// Digits 0..4 as an exact sequence (3 bits a digit, at most 20), tags
// x -> b x + c mod 5: both op and composition are non-commutative.
struct SeqAffine {
  struct S {
    ll code;
    int len;
    bool operator==(const S&) const = default;
  };
  struct F {
    int b, c;
  };
  static constexpr S e{0, 0};
  static constexpr F id{1, 0};
  static S op(S x, S y) {
    return {x.code << 3 * y.len | y.code, x.len + y.len};
  }
  static S mapping(F f, S s) {
    S r{0, s.len};
    for (int i = s.len; i--;)
      r.code = r.code << 3 | (f.b * (s.code >> 3 * i & 7) + f.c) % 5;
    return r;
  }
  static F composition(F f, F g) { // g, then f
    return {f.b * g.b % 5, (f.b * g.c + f.c) % 5};
  }
  static S leaf() { return {rnd(0, 4), 1}; }
  static F tag() { return {(int)rnd(0, 4), (int)rnd(0, 4)}; }
  static constexpr int kmin = 0; // pred(e, k) holds for k >= kmin
  static bool pred(S s, int k) { // at most k runs of equal digits
    int runs = s.len > 0;
    rep(i, 1, s.len)
      runs += (s.code >> 3 * i & 7) != (s.code >> 3 * (i - 1) & 7);
    return runs <= k;
  }
};

// Range assign with id as "none", range max: the Tweak's recipe.
struct AssignMax {
  typedef ll S, F;
  static constexpr S e = LLONG_MIN;
  static constexpr F id = LLONG_MIN;
  static S op(S x, S y) { return max(x, y); }
  static S mapping(F f, S s) { return f == id ? s : f; }
  static F composition(F f, F g) { return f == id ? g : f; }
  static S leaf() { return rnd(-5, 5); }
  static F tag() { return rnd(0, 4) ? rnd(-5, 5) : id; }
  static constexpr int kmin = -6;
  static bool pred(S s, int k) { return s < k; }
};

template<class P> using Tree = Generic<W<P>>::LazySegTree;

template<class P> W<P> node(typename P::S s) {
  W<P> w;
  w.s = s;
  return w;
}

template<class P> W<P> tag(typename P::F f) {
  W<P> w;
  w.f = f;
  return w;
}

template<class P> P::S fold(const vector<W<P>>& a, int l, int r) {
  auto s = P::e;
  rep(i, l, r) s = P::op(s, a[i].s);
  return s;
}

template<class P> void stress_policy(int rounds) {
  rep(it, 0, iters(rounds)) {
    int n = (int)rnd(0, 20);
    vector<W<P>> a(n);
    for (auto& x : a) x = node<P>(P::leaf());
    Tree<P> st(a);
    rep(op, 0, 30) {
      if (n && rnd(0, 3) == 0) {
        int i = (int)rnd(0, n - 1);
        st.set(i, a[i] = node<P>(P::leaf()));
      } else { // several applies between queries stack up tags
        int l = (int)rnd(0, n), r = (int)rnd(l, n);
        auto f = P::tag();
        rep(i, l, r) a[i].s = P::mapping(f, a[i].s);
        st.apply(l, r, tag<P>(f));
      }
      int l = (int)rnd(0, n), r = (int)rnd(l, n), i = (int)rnd(0, n);
      CHECK(st.query(l, r).s == fold(a, l, r));
      CHECK(st.query_all().s == fold(a, 0, n));
      if (n) CHECK(st.get(i % n).s == a[i % n].s);
      // The walks may only test folds of [i, j) and [j, i), in order.
      int k = (int)rnd(P::kmin, 6), want;
      auto right = [&](W<P> w) {
        bool ok = false;
        rep(j, i, n + 1) ok |= w.s == fold(a, i, j);
        CHECK(ok);
        return P::pred(w.s, k);
      };
      auto left = [&](W<P> w) {
        bool ok = false;
        rep(j, 0, i + 1) ok |= w.s == fold(a, j, i);
        CHECK(ok);
        return P::pred(w.s, k);
      };
      for (want = i; want < n && P::pred(fold(a, i, want + 1), k);) want++;
      CHECK_EQ(st.max_right(i, right), want);
      for (want = i; want > 0 && P::pred(fold(a, want - 1, i), k);) want--;
      CHECK_EQ(st.min_left(i, left), want);
    }
  }
}

ll brute_min(const vector<ll>& a, int l, int r) {
  ll m = LLONG_MAX;
  rep(i, l, r) m = min(m, a[i]);
  return m;
}

int brute_max_right(const vector<ll>& a, int l, auto pred) {
  ll acc = LLONG_MAX;
  int r = l;
  while (r < sz(a) && pred(min(acc, a[r]))) acc = min(acc, a[r++]);
  return r;
}

int brute_min_left(const vector<ll>& a, int r, auto pred) {
  ll acc = LLONG_MAX;
  int l = r;
  while (l > 0 && pred(min(acc, a[l - 1]))) acc = min(acc, a[--l]);
  return l;
}

ll rnd_value() {
  if (rnd(0, 30) == 0) return rnd(0, 1) ? ll(1e18) : ll(-1e18);
  return rnd(-5, 5);
}

void check_random(LazySegTree& st, const vector<ll>& a) {
  int n = sz(a);
  int l = (int)rnd(0, n), r = (int)rnd(l, n), i = (int)rnd(0, n);
  CHECK_EQ(st.query(l, r), brute_min(a, l, r));
  CHECK_EQ(st.query_all(), brute_min(a, 0, n));
  if (n) CHECK_EQ(st.get(i % n), a[i % n]);
  ll x = rnd(-8, 8);
  auto gt = [&](ll v) { return v > x; };
  auto ge = [&](ll v) { return v >= x; };
  CHECK_EQ(st.max_right(i, gt), brute_max_right(a, i, gt));
  CHECK_EQ(st.max_right(i, ge), brute_max_right(a, i, ge));
  CHECK_EQ(st.min_left(i, gt), brute_min_left(a, i, gt));
  CHECK_EQ(st.min_left(i, ge), brute_min_left(a, i, ge));
}

int main() {
  { // assign: push hands id to mapping, which must return s then
    vector<W<AssignMax>> a;
    for (ll v : {1, 2, 3, 4}) a.push_back(node<AssignMax>(v));
    Tree<AssignMax> st(a);
    st.apply(0, 2, tag<AssignMax>(7));
    CHECK_EQ(st.query(0, 4).s, 7);
    CHECK_EQ(st.get(3).s, 4);
    CHECK_EQ(st.query(2, 4).s, 4);
  }
  stress_policy<SeqAffine>(600);
  stress_policy<AssignMax>(600);
  rep(it, 0, iters(1500)) {
    int n = (int)rnd(0, 20);
    vector<ll> a(n);
    for (auto& x : a) x = rnd_value();
    LazySegTree st(a);
    rep(op, 0, 30) {
      int kind = (int)rnd(0, 3);
      if (n && kind == 0) {
        int i = (int)rnd(0, n - 1);
        a[i] = rnd_value();
        st.set(i, a[i]);
      } else if (n && kind == 1) {
        int i = (int)rnd(0, n - 1);
        ll f = rnd(-3, 3);
        a[i] += f, st.apply(i, f);
      } else {
        int l = (int)rnd(0, n), r = (int)rnd(l, n);
        ll f = rnd(-3, 3);
        rep(i, l, r) a[i] += f;
        st.apply(l, r, f);
      }
      check_random(st, a);
    }
  }
  // 1e5 range adds and queries; a Fenwick tree of differences
  // gives point values, the final array is checked in full.
  int n = 100'000;
  vector<ll> a(n), bit(n + 1);
  auto bit_add = [&](int i, ll f) {
    for (i++; i <= n; i += i & -i) bit[i] += f;
  };
  auto value = [&](int i) {
    ll s = a[i];
    for (i++; i; i -= i & -i) s += bit[i];
    return s;
  };
  rep(i, 0, n) a[i] = rnd(-1e9, 1e9);
  LazySegTree st(a);
  rep(it, 0, n) {
    int l = (int)rnd(0, n - 1), r = (int)rnd(l + 1, n);
    ll f = rnd(-1e9, 1e9);
    st.apply(l, r, f), bit_add(l, f), bit_add(r, -f);
    int i = (int)rnd(0, n - 1);
    CHECK_EQ(st.query(i, i + 1), value(i));
    l = (int)rnd(0, n - 1), r = min(n, l + (int)rnd(1, 8));
    ll want = LLONG_MAX;
    rep(j, l, r) want = min(want, value(j));
    CHECK_EQ(st.query(l, r), want);
  }
  vector<ll> b(n);
  rep(i, 0, n) b[i] = value(i);
  rep(i, 0, n) CHECK_EQ(st.get(i), b[i]);
  CHECK_EQ(st.query_all(), *min_element(all(b)));
  puts("OK ds/lazy-segtree");
}
