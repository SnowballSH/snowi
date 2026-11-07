#include "../testlib.hpp"
#include "../../content/poly/subset-transforms.hpp"

// Definitions, O(4^n).
vl brute_conv(const vl& a, const vl& b, char op) {
  int n = sz(a);
  vl c(n);
  rep(i, 0, n)
    rep(j, 0, n) {
      int k = op == '|' ? i | j : op == '&' ? i & j : i ^ j;
      c[k] = (c[k] + a[i] * b[j]) % MOD;
    }
  return c;
}
vl brute_subset(const vl& a, const vl& b) {
  int n = sz(a);
  vl c(n);
  rep(i, 0, n)
    rep(j, 0, n)
      if (!(i & j)) c[i | j] = (c[i | j] + a[i] * b[j]) % MOD;
  return c;
}
// sum over T with T subset of S ('|'), superset ('&'), or
// (-1)^|S & T| a_T ('^').
vl brute_transform(const vl& a, char op) {
  int n = sz(a);
  vl r(n);
  rep(S, 0, n)
    rep(T, 0, n) {
      if (op == '|' && (T & ~S)) continue;
      if (op == '&' && (S & ~T)) continue;
      ll x = op == '^' && popcount(unsigned(S & T)) % 2 ? MOD - a[T] : a[T];
      r[S] = (r[S] + x) % MOD;
    }
  return r;
}

// bit_transform with the Tweak applied: exact integers (the last
// two lines replaced by the XOR division) and max over subsets or
// supersets (the last two lines deleted).
void exact_transform(vl& a, char op, bool inv) {
  int n = sz(a);
  for (int s = 1; s < n; s *= 2)
    for (int i = 0; i < n; i += 2 * s)
      rep(j, i, i + s) {
        ll &u = a[j], &v = a[j + s];
        if (op == '|') v += inv ? -u : u;
        else if (op == '&') u += inv ? -v : v;
        else tie(u, v) = pair(u + v, u - v);
      }
  if (inv && op == '^')
    for (ll& x : a) x /= n;
}
void max_transform(vl& a, char op) {
  int n = sz(a);
  for (int s = 1; s < n; s *= 2)
    for (int i = 0; i < n; i += 2 * s)
      rep(j, i, i + s) {
        ll &u = a[j], &v = a[j + s];
        if (op == '|') v = max(v, u);
        else u = max(u, v);
      }
}

int main() {
  rep(it, 0, iters(2000)) {
    int L = (int)rnd(0, it % 5 ? 4 : 7), n = 1 << L;
    ll hi = it % 3 ? MOD - 1 : 2;
    vl a = rnd_vec(n, 0, hi), b = rnd_vec(n, 0, hi);
    if (it % 11 == 0) fill(all(a), MOD - 1), fill(all(b), MOD - 1);
    for (char op : {'|', '&', '^'}) {
      vl t = a;
      bit_transform(t, op, 0);
      CHECK(t == brute_transform(a, op));
      bit_transform(t, op, 1);
      CHECK(t == a);
      CHECK(bit_conv(a, b, op) == brute_conv(a, b, op));
    }
    CHECK(subset_conv(a, b) == brute_subset(a, b));
  }
  // The Tweak: exact integers, negative and above MOD, and max.
  rep(it, 0, iters(500)) {
    int n = 1 << rnd(0, 6);
    ll hi = it % 2 ? (ll)4e9 : 3;
    vl a = rnd_vec(n, -hi, hi), b = rnd_vec(n, -1e5, 1e5);
    for (char op : {'|', '&', '^'}) {
      vl t = a, mx = a, x = a, y = b, c(n);
      exact_transform(t, op, 0);
      if (op != '^') max_transform(mx, op);
      rep(S, 0, n) {
        ll want = 0, best = LLONG_MIN;
        rep(T, 0, n) {
          bool in = op == '|' ? !(T & ~S) : !(S & ~T);
          if (op == '^') want += popcount(unsigned(S & T)) % 2 ? -a[T] : a[T];
          else if (in) want += a[T], best = max(best, a[T]);
        }
        CHECK_EQ(t[S], want);
        if (op != '^') CHECK_EQ(mx[S], best);
      }
      exact_transform(t, op, 1);
      CHECK(t == a);
      exact_transform(x, op, 0), exact_transform(y, op, 0);
      rep(i, 0, n) x[i] *= y[i]; // bit_conv without % MOD
      exact_transform(x, op, 1);
      rep(i, 0, n)
        rep(j, 0, n) {
          int k = op == '|' ? i | j : op == '&' ? i & j : i ^ j;
          c[k] += a[i] * b[j];
        }
      CHECK(x == c);
    }
  }
  // Identities from the Notes: [S = all] for AND, [S = 0] for OR.
  vl f(8), e(8), g = rnd_vec(8, 1, MOD - 1);
  f[7] = e[0] = 1;
  CHECK(bit_conv(f, g, '&') == g && bit_conv(e, g, '|') == g);
  CHECK(bit_conv(f, g, '|') != g); // OR is not AND
  // Large: n = 2^16 transforms, 2^14 subset convolution, checked
  // at sampled masks.
  int L = 16, n = 1 << L;
  vl a = rnd_vec(n, 0, MOD - 1), b = rnd_vec(n, 0, MOD - 1);
  vl cx = bit_conv(a, b, '^'), co = bit_conv(a, b, '|');
  vl ca = bit_conv(a, b, '&');
  rep(t, 0, 20) {
    int S = (int)rnd(0, n - 1);
    ll want = 0;
    rep(i, 0, n) want = (want + a[i] * b[i ^ S]) % MOD;
    CHECK_EQ(cx[S], want);
    S &= (int)rnd(0, n - 1) & (int)rnd(0, n - 1); // few bits
    ll wor = 0, wand = 0;
    int C = (n - 1) ^ S;
    // i | j = S means i in S and j = (S ^ i) | u with u in i. By
    // complements, i & j = C for C = ~S has ~i | ~j = S.
    for (int i = S;; i = (i - 1) & S) {
      for (int u = i;; u = (u - 1) & i) {
        int j = (S ^ i) | u;
        wor = (wor + a[i] * b[j]) % MOD;
        wand = (wand + a[(n - 1) ^ i] * b[(n - 1) ^ j]) % MOD;
        if (!u) break;
      }
      if (!i) break;
    }
    CHECK_EQ(co[S], wor);
    CHECK_EQ(ca[C], wand);
  }
  int m = 1 << 14;
  a.resize(m), b.resize(m);
  vl cs = subset_conv(a, b);
  rep(t, 0, 20) {
    int S = (int)rnd(0, m - 1);
    ll want = 0;
    for (int T = S;; T = (T - 1) & S) {
      want = (want + a[T] * b[S ^ T]) % MOD;
      if (!T) break;
    }
    CHECK_EQ(cs[S], want);
  }
  puts("OK poly/subset-transforms");
}
