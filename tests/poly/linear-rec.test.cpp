#include "../testlib.hpp"
#include "../../content/poly/linear-rec.hpp"

// Does some recurrence of order L fit all of s? Gaussian elimination
// over F_MOD on the unknowns c_1..c_L.
bool fits(const vl& s, int L) {
  vector<vl> A;
  rep(i, L, sz(s)) {
    vl row(L + 1);
    rep(j, 0, L) row[j] = s[i - 1 - j];
    row[L] = s[i];
    A.push_back(row);
  }
  int r = 0;
  rep(col, 0, L) {
    int p = r;
    while (p < sz(A) && !A[p][col]) p++;
    if (p == sz(A)) continue;
    swap(A[r], A[p]);
    ll inv = modinv(A[r][col]);
    for (ll& x : A[r]) x = x * inv % MOD;
    rep(i, 0, sz(A))
      if (i != r && A[i][col]) {
        ll f = A[i][col];
        rep(j, 0, L + 1) A[i][j] = ((A[i][j] - f * A[r][j]) % MOD + MOD) % MOD;
      }
    r++;
  }
  rep(i, r, sz(A))
    if (A[i][L]) return false;
  return true;
}

int min_order(const vl& s) {
  int L = 0;
  while (!fits(s, L)) L++;
  return L;
}

bool generates(const vl& s, const vl& c) {
  rep(i, sz(c), sz(s)) {
    ll x = 0;
    rep(j, 0, sz(c)) x = (x + c[j] * s[i - 1 - j]) % MOD;
    if (x != s[i]) return false;
  }
  return true;
}

// s_k from the companion matrix power, independent of linear_rec.
ll matrix_kth(const vl& s, const vl& c, ll k) {
  int L = sz(c);
  if (k < L) return s[k];
  auto mul = [L](const vl& A, const vl& B) { // L x L, row-major
    vl C(L * L);
    rep(i, 0, L)
      rep(t, 0, L)
        rep(j, 0, L)
          C[i * L + j] = (C[i * L + j] + A[i * L + t] * B[t * L + j]) % MOD;
    return C;
  };
  vl M(L * L), P(L * L);
  rep(j, 0, L) M[j] = c[j];
  rep(i, 1, L) M[i * L + i - 1] = 1;
  rep(i, 0, L) P[i * L + i] = 1;
  for (ll e = k - L + 1; e; e /= 2, M = mul(M, M))
    if (e & 1) P = mul(P, M);
  ll res = 0; // state (s_(L-1), ..., s_0) maps to (s_k, ...)
  rep(j, 0, L) res = (res + P[j] * s[L - 1 - j]) % MOD;
  return res;
}

vl extend(vl s, const vl& c, int n) {
  while (sz(s) < n) {
    ll x = 0;
    rep(j, 0, sz(c)) x = (x + c[j] * s[sz(s) - 1 - j]) % MOD;
    s.push_back(x);
  }
  return s;
}

int main() {
  // Arbitrary short sequences: BM is as short as the brute force.
  rep(it, 0, iters(1000)) {
    int n = (int)rnd(0, 12);
    vl s = rnd_vec(n, 0, it % 2 ? 2 : MOD - 1);
    vl c = berlekamp_massey(s);
    CHECK_EQ(sz(c), min_order(s));
    CHECK(generates(s, c));
  }
  // Planted recurrences, k-th term against generation and matrices.
  rep(it, 0, iters(400)) {
    int L = (int)rnd(0, 8);
    vl c = rnd_vec(L, 0, it % 3 ? MOD - 1 : 2);
    vl s = extend(rnd_vec(L, 0, MOD - 1), c, 2 * L + (int)rnd(0, 4));
    vl got = berlekamp_massey(s);
    CHECK(sz(got) <= L && generates(s, got));
    vl full = extend(s, c, 300);
    CHECK(generates(full, got));
    rep(t, 0, 20) {
      int k = t < 10 ? t : (int)rnd(0, 299);
      CHECK_EQ(linear_rec(s, got, k), full[k]);
    }
    ll k = rnd(0, (ll)1e18);
    CHECK_EQ(linear_rec(s, got, k), matrix_kth(s, got, k));
    CHECK_EQ(linear_rec(s, c, k), matrix_kth(s, c, k));
    // Hand-written input may use negative representatives.
    vl sn = s, cn = c;
    for (ll& x : sn) x -= rnd(0, 1) ? MOD : 0;
    for (ll& x : cn) x -= rnd(0, 1) ? MOD : 0;
    CHECK(berlekamp_massey(sn) == got);
    rep(t, 0, 5) {
      int kk = (int)rnd(0, 299);
      CHECK_EQ(linear_rec(sn, cn, kk), full[kk]);
    }
    CHECK_EQ(linear_rec(sn, cn, k), matrix_kth(s, c, k));
  }
  CHECK(berlekamp_massey({}).empty());
  CHECK(berlekamp_massey({0, 0, 0}).empty());
  CHECK_EQ(sz(berlekamp_massey({0, 0, 1})), 3);
  CHECK(berlekamp_massey({1, 1, 2, 3, 5, 8}) == vl({1, 1}));
  vl pw(10); // 2^k + 3^k: orders add under sums
  rep(k, 0, 10) pw[k] = (modpow(2, k) + modpow(3, k)) % MOD;
  CHECK(berlekamp_massey(pw) == vl({5, MOD - 6}));
  CHECK_EQ(linear_rec({}, {}, (ll)1e18), 0LL);
  CHECK_EQ(linear_rec({7}, {0}, 0), 7LL);
  CHECK_EQ(linear_rec({7}, {0}, 1), 0LL);
  CHECK_EQ(linear_rec({0, 1}, {1, 1}, 90), 2'880'067'194'370'816'120 % MOD);
  CHECK_EQ(linear_rec({0, 1}, {2, -1}, (ll)1e18), (ll)1e18 % MOD);
  CHECK_EQ(linear_rec({1, 1}, {-1, 0}, 10), MOD - 1);
  CHECK_EQ(linear_rec({-3}, {1}, 0), MOD - 3);
  // Large: order 500 from 1010 terms, then far ahead.
  int L = 500;
  vl c = rnd_vec(L, 0, MOD - 1);
  vl s = extend(rnd_vec(L, 0, MOD - 1), c, 2 * L + 10);
  CHECK(berlekamp_massey(s) == c);
  ll want = extend(s, c, 3000)[2999];
  CHECK_EQ(linear_rec(s, c, 2999), want);
  vl period(200); // s_i = s_(i - 200)
  period.back() = 1;
  CHECK_EQ(linear_rec(s, period, (ll)1e18), s[(ll)1e18 % 200]);
  puts("OK poly/linear-rec");
}
