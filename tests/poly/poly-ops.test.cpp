#include "../testlib.hpp"
#include "../../content/poly/poly-ops.hpp"

// Independent O(n^2) series arithmetic modulo x^n.
vl trunc(vl a, int n) {
  a.resize(n);
  return a;
}
ll at(const vl& a, int i) { return i < sz(a) ? a[i] : 0; }
vl mul_trunc(const vl& a, const vl& b, int n) {
  vl c(n);
  rep(i, 0, min(sz(a), n))
    rep(j, 0, min(sz(b), n - i)) c[i + j] = (c[i + j] + a[i] * b[j]) % MOD;
  return c;
}
vl brute_inv(const vl& a, int n) {
  vl b(n);
  rep(i, 0, n) {
    ll s = i == 0;
    rep(j, 1, i + 1) s = (s - at(a, j) * b[i - j]) % MOD;
    b[i] = (s + MOD) * modinv(a[0]) % MOD;
  }
  return b;
}
// i L_i = i a_i - sum_(0 < j < i) j L_j a_(i - j), from a' = a L'.
vl brute_log(const vl& a, int n) {
  vl L(n);
  rep(i, 1, n) {
    ll s = at(a, i) * i % MOD;
    rep(j, 1, i) s = (s - j * L[j] % MOD * at(a, i - j)) % MOD;
    L[i] = (s + MOD) * modinv(i) % MOD;
  }
  return L;
}
// i e_i = sum_(0 < j <= i) j a_j e_(i - j), from E' = A' E.
vl brute_exp(const vl& a, int n) {
  vl e(n);
  rep(i, 0, n) {
    ll s = i == 0;
    rep(j, 1, i + 1) s = (s + j * at(a, j) % MOD * e[i - j]) % MOD;
    e[i] = i ? s * modinv(i) % MOD : s;
  }
  return e;
}
vl brute_pow(vl a, ll k, int n) {
  vl r = trunc({1}, n);
  for (; k; k >>= 1, a = mul_trunc(a, a, n))
    if (k & 1) r = mul_trunc(r, a, n);
  return r;
}
pair<vl, vl> brute_divmod(vl a, const vl& b) { // long division
  int n = sz(a), m = sz(b);
  vl q(max(0, n - m + 1));
  for (int i = n - m; i >= 0; i--) {
    q[i] = a[i + m - 1] * modinv(b.back()) % MOD;
    rep(j, 0, m) a[i + j] = (a[i + j] - q[i] * b[j] % MOD + MOD) % MOD;
  }
  return {q, trunc(a, m - 1)};
}

// The online convolution of the techniques note (cdq as printed
// there), with a linear step f_t = h_t + acc_t or a non-linear one.
namespace relaxed {
vl f, g, acc, h;
bool squared;
ll finish(int t, ll a) { return (h[t] + (squared ? a * a : a)) % MOD; }
void cdq(int l, int r) { // f, g, acc: vl of length n
  if (r - l == 1) return void(f[l] = finish(l, acc[l]));
  int m = (l + r) / 2;
  cdq(l, m);
  vl c =
      conv(vl(f.begin() + l, f.begin() + m), vl(g.begin(), g.begin() + r - l));
  rep(t, m, r) acc[t] = (acc[t] + c[t - l]) % MOD;
  cdq(m, r);
}
vl run(const vl& h_, const vl& g_, bool sq) {
  h = h_, g = g_, squared = sq;
  f = acc = vl(sz(h));
  if (sz(h)) cdq(0, sz(h));
  return f;
}
vl brute(const vl& h_, const vl& g_, bool sq) {
  h = h_, squared = sq;
  int n = sz(h);
  vl r(n);
  rep(t, 0, n) {
    ll a = 0;
    rep(i, 1, t + 1) a = (a + g_[i] * r[t - i]) % MOD;
    r[t] = finish(t, a);
  }
  return r;
}
// The note's closed form for the linear step: f = h / (1 - g) with
// g_0 set to 0, since acc never reads g_0.
vl closed(const vl& h_, const vl& g_) {
  int n = sz(h_);
  vl d{1};
  rep(i, 1, n) d.push_back((MOD - g_[i]) % MOD);
  return trunc(conv(h_, poly_inv(d, n)), n);
}
} // namespace relaxed

vl random_poly(int len, int it) {
  vl a = rnd_vec(len, 0, it % 4 ? MOD - 1 : 2);
  if (it % 9 == 0) fill(all(a), MOD - 1);
  return a;
}

int main() {
  rep(it, 0, iters(1200)) {
    int n = (int)rnd(0, it % 5 ? 40 : 140);
    vl a = random_poly((int)rnd(1, n + 3), it);
    a[0] = rnd(1, MOD - 1);
    CHECK(poly_inv(a, n) == brute_inv(a, n));
    a[0] = 1;
    vl L = poly_log(a, n);
    CHECK(L == brute_log(a, n));
    CHECK(poly_exp(L, n) == trunc(a, n));
    // Square roots from the Tweak: a_0 = 1, and r^2 x^(2w) a with a
    // known r (as mod_sqrt would return).
    vl s = poly_pow(a, (MOD + 1) / 2, n);
    CHECK(mul_trunc(s, s, n) == trunc(a, n));
    int w = (int)rnd(0, 3);
    ll r = rnd(1, MOD - 1), c = r * r % MOD;
    vl sq(w * 2);
    for (ll x : a) sq.push_back(x * c % MOD);
    vl root(w);
    for (ll x : poly_pow(trunc(a, n), (MOD + 1) / 2, n))
      root.push_back(x * r % MOD);
    CHECK(mul_trunc(root, root, n) == trunc(sq, n));
    a[0] = 0;
    CHECK(poly_exp(a, n) == brute_exp(a, n));
    // pow: leading zeros, any leading coefficient, small and huge k.
    int m = min(n, 40), zeros = (int)rnd(0, 4);
    vl p = random_poly((int)rnd(1, m + 3), it);
    rep(i, 0, min(zeros, sz(p))) p[i] = 0;
    if (it % 11 == 0) fill(all(p), 0);
    ll k = it % 3 ? rnd(0, 6) : rnd(0, (ll)1e18);
    if (it % 7 == 0 && zeros && m) k = (m - 1) / zeros + rnd(0, 1);
    CHECK(poly_pow(p, k, m) == brute_pow(p, k, m));
    // divmod, including deg a < deg b and constant b.
    vl x = random_poly((int)rnd(0, 60), it);
    vl y = random_poly((int)rnd(1, it % 4 ? 30 : 70), it);
    if (y.back() == 0) y.back() = 1;
    CHECK(poly_divmod(x, y) == brute_divmod(x, y));
  }
  // Exactly n terms back, also for n = 0 and n = 1.
  rep(n, 0, 3) {
    CHECK(poly_inv({5}, n) == brute_inv({5}, n));
    CHECK(poly_log({1}, n) == vl(n));
    CHECK(poly_exp({0}, n) == trunc({1}, n));
    CHECK(poly_pow({0}, 0, n) == trunc({1}, n));
    CHECK(poly_pow({0, 3}, 5, n) == vl(n));
    vl late(12);
    late[11] = 7; // v = 11: v k overflows ll for k = 1e18
    CHECK(poly_pow(late, (ll)1e18, n) == vl(n));
  }
  // The header's partition example against Euler's pentagonal
  // recurrence p(m) = sum_(k >= 1) (-1)^(k + 1) (p(m - k(3k - 1)/2)
  // + p(m - k(3k + 1)/2)).
  int N = 1 << 14;
  vl Lp(N), inv(N, 1), want(N);
  rep(i, 2, N) inv[i] = (MOD - MOD / i) * inv[MOD % i] % MOD;
  rep(k, 1, N)
    for (int j = 1; j * k < N; j++) Lp[j * k] = (Lp[j * k] + inv[j]) % MOD;
  want[0] = 1;
  rep(i, 1, N)
    for (int k = 1; k * (3 * k - 1) / 2 <= i; k++) {
      int g = k * (3 * k - 1) / 2, h = g + k;
      ll s = want[i - g] + (h <= i ? want[i - h] : 0);
      want[i] = (want[i] + (k % 2 ? s : 2 * MOD - s)) % MOD;
    }
  CHECK(poly_exp(Lp, N) == want);
  // Online convolution (techniques note), g_0 != 0 included: the
  // reviewer's case g = 5 1 2 3 4 5 gives 1 1 3 8 21 55.
  vl h0{1, 0, 0, 0, 0, 0}, g0{5, 1, 2, 3, 4, 5};
  CHECK(relaxed::run(h0, g0, 0) == vl({1, 1, 3, 8, 21, 55}));
  CHECK(relaxed::closed(h0, g0) == vl({1, 1, 3, 8, 21, 55}));
  rep(it, 0, iters(300)) {
    int n = (int)rnd(0, it % 10 ? 40 : 300);
    vl h = random_poly(n, it), g = random_poly(n, it + 1);
    if (n && it % 3 == 0) g[0] = 1;
    CHECK(relaxed::run(h, g, 0) == relaxed::brute(h, g, 0));
    CHECK(relaxed::closed(h, g) == relaxed::brute(h, g, 0));
    CHECK(relaxed::run(h, g, 1) == relaxed::brute(h, g, 1));
  }
  // Large: identities at n = 2^16, pow and divmod at 3e4.
  int n = 1 << 16;
  vl a = rnd_vec(n, 0, MOD - 1);
  a[0] = 1;
  double t0 = elapsed();
  CHECK(trunc(conv(a, poly_inv(a, n)), n) == trunc({1}, n));
  CHECK(poly_exp(poly_log(a, n), n) == a);
  a.resize(n = 30000);
  vl sq = trunc(conv(a, a), n);
  CHECK(poly_pow(a, 2, n) == sq);
  CHECK(poly_pow(sq, (MOD + 1) / 2, n) == a);
  vl b = trunc(a, n / 2);
  if (!b.back()) b.back() = 1;
  auto [q, r] = poly_divmod(a, b);
  vl back = conv(q, b);
  rep(i, 0, sz(r)) back[i] = (back[i] + r[i]) % MOD;
  CHECK(back == a);
  fprintf(stderr, "  large: %.2fs\n", elapsed() - t0);
  puts("OK poly/poly-ops");
}
