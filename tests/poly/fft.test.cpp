#include "../testlib.hpp"
#include "../../content/poly/fft.hpp"

vector<ll> brute_conv(const vector<ll>& a, const vector<ll>& b) {
  if (a.empty() || b.empty()) return {};
  vector<ll> c(sz(a) + sz(b) - 1);
  rep(i, 0, sz(a))
    rep(j, 0, sz(b)) c[i + j] += a[i] * b[j];
  return c;
}

vector<ll> rounded(const vd& c) {
  vector<ll> r;
  for (double x : c) r.push_back(llround(x));
  return r;
}

// a[k] <- sum_j a_j w^(jk), w = e^(2 pi i / n), natural order.
vector<cd> brute_dft(const vector<cd>& a) {
  int n = sz(a);
  vector<cd> res(n);
  rep(k, 0, n)
    rep(j, 0, n) res[k] += a[j] * polar(1.0, 2 * M_PI * ((ll)j * k % n) / n);
  return res;
}

// a = n copies of v, b = m copies of w:
// c[k] = min(k + 1, n, m, n + m - 1 - k) v w.
void check_constant(int n, int m, ll v, ll w) {
  vd c = fft_conv(vd(n, double(v)), vd(m, double(w)));
  CHECK_EQ(sz(c), n + m - 1);
  rep(k, 0, sz(c))
    CHECK_EQ(llround(c[k]), min({k + 1, n, m, n + m - 1 - k}) * v * w);
}

// a_i = (-1)^i v, b = m copies of w: the results are 0 or +-v w,
// but the rounding error follows the inputs, not the results.
void check_alternating(int n, int m, ll v, ll w) {
  vd a(n, double(v));
  for (int i = 1; i < n; i += 2) a[i] = -double(v);
  vd c = fft_conv(a, vd(m, double(w)));
  CHECK_EQ(sz(c), n + m - 1);
  rep(k, 0, sz(c)) {
    int lo = max(0, k - m + 1), hi = min(k, n - 1);
    ll want = (hi - lo) % 2 ? 0 : lo % 2 ? -v * w : v * w;
    CHECK_EQ(llround(c[k]), want);
  }
}

int main() {
  rep(it, 0, iters(1500)) {
    int n = (int)rnd(0, it % 4 ? 20 : 100);
    int m = (int)rnd(0, it % 4 ? 20 : 100);
    // Magnitudes may differ a lot, e.g. weights times a 0/1 mask.
    ll ha = it % 3 ? 1000 : 1'000'000, hb = ha;
    if (it % 4 == 1) ha = 1'000'000'000, hb = 1;
    if (it % 4 == 3) ha = 1, hb = 1'000'000'000;
    vector<ll> a = rnd_vec(n, -ha, ha), b = rnd_vec(m, -hb, hb);
    CHECK(rounded(fft_conv(vd(all(a)), vd(all(b)))) == brute_conv(a, b));
    vd x(n), y(m);
    for (double& v : x) v = rnd_real(-1, 1);
    for (double& v : y) v = rnd_real(-1, 1);
    vd c = fft_conv(x, y);
    CHECK_EQ(sz(c), n && m ? n + m - 1 : 0);
    rep(k, 0, sz(c)) {
      double want = 0;
      rep(i, max(0, k - m + 1), min(k, n - 1) + 1) want += x[i] * y[k - i];
      CHECK(fabs(c[k] - want) < 1e-9);
    }
  }
  rep(lg, 0, 8) {
    int n = 1 << lg;
    vector<cd> a(n);
    for (cd& v : a) v = cd(rnd_real(-1, 1), rnd_real(-1, 1));
    vector<cd> f = a, want = brute_dft(a);
    fft(f);
    rep(k, 0, n) CHECK(abs(f[k] - want[k]) < 1e-9);
    fft(f), reverse(f.begin() + 1, f.end());
    rep(k, 0, n) CHECK(abs(f[k] / double(n) - a[k]) < 1e-9);
  }
  // Precision bound from the header: 2^19 values up to 3e4.
  check_constant(1 << 19, 1 << 19, 30'000, 30'000);
  check_alternating(1 << 19, 1 << 19, 30'000, 30'000);
  check_constant(1, 1 << 19, 30'000, 30'000);
  check_constant(1, 1, 30'000, 30'000);
  check_constant(1 << 10, 1, 1'000'000'000, 0); // exact zeros
  // Unbalanced inputs inside the bound |a| |b| < 5e14.
  check_constant(1 << 17, 1 << 17, 1'000'000'000, 1);
  check_constant(1 << 17, 1 << 17, 1, 1'000'000'000);
  check_constant(1 << 18, 1 << 10, 141'000, 141'000);
  check_constant(1 << 10, 1 << 18, 141'000, 141'000);
  check_alternating(1 << 18, 1 << 10, 141'000, 141'000);
  int n = 100'000;
  rep(t, 0, 2) {
    vector<ll> a = rnd_vec(n, 0, t ? 30'000 : 1'000'000);
    vector<ll> b = rnd_vec(n, 0, t ? 30'000 : 1);
    vd c = fft_conv(vd(all(a)), vd(all(b)));
    rep(s, 0, 30) {
      int k = (int)rnd(0, 2 * n - 2);
      ll want = 0;
      rep(i, max(0, k - n + 1), min(k, n - 1) + 1) want += a[i] * b[k - i];
      CHECK_EQ(llround(c[k]), want);
    }
  }
  puts("OK poly/fft");
}
