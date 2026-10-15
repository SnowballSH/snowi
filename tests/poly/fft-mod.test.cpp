#include "../testlib.hpp"
#include "../../content/poly/fft-mod.hpp"

vl brute_conv(const vl& a, const vl& b, ll m) {
  if (a.empty() || b.empty()) return {};
  vl c(sz(a) + sz(b) - 1);
  rep(i, 0, sz(a))
    rep(j, 0, sz(b)) c[i + j] = (c[i + j] + a[i] * b[j]) % m;
  return c;
}

// All entries v: c[k] = min(k + 1, len - k) v^2 mod m. The worst case
// for rounding is v near m / 2; v = m - 1 needs the balanced split.
void check_constant(int n, ll m, ll v) {
  vl c = conv_mod(vl(n, v), vl(n, v), m);
  ll sq = v * v % m;
  CHECK_EQ(sz(c), 2 * n - 1);
  rep(k, 0, 2 * n - 1) CHECK_EQ(c[k], min(k + 1, 2 * n - 1 - k) * sq % m);
}

int main() {
  const vl mods = {INT_MAX,     (1 << 30) - 1, 1 << 30, 1'000'000'007,
                   998'244'353, 1009,          3,       2};
  rep(it, 0, iters(1500)) {
    ll m = it % 3 ? mods[rnd(0, sz(mods) - 1)] : rnd(2, (1LL << 31) - 1);
    int n = (int)rnd(0, it % 3 ? 40 : 120);
    int k = (int)rnd(0, it % 3 ? 40 : 120);
    vl a = rnd_vec(n, 0, m - 1), b = rnd_vec(k, 0, m - 1);
    ll v[] = {0, m - 1, m / 2, m / 2 + 1};
    if (it % 5 == 0) fill(all(a), v[it % 4]), fill(all(b), v[it / 5 % 4]);
    CHECK(conv_mod(a, b, m) == brute_conv(a, b, m));
  }
  // Limits from the header: length 2^20 below 2^30, 2^19 below 2^31.
  check_constant(1 << 19, (1 << 30) - 1, (1 << 29) - 1);
  check_constant(1 << 18, (1LL << 31) - 1, 1 << 30);
  check_constant(1 << 18, (1LL << 31) - 1, (1LL << 31) - 2);
  int n = 1 << 16;
  ll m = 1'000'000'007;
  vl a = rnd_vec(n, 0, m - 1), b = rnd_vec(n, 0, m - 1);
  vl c = conv_mod(a, b, m);
  rep(t, 0, 20) {
    int k = (int)rnd(0, 2 * n - 2);
    ll want = 0;
    rep(i, max(0, k - n + 1), min(k, n - 1) + 1)
      want = (want + a[i] * b[k - i]) % m;
    CHECK_EQ(c[k], want);
  }
  puts("OK poly/fft-mod");
}
