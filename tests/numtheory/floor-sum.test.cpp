#include "../testlib.hpp"
#include "../../content/numtheory/floor-sum.hpp"

typedef __int128 L;

L floor_div(L x, L m) { return x / m - (x % m < 0); }

L brute(ll n, ll m, ll a, ll b) {
  L s = 0;
  rep(i, 0, n) s += floor_div((L)a * i + b, m);
  return s;
}

int main() {
  rep(it, 0, iters(30000)) {
    ll n = rnd(0, 30), m = rnd(1, 30);
    ll a = rnd(-60, 60), b = rnd(-60, 60);
    CHECK((L)floor_sum(n, m, a, b) == brute(n, m, a, b));
  }
  rep(it, 0, iters(3000)) {
    ll n = rnd(0, 2000), m = rnd(1, 1LL << rnd(0, 31));
    ll a = rnd(-(1LL << 31), 1LL << 31), b = rnd(-(1LL << 40), 1LL << 40);
    CHECK((L)floor_sum(n, m, a, b) == brute(n, m, a, b));
  }
  // Extremes, n and m near 2^32: if gcd(a, m) = 1 then
  // floor_sum(m, m, a, b) = (a - 1)(m - 1) / 2 + b exactly.
  ll M = (1LL << 32) - 5;
  CHECK_EQ(floor_sum(M, M, M - 2, 0), (M - 3) / 2 * (M - 1));
  CHECK_EQ(floor_sum(M, M, 2 - M, 7), -(M - 1) / 2 * (M - 1) + 7);
  rep(it, 0, iters(2000)) {
    ll m = rnd(1, (1LL << 32) - 1), a = rnd(-(1LL << 31), 1LL << 31);
    if (gcd(a, m) != 1) continue;
    ll b = rnd(-(1LL << 40), 1LL << 40);
    CHECK_EQ(floor_sum(m, m, a, b), (a - 1) * (m - 1) / 2 + b);
  }
  CHECK_EQ(floor_sum(M + 4, M, 0, LLONG_MAX / (M + 4) * M),
           LLONG_MAX / (M + 4) * (M + 4));
  rep(it, 0, iters(1000)) {
    ll n = rnd(0, (1LL << 32) - 1), m = rnd(1, (1LL << 32) - 1);
    ll a = rnd(-(1LL << 20), 1LL << 20), b = rnd(-m, m);
    // sum_i (a i + b) = m * floor_sum + sum_i ((a i + b) mod m)
    L lin = (L)a * n * (n - 1) / 2 + (L)b * n;
    L mods = lin - (L)m * floor_sum(n, m, a, b);
    CHECK(0 <= mods && mods <= (L)(m - 1) * n);
  }
  CHECK_EQ(floor_sum(0, 5, 3, 2), 0LL);
  CHECK_EQ(floor_sum(1, 1, -7, -9), -9LL);

  // Example formulas.
  rep(it, 0, iters(3000)) {
    ll n = rnd(0, 20), m = rnd(1, 20), a = rnd(0, 20), b = rnd(0, 20);
    ll c = rnd(0, m), pts = 0, rem = 0, cnt = 0, F = floor_sum(n, m, a, b);
    rep(x, 0, n) {
      rep(y, 1, 401) pts += y * m <= a * x + b;
      rem += (a * x + b) % m, cnt += (a * x + b) % m < c;
    }
    CHECK_EQ(F, pts);
    CHECK_EQ(a * (n * (n - 1) / 2) + b * n - m * F, rem);
    CHECK_EQ(F - floor_sum(n, m, a, b - c), cnt);
  }

  elapsed();
  ll acc = 0;
  rep(it, 0, 100000) {
    ll m = rnd(1, (1LL << 32) - 1), n = rnd(0, (1LL << 32) - 1);
    acc ^= floor_sum(n, m, rnd(-m, m), rnd(-m, m));
  }
  CHECK(acc != 42);
  fprintf(stderr, "  1e5 calls, n, m < 2^32: %.2fs\n", elapsed());
  puts("OK numtheory/floor-sum");
}
