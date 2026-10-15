#include "../testlib.hpp"
#include "../../content/numtheory/euclid.hpp"

typedef __int128 L;

ll brute_inv(ll a, ll m) {
  rep(x, 0, m)
    if ((L)a * x % m == 1 % m) return x;
  return -1;
}

// Smallest x >= 0 meeting all congruences, and the lcm, by search.
pair<ll, ll> brute_crt(const vector<ll>& r, const vector<ll>& m) {
  ll l = 1;
  for (ll v : m) l = lcm(l, v);
  rep(x, 0, l) {
    bool ok = true;
    rep(i, 0, sz(r)) ok &= x % m[i] == r[i];
    if (ok) return {x, l};
  }
  return {-1, -1};
}

void check_euclid(ll a, ll b) {
  ll x, y, g = euclid(a, b, x, y);
  CHECK_EQ(g, gcd(a, b));
  CHECK((L)a * x + (L)b * y == g);
  if (a && b) CHECK(abs(x) <= b / g && abs(y) <= a / g);
}

// Residues reduced modulo m and n, lcm(m, n) < 2^63.
void check_crt_big(ll m, ll n) {
  ll g = gcd(m, n);
  if ((L)(m / g) * n > LLONG_MAX) return;
  ll a = rnd(0, m - 1), b = rnd(0, n - 1);
  if (rnd(0, 1)) b = (ll)(((L)a + (L)g * rnd(0, n)) % n);
  auto [x, l] = crt(a, m, b, n);
  if ((b - a) % g) {
    CHECK(x == -1 && l == -1);
    return;
  }
  CHECK_EQ(l, m / g * n);
  CHECK(0 <= x && x < l && x % m == a && x % n == b);
}

int main() {
  elapsed();
  ll acc = 0;
  rep(i, 0, 50000) {
    ll m = rnd(1, 1LL << 31), n = rnd(1, 1LL << 31);
    acc ^= crt(rnd(0, m - 1), m, rnd(0, n - 1), n).second;
    acc ^= inv_mod(rnd(0, LLONG_MAX), LLONG_MAX);
  }
  CHECK(acc != 42);
  fprintf(stderr, "  5e4 crt + inv_mod: %.2fs\n", elapsed());
  rep(a, 0, 300)
    rep(b, 0, 300) check_euclid(a, b);
  rep(it, 0, iters(10000)) {
    check_euclid(rnd(0, LLONG_MAX), rnd(0, LLONG_MAX));
    ll g = rnd(1, 1000);
    check_euclid(g * rnd(0, LLONG_MAX / g), g * rnd(0, LLONG_MAX / g));
  }
  check_euclid(0, 0), check_euclid(LLONG_MAX, LLONG_MAX - 1);
  rep(m, 1, 60)
    rep(a, -m, 2 * m) CHECK_EQ(inv_mod(a, m), brute_inv((a % m + m) % m, m));
  rep(it, 0, iters(10000)) {
    ll m = rnd(1, LLONG_MAX), a = rnd(LLONG_MIN, LLONG_MAX);
    ll v = inv_mod(a, m), r = a % m < 0 ? a % m + m : a % m;
    if (gcd(r, m) != 1) CHECK_EQ(v, -1LL);
    else CHECK(0 <= v && v < m && (L)r * v % m == 1 % m);
  }
  rep(m, 1, 16)
    rep(n, 1, 16)
      rep(a, 0, m)
        rep(b, 0, n) {
          auto [x, l] = crt(a, m, b, n);
          CHECK(make_pair(x, l) == brute_crt({a, b}, {m, n}));
        }
  rep(it, 0, iters(3000)) {
    int k = (int)rnd(0, 4);
    vector<ll> r(k), m(k);
    rep(i, 0, k) m[i] = rnd(1, 10), r[i] = rnd(0, m[i] - 1);
    CHECK(crt_all(r, m) == brute_crt(r, m));
  }
  rep(it, 0, iters(10000)) {
    check_crt_big(rnd(1, LLONG_MAX), rnd(1, LLONG_MAX));
    ll g = rnd(1, 1LL << rnd(0, 40));
    check_crt_big(g * rnd(1, (1LL << 22)), g * rnd(1, (1LL << 22)));
  }
  ll p = (1LL << 61) - 1;
  CHECK(crt(p - 1, p, 2, 4) == make_pair(p - 1, 4 * p));
  CHECK(crt(0, p, 1, 4) == make_pair(3 * p, 4 * p));
  CHECK(crt(5, 6, 1, 4) == make_pair(5LL, 12LL));
  CHECK(inv_mod(3, 10) == 7 && crt(2, 4, 4, 6) == make_pair(10LL, 12LL));
  CHECK(crt_all({1, 2, 3}, {2, 3, 5}) == make_pair(23LL, 30LL));
  CHECK(crt(0, 6, 1, 4) == make_pair(-1LL, -1LL));
  CHECK(crt_all({2, 3}, {4, 6}) == make_pair(-1LL, -1LL));

  // Example: all solutions of a x + b y = c.
  rep(it, 0, iters(5000)) {
    ll a = rnd(-50, 50), b = rnd(-50, 50), c = rnd(-200, 200);
    if (!a && !b) continue;
    ll x, y, g = euclid(abs(a), abs(b), x, y);
    bool brute = false;
    rep(u, -200, 201) brute |= b ? (c - a * u) % b == 0 : a * u == c;
    CHECK_EQ(c % g == 0, brute);
    if (c % g) continue;
    x *= c / g * (a < 0 ? -1 : 1);
    y *= c / g * (b < 0 ? -1 : 1);
    ll k = rnd(-5, 5);
    CHECK_EQ(a * (x + k * (b / g)) + b * (y - k * (a / g)), c);
  }

  puts("OK numtheory/euclid");
}
