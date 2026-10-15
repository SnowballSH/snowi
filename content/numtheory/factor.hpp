/**
 * Title: Primality and factoring
 * Description: `is_prime(n)`: deterministic Miller-Rabin for all
 *   64-bit `n`. `factor(n)`: the prime factors of `n` in increasing
 *   order, with multiplicity, by Pollard-Brent rho.
 * Precondition: $1 <= n < 2^63$; `factor(1)` is empty.
 * Time: `is_prime`: 7 modular powers. `factor`: $O(n^(1/4))$
 *   expected: 0.07 ms on average for a random $n <= 10^18$, 3 ms in
 *   the worst case (two primes near $3 dot 10^9$).
 * Usage:
 *   bool p = is_prime(1'000'000'007);
 *   vector<ll> f = factor(360); // {2, 2, 2, 3, 3, 5}
 * Tweak: `mulmod` and `powmod` take nonnegative arguments and any
 *   modulus $1 <= m < 2^63$. For $n <= 10^12$ trial division up to
 *   $sqrt(n)$ is simpler.
 * Example:
 *   vector<ll> divs = {1}; // all divisors of n, unsorted
 *   ll last = 0;
 *   for (int from = 0, to; ll p : factor(n)) {
 *     if (p != last) from = 0, last = p;
 *     to = sz(divs);
 *     rep(i, from, to) divs.push_back(divs[i] * p);
 *     from = to;
 *   }
 * Notes: Factor once, then derive what you need: the divisors (at
 *   most 103680 for $n <= 10^18$), $phi(n) = n product_p (1 - 1/p)$,
 *   or the multiplicative order of `a` modulo a prime `p`: start from
 *   $o = p - 1$ and divide out each prime factor `q` while
 *   $a^(o/q) = 1$. Rho splits off the smallest prime factor `r` in
 *   about $sqrt(r)$ steps, so the worst case is the largest balanced
 *   semiprime the input allows: budget for it, not for a random `n`.
 * Source: Miller-Rabin bases by Jim Sinclair; Brent (1980); kactl
 *   MillerRabin and Factor (CC0).
 * Status: stress-tested against trial division and random-base
 *   Miller-Rabin
 */
#pragma once

ll mulmod(ll a, ll b, ll m) { return (__int128)a * b % m; }
ll powmod(ll b, ll e, ll m) {
  ll r = 1 % m;
  for (b %= m; e; e >>= 1, b = mulmod(b, b, m))
    if (e & 1) r = mulmod(r, b, m);
  return r;
}
bool is_prime(ll n) {
  if (n < 2 || n % 2 == 0) return n == 2;
  ll A[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
  int s = __builtin_ctzll(n - 1);
  for (ll a : A) {
    ll x = powmod(a, (n - 1) >> s, n);
    if (a % n == 0 || x == 1) continue;
    rep(i, 1, s)
      if (x != n - 1) x = mulmod(x, x, n);
    if (x != n - 1) return false;
  }
  return true;
}
ll rho(ll n) { // a proper divisor of an odd composite n
  for (ll c = 1;; c++) {
    auto f = [&](ll v) { return ((__int128)v * v + c) % n; };
    ll x, y = 2, g = 1, q = 1, ys = 2;
    for (ll r = 1; g == 1; r *= 2) {
      x = y;
      rep(i, 0, r) y = f(y);
      for (ll k = 0; k < r && g == 1; k += 128) {
        ys = y;
        for (ll i = 0; i < 128 && i < r - k; i++)
          y = f(y), q = mulmod(q, abs(x - y), n);
        g = gcd(q, n);
      }
    }
    if (g == n) // redo the last batch step by step
      for (g = 1; g == 1;) ys = f(ys), g = gcd(abs(x - ys), n);
    if (g != n) return g;
  }
}
vector<ll> factor(ll n) {
  if (n == 1) return {};
  if (is_prime(n)) return {n};
  ll d = n % 2 ? rho(n) : 2;
  vector<ll> a = factor(d), b = factor(n / d);
  a.insert(a.end(), all(b));
  inplace_merge(a.begin(), a.end() - sz(b), a.end());
  return a;
}
