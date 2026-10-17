/**
 * Title: Extended Euclid and CRT
 * Description: `euclid(a, b, x, y)` returns $g = gcd(a, b)$ and sets
 *   `x`, `y` with $a x + b y = g$; if $a, b > 0$ then $|x| <= b/g$
 *   and $|y| <= a/g$. `inv_mod(a, m)` is the inverse of `a` modulo
 *   any `m`, or -1 if $gcd(a, m) != 1$. `crt(a, m, b, n)` solves
 *   $x equiv a thick (mod m)$, $x equiv b thick (mod n)$ for any
 *   moduli, coprime or not: it returns $(x, op("lcm")(m, n))$ with
 *   $0 <= x < op("lcm")$, or $(-1, -1)$ if there is no solution.
 *   `crt_all` folds a list of congruences the same way.
 * Precondition: `euclid`: $a, b >= 0$. Moduli are positive. `crt`
 *   needs reduced residues ($0 <= a < m$, $0 <= b < n$) and the lcm
 *   to fit in `ll`.
 * Time: $O(log min(a, b))$ per call.
 * Usage:
 *   ll x, y, g = euclid(a, b, x, y);
 *   ll inv = inv_mod(3, 10); // 7
 *   auto [r, l] = crt(2, 4, 4, 6); // (10, 12)
 *   auto [r2, l2] = crt_all({1, 2, 3}, {2, 3, 5}); // (23, 30)
 * Example:
 *   // a x + b y = c, any signs, a or b nonzero: solvable iff g | c
 *   ll x, y, g = euclid(abs(a), abs(b), x, y);
 *   if (c % g == 0) {
 *     x *= c / g * (a < 0 ? -1 : 1);
 *     y *= c / g * (b < 0 ? -1 : 1);
 *   } // all solutions: x + k * (b / g), y - k * (a / g)
 * Notes: A linear congruence $a x equiv c thick (mod m)$ is solvable
 *   iff $g = gcd(a, m)$ divides `c`, with `g` solutions modulo `m`:
 *   divide everything by `g` and multiply by `inv_mod(a/g, m/g)`.
 *   Use CRT when the modulus is composite (solve modulo each prime
 *   power, then combine) or when a quantity is known modulo several
 *   periods. Multiplying `x` by `c / g` above can overflow `ll`; use
 *   `__int128` or reduce `x` modulo $b/g$ first.
 * Source: kactl Euclid and CRT (CC0); AtCoder Library inv_mod/crt
 *   (CC0).
 * Status: stress-tested against brute-force search
 */
#pragma once

ll euclid(ll a, ll b, ll& x, ll& y) {
  if (!b) return x = 1, y = 0, a;
  ll g = euclid(b, a % b, y, x);
  return y -= a / b * x, g;
}
ll inv_mod(ll a, ll m) {
  ll x, y, g = euclid((a %= m) < 0 ? a + m : a, m, x, y);
  return g == 1 ? (x < 0 ? x + m : x) : -1;
}
pair<ll, ll> crt(ll a, ll m, ll b, ll n) {
  ll x, y, g = euclid(m, n, x, y), k = n / g;
  if ((b - a) % g) return {-1, -1};
  ll t = (ll)((__int128)((b - a) / g) * x % k);
  return {a + (t < 0 ? t + k : t) * m, m * k};
}
pair<ll, ll> crt_all(const vector<ll>& r, const vector<ll>& m) {
  pair<ll, ll> res{0, 1};
  rep(i, 0, sz(r)) {
    res = crt(res.first, res.second, r[i], m[i]);
    if (res.second < 0) break;
  }
  return res;
}
