/**
 * Title: Lagrange interpolation
 * Description: `lagrange(y, x)` returns $f(x)$, where $f$ is the
 *   polynomial of degree below `sz(y)` with $f(i) = y_i$ for every
 *   index $i$. `interpolate(xs, ys)` returns the coefficients of the
 *   polynomial of degree below $n$ through $n$ points $(x_i, y_i)$,
 *   lowest degree first.
 * Precondition: Arithmetic modulo the prime `MOD`; values in
 *   `[0, MOD)`. `lagrange`: `y` not empty, `sz(y)` $<=$ `MOD`, any
 *   `x` (reduced inside, so $x >= "MOD"$ is fine). `interpolate`:
 *   the `xs` are distinct and in `[0, MOD)`.
 * Time: $O(n)$ for `lagrange`, 0.06 s for $n = 10^6$; $O(n^2)$ for
 *   `interpolate`, 0.5 s for $n = 5000$.
 * Usage:
 *   ll s = lagrange({0, 1, 5, 14}, 10); // sum_(i <= 10) i^2: 385
 *   vl c = interpolate(xs, ys); // f(t) = sum_k c[k] t^k
 * Tweak: Nodes $a, a + h, a + 2h, dots$: evaluate at $(x - a) / h$
 *   (division mod `MOD`). Many queries at arbitrary nodes:
 *   `interpolate` once, then Horner. An exact integer answer without
 *   a modulus: interpolate modulo two or three primes and combine by
 *   CRT; never in floating point.
 * Example:
 *   // sum_(i = 1)^n i^k for n up to 1e18: degree k + 1, so k + 2
 *   // samples at 0, 1, ..., k + 1.
 *   vl y(k + 2);
 *   rep(i, 1, k + 2) y[i] = (y[i - 1] + modpow(i, k)) % MOD;
 *   ll s = lagrange(y, n);
 * Notes: The answer is a polynomial in a huge parameter $K$ when a
 *   DP over values $0, dots, K$ only adds, multiplies and takes
 *   prefix sums: a prefix sum raises the degree by one, a product
 *   adds degrees. Compute the answer for $K = 0, dots, d$ by brute
 *   force or the slow DP and call `lagrange`. Take two extra samples
 *   and check that they are reproduced: a wrong degree bound fails
 *   silently. With $min(K, a_i)$ or thresholds the answer is only
 *   piecewise polynomial; sample beyond the last breakpoint (or per
 *   piece), and above any small $K$ the statement treats specially.
 *   Cues: "$K$ up to $10^9$" next to $n <= 3000$, "values in
 *   $[1, K]$", "$sum i^k$", child $<=$ parent labellings of a tree
 *   (degree = size). A DP dimension that only adds under merges can
 *   be replaced by a variable $z$: run the DP at $n + 1$ numbers $z$
 *   and interpolate.
 * Source: Classical Lagrange formula with prefix and suffix
 *   products.
 * Status: stress-tested against Horner evaluation of random
 *   polynomials
 */
#pragma once
#include "../numtheory/modpow.hpp"

typedef vector<ll> vl;
ll lagrange(const vl& y, ll x) {
  int d = sz(y) - 1;
  x = (x % MOD + MOD) % MOD;
  vl suf(d + 2, 1), ifact(d + 1, 1);
  for (int i = d; i >= 0; i--)
    suf[i] = suf[i + 1] * ((x - i + MOD) % MOD) % MOD;
  rep(i, 1, d + 1) ifact[d] = ifact[d] * i % MOD;
  ifact[d] = modinv(ifact[d]);
  for (int i = d; i > 0; i--) ifact[i - 1] = ifact[i] * i % MOD;
  ll res = 0, pre = 1; // pre = prod_(j < i) (x - j)
  rep(i, 0, d + 1) {
    ll t = y[i] * pre % MOD * suf[i + 1] % MOD * ifact[i] % MOD;
    t = t * ifact[d - i] % MOD;
    res = (res + ((d - i) % 2 ? MOD - t : t)) % MOD;
    pre = pre * ((x - i + MOD) % MOD) % MOD;
  }
  return res;
}
vl interpolate(const vl& x, const vl& y) {
  int n = sz(x);
  vl P(n + 1), Q(n), res(n); // P = prod_j (z - x_j)
  P[0] = 1;
  rep(j, 0, n)
    for (int k = j + 1; k >= 0; k--)
      P[k] = ((k ? P[k - 1] : 0) + (MOD - x[j]) * P[k]) % MOD;
  rep(i, 0, n) {
    ll c = 0, den = 0; // Q = P / (z - x_i), den = Q(x_i)
    for (int k = n; k > 0; k--) {
      Q[k - 1] = c = (P[k] + c * x[i]) % MOD;
      den = (den * x[i] + c) % MOD;
    }
    c = y[i] * modinv(den) % MOD;
    rep(k, 0, n) res[k] = (res[k] + c * Q[k]) % MOD;
  }
  return res;
}
