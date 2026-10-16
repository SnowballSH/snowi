/**
 * Title: Power series operations
 * Description: Formal power series modulo $x^n$ and `MOD`:
 *   `poly_inv` gives $1 / a$, `poly_log` gives $ln a$, `poly_exp`
 *   gives $exp a$ and `poly_pow` gives $a^k$, each as `n`
 *   coefficients. `poly_divmod(a, b)` divides polynomials:
 *   $a = q b + r$ with `sz(q) = sz(a) - sz(b) + 1` (empty when that
 *   is not positive) and `sz(r) = sz(b) - 1`.
 * Precondition: Entries in `[0, MOD)`, `n` $>= 0$. `poly_inv` needs
 *   $a_0 != 0$, `poly_log` needs $a_0 = 1$, `poly_exp` needs
 *   $a_0 = 0$; `poly_pow` takes any `a` and $k >= 0$. `poly_divmod`
 *   needs `b.back()` $!= 0$. Lengths up to $2^22$ (NTT size).
 * Time: $O(n log n)$. For $n = 5 dot 10^5$: inverse 0.3 s, log 0.55
 *   s, exp 0.65 s, pow 1.2 s (`conv` of that length: 0.25 s).
 *   Newton steps and NTT sizes round up to powers of two: $n = 2^18$
 *   takes half as long, $n = 2^18 + 1$ as long as $5 dot 10^5$.
 * Usage:
 *   vl b = poly_inv(a, n); // a b = 1 mod x^n
 *   vl e = poly_exp(poly_log(a, n), n); // e = a mod x^n
 *   auto [q, r] = poly_divmod(a, b);
 * Tweak: Square root: `poly_pow(a, (MOD + 1) / 2, n)` if $a_0 = 1$;
 *   else strip $x^v$ and $c = a_v$ as `poly_pow` does and multiply
 *   that root by $x^(v \/ 2)$ and `mod_sqrt(c, MOD)`
 *   (`numtheory/mod-sqrt`): none if $v$ is odd or $c$ is not a
 *   square.
 * Example:
 *   // Partitions p(0..n-1) = 1 / prod_k (1 - x^k), via
 *   // ln(1 / (1 - x^k)) = sum_j x^(jk) / j:
 *   vl L(n), inv(n, 1); // inv[j] = 1 / j
 *   rep(i, 2, n) inv[i] = (MOD - MOD / i) * inv[MOD % i] % MOD;
 *   rep(k, 1, n) for (int j = 1; j * k < n; j++)
 *     L[j * k] = (L[j * k] + inv[j]) % MOD;
 *   vl p = poly_exp(L, n); // 1 1 2 3 5 7 11 15 ...
 *   // Catalan numbers (1 - sqrt(1 - 4x)) / (2x):
 *   vl s = poly_pow({1, MOD - 4}, (MOD + 1) / 2, n + 1);
 *   vl cat(n);
 *   rep(i, 0, n) cat[i] = (MOD - s[i + 1]) * modinv(2) % MOD;
 * Notes: Newton doubles the number of correct terms per step:
 *   $b <- b (2 - a b)$ for the inverse; exp keeps $g = 1 / f$
 *   alongside and sets $f <- f (1 + a - ln f)$. The products are
 *   cyclic of length $2m$ (`cyc_mul`): the wrapped terms only land
 *   where the result is already known. $ln a = integral a' / a$.
 *   `poly_pow` strips the leading $x^v$ and coefficient $c$, then
 *   takes $exp(k ln a)$; the result is 0 when $v k >= n$. Reduce $k$
 *   mod `MOD` only inside the log, never in $c^k$ or in $v k$. Cues:
 *   labelled sets of connected pieces (exp), connected from all
 *   (log), sequences $1 / (1 - A)$ (inverse), $A^k$ for huge $k$,
 *   "remainder of a polynomial division" (divmod). Build inputs
 *   in $O(n log n)$ too: a `modinv` per term of a harmonic double
 *   loop costs more than the `poly_exp`; use the linear inverse
 *   table of the Example. Below $n approx 3000$ the $O(n^2)$
 *   recurrences in the generating functions note are shorter to
 *   write.
 * Source: cp-algorithms, Operations on polynomials and series; exp
 *   with the inverse kept alongside after Bernstein, Removing
 *   redundancy in high-precision Newton iteration.
 * Status: stress-tested against $O(n^2)$ series arithmetic
 */
#pragma once
#include "ntt.hpp"

// a b mod x^n - 1 for n = sz(a), given B = ntt(b).
vl cyc_mul(vl a, const vl& B) {
  ll k = modinv(sz(a));
  ntt(a);
  rep(i, 0, sz(a)) a[i] = a[i] * B[i] % MOD * k % MOD;
  ntt(a), reverse(a.begin() + 1, a.end());
  return a;
}
vl poly_inv_step(const vl& a, vl b) { // mod x^k -> x^2k
  int k = sz(b);
  vl f(2 * k), g(b);
  copy(a.begin(), a.begin() + min(sz(a), 2 * k), f.begin());
  g.resize(2 * k), ntt(g), f = cyc_mul(f, g);
  fill(f.begin(), f.begin() + k, 0); // a b - 1
  f = cyc_mul(f, g), b.resize(2 * k);
  rep(i, k, 2 * k) b[i] = (MOD - f[i]) % MOD;
  return b;
}
vl poly_inv(const vl& a, int n) {
  vl b{modinv(a[0])};
  while (sz(b) < n) b = poly_inv_step(a, b);
  b.resize(n);
  return b;
}
vl poly_log(const vl& a, int n) {
  if (n <= 1) return vl(n);
  vl d(n - 1), inv(n, 1), res(n);
  rep(i, 1, min(sz(a), n)) d[i - 1] = a[i] * i % MOD;
  rep(i, 2, n) inv[i] = (MOD - MOD / i) * inv[MOD % i] % MOD;
  d = conv(d, poly_inv(a, n - 1));
  rep(i, 1, n) res[i] = d[i - 1] * inv[i] % MOD;
  return res;
}
vl poly_exp(const vl& a, int n) {
  vl f{1}, g{1}, inv(2 * n + 1, 1); // f = exp a, g = 1 / f
  rep(i, 2, 2 * n + 1)
    inv[i] = (MOD - MOD / i) * inv[MOD % i] % MOD;
  for (int m = 1; m < n; m *= 2) {
    if (sz(g) < m) g = poly_inv_step(f, g);
    vl F(f), G(g), P(2 * m), T(2 * m);
    F.resize(2 * m), G.resize(2 * m), ntt(F), ntt(G);
    rep(i, 0, min(2 * m, sz(a) - 1))
      P[i] = a[i + 1] * (i + 1) % MOD;
    P = cyc_mul(P, F);
    fill(P.begin(), P.begin() + m - 1, 0); // f a' - f'
    P = cyc_mul(P, G);                     // a' - (ln f)'
    rep(i, m, 2 * m) T[i] = P[i - 1] * inv[i] % MOD;
    T = cyc_mul(T, F), f.resize(2 * m);
    rep(i, m, 2 * m) f[i] = T[i]; // f + f (a - ln f)
  }
  f.resize(n);
  return f;
}
vl poly_pow(const vl& a, ll k, int n) {
  vl res(n);
  if (!k && n) res[0] = 1;
  int v = 0;
  while (v < sz(a) && !a[v]) v++;
  if (!k || v == sz(a) || v * min(k, (ll)n) >= n) return res;
  int m = n - int(v * k);
  ll c = modinv(a[v]), ck = modpow(a[v], k);
  vl b(m);
  rep(i, 0, min(m, sz(a) - v)) b[i] = a[v + i] * c % MOD;
  b = poly_log(b, m);
  for (ll& x : b) x = x * (k % MOD) % MOD;
  b = poly_exp(b, m);
  rep(i, 0, m) res[n - m + i] = b[i] * ck % MOD;
  return res;
}
pair<vl, vl> poly_divmod(vl a, const vl& b) {
  int n = sz(a), m = sz(b), k = n - m + 1;
  vl q;
  if (k > 0) {
    vl ra(a.rbegin(), a.rbegin() + k), rb(b.rbegin(), b.rend());
    q = conv(ra, poly_inv(rb, k)), q.resize(k), reverse(all(q));
  }
  vl bq = conv(b, q);
  a.resize(m - 1);
  rep(i, 0, m - 1)
    a[i] = (a[i] - (i < sz(bq) ? bq[i] : 0) + MOD) % MOD;
  return {q, a};
}
