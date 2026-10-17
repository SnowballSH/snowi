/**
 * Title: Number theoretic transform
 * Description: `conv(a, b)` multiplies polynomials modulo `MOD`.
 *   `ntt(a)` replaces `a` by its values at $w^0, w^1, dots$ for
 *   $w = 3^(("MOD" - 1) / n)$, a primitive `n`-th root of unity, in
 *   natural order.
 * Precondition: `MOD` $= 998244353 = 119 dot 2^23 + 1$ with primitive
 *   root 3. Entries in `[0, MOD)`. `sz(a)` is a power of two, at most
 *   $2^23$, so `conv` needs `sz(a) + sz(b) - 1` $<= 2^23$.
 * Time: $O(n log n)$. A $2^19 times 2^19$ product takes about 0.25 s.
 * Usage:
 *   vl c = conv(a, b); // sz(a) + sz(b) - 1 entries
 * Tweak: Inverse transform: `ntt(a)`, then
 *   `reverse(a.begin() + 1, a.end())` and multiply by `modinv(n)`.
 *   Other primes: the table in Polynomial techniques.
 * Notes: `conv` falls back to the schoolbook product up to 40
 *   entries, which matters for divide-and-conquer products of many
 *   small polynomials. Cues: "for every $k$" with a sum over
 *   $i + j = k$, a product of many small polynomials, the modulus
 *   998244353 itself. For a correlation $sum_i a_i b_(i + d)$ reverse
 *   `a`; the answer for `d` is `c[sz(a) - 1 + d]`. Reduce negative
 *   values first.
 * Source: kactl NumberTheoreticTransform (CC0).
 * Status: stress-tested against the schoolbook product
 */
#pragma once
#include "../numtheory/modpow.hpp"

typedef vector<ll> vl;
void ntt(vl& a) {
  int n = sz(a), L = __lg(n);
  static vl rt(2, 1);
  for (static int k = 2, s = 2; k < n; k *= 2, s++) {
    rt.resize(n);
    ll z[] = {1, modpow(3, MOD >> s)};
    rep(i, k, 2 * k) rt[i] = rt[i / 2] * z[i & 1] % MOD;
  }
  vi rev(n);
  rep(i, 0, n) rev[i] = (rev[i / 2] | (i & 1) << L) / 2;
  rep(i, 0, n)
    if (i < rev[i]) swap(a[i], a[rev[i]]);
  for (int k = 1; k < n; k *= 2)
    for (int i = 0; i < n; i += 2 * k)
      rep(j, 0, k) {
        ll z = rt[j + k] * a[i + j + k] % MOD, &ai = a[i + j];
        a[i + j + k] = ai - z + (z > ai ? MOD : 0);
        ai += (ai + z >= MOD ? z - MOD : z);
      }
}
vl conv(const vl& a, const vl& b) {
  if (a.empty() || b.empty()) return {};
  int s = sz(a) + sz(b) - 1, n = bit_ceil(unsigned(s));
  vl L(a), R(b), out(n);
  if (min(sz(a), sz(b)) <= 40) {
    rep(i, 0, sz(a))
      rep(j, 0, sz(b))
        out[i + j] = (out[i + j] + a[i] * b[j]) % MOD;
  } else {
    L.resize(n), R.resize(n), ntt(L), ntt(R);
    ll inv = modinv(n);
    rep(i, 0, n)
      out[-i & (n - 1)] = L[i] * R[i] % MOD * inv % MOD;
    ntt(out);
  }
  return {out.begin(), out.begin() + s};
}
