/**
 * Title: Linear recurrences
 * Description: `berlekamp_massey(s)` returns the shortest `c` with
 *   $s_i = sum_(j = 1)^L c_j s_(i - j)$ for every $i >= L$, as
 *   `c[j - 1]` $= c_j$. `linear_rec(s, c, k)` returns $s_k$ from the
 *   first $L$ terms, via $x^k$ modulo the characteristic polynomial.
 * Precondition: Arithmetic modulo the prime `MOD`. Entries of `s`
 *   and `c` are below `MOD` in absolute value, so a hand-written
 *   `c = {2, -1}` is fine. A recurrence of order $L$ is found from
 *   $2L$ terms; give a margin, as fewer terms give a shorter wrong
 *   answer. $k >= 0$.
 * Time: $O(n^2)$ for $n$ terms; $O(L^2 log k)$ for `linear_rec`.
 *   About 0.4 s for $L = 1000$ and $k = 10^18$, 0.3 s for $n = 10^4$.
 * Usage:
 *   vl c = berlekamp_massey({1, 1, 3, 5, 11, 21}); // {1, 2}
 *   ll x = linear_rec({1, 1}, c, k); // s_k, 0-indexed
 * Tweak: Squarefree composite modulus: make `MOD` a variable, run
 *   modulo each prime factor and combine with CRT. A modulus with a
 *   square factor ($2^32$, $10^9$) needs a matrix power instead. For
 *   $L$ in the thousands with huge $k$ replace `mul` by NTT
 *   multiplication and division (or use Bostan-Mori).
 * Example:
 *   // Walks of length k in a graph with V vertices satisfy a
 *   // recurrence of order <= V: get 2V + 10 terms by DP, then
 *   vl s = first_terms(2 * V + 10);
 *   ll ans = linear_rec(s, berlekamp_massey(s), k);
 * Notes: Any sequence $s_k = u^T A^k v$ for a fixed $m times m$
 *   matrix has order at most $m$ (Cayley-Hamilton). The order of a
 *   sum is at most the sum of the orders, of a termwise product at
 *   most their product, and a prefix sum adds one; a polynomial of
 *   degree $d$ in $k$, or $k^d r^k$, has order $d + 1$. Cues: "a
 *   fixed transition repeated $k <= 10^18$ times", a DP over a small
 *   automaton, "the answer modulo a prime for huge $n$" where a
 *   matrix is awkward to write. Bound the hidden state size $m$,
 *   compute $2m + 10$ terms, run `berlekamp_massey` and check that
 *   $L$ is well below half the terms. Needs a field: no max-plus DP;
 *   for a composite modulus see Tweak.
 *   `mul` works modulo $x P(x)$, where
 *   $P(x) = x^L - sum_j c_j x^(L - j)$, so that $x$ itself needs no
 *   reduction even when $L = 1$.
 * Source: kactl BerlekampMassey and LinearRecurrence (CC0).
 * Status: stress-tested against Gaussian elimination and matrix power
 */
#pragma once
#include "../numtheory/modpow.hpp"

typedef vector<ll> vl;
vl berlekamp_massey(const vl& s) {
  int n = sz(s), L = 0, m = 0;
  vl C(n + 1), B{1}, T;
  C[0] = 1;
  ll b = 1;
  rep(i, 0, n) {
    ++m;
    ll d = s[i] % MOD;
    rep(j, 1, L + 1) d = (d + C[j] * s[i - j]) % MOD;
    if (!d) continue;
    T = vl(C.begin(), C.begin() + L + 1);
    ll coef = d * modinv(b) % MOD;
    rep(j, 0, sz(B)) C[j + m] = (C[j + m] - coef * B[j]) % MOD;
    if (2 * L > i) continue;
    L = i + 1 - L, B = T, b = d, m = 0;
  }
  C.resize(L + 1), C.erase(C.begin());
  for (ll& x : C) x = (MOD - x) % MOD;
  return C;
}
ll linear_rec(const vl& s, const vl& c, ll k) {
  int n = sz(c);
  auto mul = [&](const vl& a, const vl& b) { // a b mod x P(x)
    vl r(sz(a) + sz(b) - 1);
    rep(i, 0, sz(a))
      rep(j, 0, sz(b))
        r[i + j] = (r[i + j] + a[i] * b[j]) % MOD;
    for (int i = sz(r) - 1; i > n; i--)
      rep(j, 0, n)
        r[i - 1 - j] = (r[i - 1 - j] + r[i] * c[j]) % MOD;
    r.resize(min(sz(r), n + 1));
    return r;
  };
  vl r{1}; // ends as x^(k + 1) mod x P(x) = x (x^k mod P(x))
  for (int t = __lg(++k); t >= 0; t--) {
    r = mul(r, r);
    if (k >> t & 1) r = mul(r, {0, 1});
  }
  ll res = 0;
  rep(i, 1, sz(r)) res = (res + r[i] * s[i - 1]) % MOD;
  return (res + MOD) % MOD;
}
