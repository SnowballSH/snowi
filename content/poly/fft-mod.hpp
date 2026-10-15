/**
 * Title: Convolution modulo any m
 * Description: `conv_mod(a, b, m)` multiplies polynomials modulo any
 *   `m`, prime or not: `c[k]` is $sum_(i + j = k) a_i b_j mod m$. It
 *   splits each value into two 15-bit halves and needs four FFTs.
 * Precondition: $2 <= m < 2^31$, entries in `[0, m)`. Exact for
 *   result lengths up to $2^20$ when $m < 2^30$ (e.g. $10^9 + 7$),
 *   up to $2^19$ when $m < 2^31$; random inputs allow more.
 * Time: $O(n log n)$, about twice `fft_conv`. A $2^19 times 2^19$
 *   product takes about 0.3 s.
 * Usage:
 *   vl c = conv_mod(a, b, 1'000'000'007);
 * Tweak: For $m = 998244353$ use the faster `conv`. Three NTT primes
 *   with CRT are an exact alternative while $n m^2 < 2^86$.
 * Notes: Values are first shifted into $(-m/2, m/2]$ and split into
 *   balanced digits in $[-2^14, 2^14)$. Their products then cancel in
 *   the sums, which makes the doubles hundreds of times more accurate
 *   on random inputs than unsigned halves; values near 0 or $m$
 *   become small. The worst case is all values near $m/2$. Up to 40
 *   entries it multiplies directly.
 *   Cues: an answer modulo $10^9 + 7$ that needs a convolution, or a
 *   small modulus like 1009 in a divide-and-conquer product.
 * Source: kactl FastFourierTransformMod (CC0); balanced digits are a
 *   standard precision trick.
 * Status: stress-tested against the schoolbook product
 */
#pragma once
#include "fft.hpp"

typedef vector<ll> vl;
vl conv_mod(const vl& a, const vl& b, ll m) {
  if (a.empty() || b.empty()) return {};
  vl res(sz(a) + sz(b) - 1);
  if (min(sz(a), sz(b)) <= 40) {
    rep(i, 0, sz(a))
      rep(j, 0, sz(b))
        res[i + j] = (res[i + j] + a[i] * b[j]) % m;
    return res;
  }
  int n = bit_ceil(unsigned(sz(res))), h = 1 << 14, c = 2 * h;
  vector<cd> L(n), R(n), outl(n), outs(n);
  auto split = [=](ll x) {
    x -= x > m / 2 ? m : 0;
    ll lo = ((x + h) & (c - 1)) - h;
    return cd(double((x - lo) / c), double(lo));
  };
  rep(i, 0, sz(a)) L[i] = split(a[i]);
  rep(i, 0, sz(b)) R[i] = split(b[i]);
  fft(L), fft(R);
  rep(i, 0, n) {
    int j = -i & (n - 1);
    outl[j] = (L[i] + conj(L[j])) * R[i] / (2.0 * n);
    outs[j] = (L[i] - conj(L[j])) * R[i] / (2.0 * n) / 1i;
  }
  fft(outl), fft(outs);
  rep(i, 0, sz(res)) {
    ll x = llround(real(outl[i])), z = llround(imag(outs[i]));
    ll y = llround(imag(outl[i])) + llround(real(outs[i]));
    res[i] = ((x % m * c + y) % m * c + z) % m;
    res[i] += res[i] < 0 ? m : 0;
  }
  return res;
}
