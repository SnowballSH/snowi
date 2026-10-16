/**
 * Title: Fast Fourier transform
 * Description: `fft_conv(a, b)` multiplies polynomials with real
 *   coefficients in doubles. `fft(a)` replaces `a` in place by its
 *   values at the powers of $omega = e^((2 pi i) / n)$: `a[k]`
 *   becomes $sum_j a_j omega^(j k)$.
 * Precondition: `sz(a)` is a power of two. Results round exactly
 *   while $norm(a) norm(b) < 5 dot 10^14$ (Euclidean norms), which
 *   holds if $n dot max abs(a_i) dot max abs(b_j) < 5 dot 10^14$ for
 *   the longer length $n$, e.g. $2^19$ values up to $3 dot 10^4$.
 *   Small results do not help: with mixed signs the error still
 *   follows the inputs. Beyond, use `conv_mod`.
 * Time: $O(n log n)$. A $2^19 times 2^19$ product takes about 0.15 s.
 * Usage:
 *   vd c = fft_conv(a, b); // sz(a) + sz(b) - 1 entries
 *   ll x = llround(c[k]);
 * Tweak: Inverse transform: `fft(a)`, then
 *   `reverse(a.begin() + 1, a.end())` and divide by `n`.
 * Notes: `fft_conv` packs `a` and $2^e b$ into one complex vector
 *   and squares it, so a product costs two transforms. The packed
 *   error grows with $norm(a)^2 + norm(b)^2$, so $2^e$ (exact in
 *   floating point) evens out the norms; without it weights times a
 *   0/1 mask round wrongly. Roots come from `polar` directly rather
 *   than repeated multiplication, which keeps them accurate. Cues:
 *   integer convolutions whose exact values exceed the NTT modulus
 *   (counts of pair sums, wildcard string matching via
 *   $sum (x - y)^2 x y$), real-valued signals.
 * Source: kactl FastFourierTransform (CC0).
 * Status: stress-tested against the schoolbook product
 */
#pragma once

typedef complex<double> cd;
typedef vector<double> vd;
void fft(vector<cd>& a) {
  int n = sz(a), L = __lg(n);
  static vector<cd> rt(2, 1);
  for (static int k = 2; k < n; k *= 2) {
    rt.resize(n);
    rep(i, k, 2 * k) rt[i] = polar(1.0, M_PI * (i - k) / k);
  }
  vi rev(n);
  rep(i, 0, n) rev[i] = (rev[i / 2] | (i & 1) << L) / 2;
  rep(i, 0, n)
    if (i < rev[i]) swap(a[i], a[rev[i]]);
  for (int k = 1; k < n; k *= 2)
    for (int i = 0; i < n; i += 2 * k)
      rep(j, 0, k) {
        cd z = rt[j + k] * a[i + j + k];
        a[i + j + k] = a[i + j] - z;
        a[i + j] += z;
      }
}
vd fft_conv(const vd& a, const vd& b) {
  if (a.empty() || b.empty()) return {};
  vd res(sz(a) + sz(b) - 1);
  double na = inner_product(all(a), begin(a), 0.0);
  double nb = inner_product(all(b), begin(b), 0.0);
  if (!na || !nb) return res;
  int n = bit_ceil(unsigned(sz(res)));
  int e = (ilogb(na) - ilogb(nb)) / 2; // |a| ~ |2^e b|
  vector<cd> in(n), out(n);
  copy(all(a), begin(in));
  rep(i, 0, sz(b)) in[i].imag(ldexp(b[i], e));
  fft(in);
  for (cd& x : in) x *= x;
  rep(i, 0, n) out[i] = in[-i & (n - 1)] - conj(in[i]);
  fft(out);
  rep(i, 0, sz(res)) res[i] = ldexp(imag(out[i]), -e) / (4 * n);
  return res;
}
