/**
 * Title: Subset transforms
 * Description: On arrays indexed by masks, `bit_transform` with
 *   `'|'` gives subset sums $sum_(T subset.eq S) a_T$ (zeta; `inv`:
 *   Möbius), `'&'` superset sums, `'^'` Walsh-Hadamard
 *   $sum_T (-1)^(|S inter T|) a_T$. `bit_conv`: $c_S = sum a_i b_j$
 *   over `i op j == S`; `subset_conv` only over disjoint `i`, `j`.
 * Precondition: `a` and `b` have the same length $2^n$ (pad with
 *   zeros), entries in `[0, MOD)`; results are modulo `MOD`.
 * Time: $O(n 2^n)$ for the transforms and `bit_conv`; $O(n^2 2^n)$
 *   for `subset_conv`, 1 s and 200 MB for $n = 20$.
 * Usage:
 *   bit_transform(a, '|', 0); // a[S] = sum of a[T], T in S
 *   vl c = bit_conv(a, b, '^'); // c[i ^ j] += a[i] b[j]
 *   vl d = subset_conv(a, b); // i | j = S with i & j = 0
 * Tweak: Exact integers without a modulus: replace the last two
 *   lines of `bit_transform` by
 *   `if (inv && op == '^') for (ll& x : a) x /= n;` and drop
 *   `% MOD` in `bit_conv`. Max, min or OR over subsets instead of a
 *   sum: `v = max(v, u)` in the `'|'` line (`u = max(u, v)` for
 *   `'&'`), delete the last two lines, no inverse.
 * Example:
 *   // Number of k-tuples of masks whose OR is exactly S:
 *   bit_transform(cnt, '|', 0);
 *   for (ll& x : cnt) x = modpow(x, k);
 *   bit_transform(cnt, '|', 1);
 * Notes: `bit_transform` leaves values unreduced (at most $2^n$
 *   times larger) and reduces once at the end. `subset_conv` keeps
 *   for each mask a polynomial in the rank $k$ (the size of the
 *   subset); a pair is disjoint exactly when the ranks add up to
 *   $|S|$. Ranks above $|S|$ are zero after the zeta transform, and
 *   only ranks from $|S|$ up are read back, so the loops skip the
 *   rest. Cues: "for every mask, sum over its submasks" (zeta);
 *   "exactly $S$" from "within $S$" (Möbius); pairs combined by OR,
 *   AND or XOR; "split the set into disjoint parts" or covers by
 *   non-overlapping pieces (subset convolution, $n <= 20$); the
 *   same random XOR move repeated $k$ times (power in the
 *   Walsh-Hadamard domain). OR is subsets, AND is supersets; check
 *   with $f = g = [S = "all"]$, whose AND convolution is $f$.
 * Source: kactl FastSubsetTransform (CC0); ranked zeta transform for
 *   subset convolution (folklore).
 * Status: stress-tested against $O(4^n)$ definitions
 */
#pragma once
#include "../numtheory/modpow.hpp"

typedef vector<ll> vl;
void bit_transform(vl& a, char op, bool inv) {
  int n = sz(a);
  for (int s = 1; s < n; s *= 2)
    for (int i = 0; i < n; i += 2 * s)
      rep(j, i, i + s) {
        ll &u = a[j], &v = a[j + s];
        if (op == '|') v += inv ? -u : u;
        else if (op == '&') u += inv ? -v : v;
        else tie(u, v) = pair(u + v, u - v);
      }
  ll k = inv && op == '^' ? modinv(n) : 1;
  for (ll& x : a) x = (x % MOD + MOD) * k % MOD;
}
vl bit_conv(vl a, vl b, char op) {
  bit_transform(a, op, 0), bit_transform(b, op, 0);
  rep(i, 0, sz(a)) a[i] = a[i] * b[i] % MOD;
  bit_transform(a, op, 1);
  return a;
}
vl subset_conv(const vl& a, const vl& b) {
  int n = sz(a), W = __lg(n) + 1;
  vector<unsigned> A(n * W), B(n * W); // [mask * W + rank]
  vi pc(n);
  rep(S, 1, n) pc[S] = pc[S / 2] + S % 2;
  rep(S, 0, n) A[S * W + pc[S]] = a[S], B[S * W + pc[S]] = b[S];
  auto add = [](unsigned& x, unsigned y) {
    if ((x += y) >= MOD) x -= MOD;
  };
  // Zeta per rank; ranks above |j| stay 0.
  for (int s = 1; s < n; s *= 2)
    for (int j = s; j < n; j = (j + 1) | s)
      rep(k, 0, pc[j]) {
        add(A[j * W + k], A[(j ^ s) * W + k]);
        add(B[j * W + k], B[(j ^ s) * W + k]);
      }
  rep(S, 0, n) { // rank polynomial product, ranks >= |S|
    int p = pc[S];
    unsigned *x = &A[S * W], *y = &B[S * W];
    for (int k = min(2 * p, W - 1); k >= p; k--) {
      ll t = 0;
      rep(i, k - p, p + 1) t = (t + (ll)x[i] * y[k - i]) % MOD;
      x[k] = unsigned(t);
    }
  }
  // Mobius on the ranks that are read: k >= |j|.
  for (int s = 1; s < n; s *= 2)
    for (int j = s; j < n; j = (j + 1) | s)
      rep(k, pc[j], W)
        add(A[j * W + k], MOD - A[(j ^ s) * W + k]);
  vl c(n);
  rep(S, 0, n) c[S] = A[S * W + pc[S]];
  return c;
}
