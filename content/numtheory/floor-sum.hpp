/**
 * Title: Floor sum
 * Description: `floor_sum(n, m, a, b)` is
 *   $sum_(i=0)^(n-1) floor((a i + b) / m)$, for any signs of `a` and
 *   `b`, by a Euclid-like reduction.
 * Precondition: $0 <= n < 2^32$, $1 <= m < 2^32$. The answer must
 *   fit in `ll`; intermediate values wrap modulo $2^64$ harmlessly.
 * Time: $O(log m)$, about 0.3 µs per call.
 * Usage:
 *   ll s = floor_sum(n, m, a, b);
 * Tweak: Sum over $i in [l, r)$: `floor_sum(r - l, m, a, a * l + b)`.
 *   For the answer modulo a prime, keep the same loop and reduce the
 *   products instead of relying on the wrap.
 * Example:
 *   ll F = floor_sum(n, m, a, b);
 *   // a, b >= 0: F = #points 0 <= x < n, 0 < y <= (a x + b) / m
 *   // sum of (a i + b) mod m over i in [0, n):
 *   ll rem = a * (n * (n - 1) / 2) + b * n - m * F;
 *   // #i in [0, n) with (a i + b) mod m < c, for 0 <= c <= m:
 *   ll cnt = F - floor_sum(n, m, a, b - c);
 * Notes: Counts lattice points under the line $y = (a x + b) / m$.
 *   Typical cues: a sum of `floor` or of `mod` over an arithmetic
 *   progression, counting `i` with $(a i + b) mod m < c$ (indicator
 *   $[x mod m < c] = floor(x / m) - floor((x - c) / m)$), or comparing
 *   residues of two progressions, as in "Planting Trees" (UCup
 *   S4-X2 A). Each step swaps the roles of `a` and `m` like Euclid.
 * Source: AtCoder Library floor_sum (CC0).
 * Status: stress-tested against direct summation
 */
#pragma once

ll floor_sum(ll n, ll m, ll a, ll b) {
  typedef unsigned long long ull;
  ull N = n, M = m, A = (a % m + m) % m, B = (b % m + m) % m;
  ull res = (a / m - (a % m < 0)) * (N * (N - 1) / 2) +
            (b / m - (b % m < 0)) * N;
  for (ull y; (y = A * N + B) >= M;) {
    N = y / M, B = y % M, swap(A, M);
    res += N * (N - 1) / 2 * (A / M) + N * (B / M);
    A %= M, B %= M;
  }
  return res;
}
