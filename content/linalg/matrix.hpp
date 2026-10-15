/**
 * Title: Matrix power
 * Description: Square matrices modulo `MOD` with product and fast
 *   power. $(M^e)_(i j)$ counts the walks of exactly `e` steps from
 *   `i` to `j` when $M_(i j)$ counts the edges $i -> j$.
 * Precondition: Entries in $[0, "MOD")$, $e >= 0$.
 * Time: $O(n^3 log e)$; $n = 100$ takes up to about 0.4 s for
 *   $e < 2^60$ (0.3 s for $e = 10^18$, which has fewer one bits).
 * Usage:
 *   Matrix m(n);
 *   m.a[i][j] = 1;
 *   ll walks = m.pow(e).a[s][t];
 * Tweak: (min, +) for the shortest walk with exactly `e` edges: set
 *   `ZERO = LLONG_MAX / 4` (no walk), `ONE = 0` and the update to
 *   `r.a[i][j] = min(r.a[i][j], x + o.a[k][j])`. An entry
 *   $>= "ZERO"\/2$ then means no walk, as long as every walk
 *   weighs less than $"ZERO" \/ 4$ in absolute value. (max, +):
 *   `ZERO = LLONG_MIN / 4` and `max`; an entry $<= "ZERO"\/2$
 *   then means no walk.
 * Example:
 *   // Linear recurrence f_n = f_(n-1) + 2 f_(n-3) as a state
 *   // vector (f_n, f_(n-1), f_(n-2)) times a companion matrix:
 *   Matrix m(3);
 *   m.a[0] = {1, 0, 2};
 *   m.a[1][0] = m.a[2][1] = 1;
 *   Matrix p = m.pow(n - 2); // from (f_2, f_1, f_0)
 *   ll fn = (p.a[0][0] * f2 + p.a[0][1] * f1 + p.a[0][2] * f0)
 *           % MOD;
 * Notes: State machine DP repeated $10^18$ times: the state
 *   vector times the transition matrix. Walks of length at most
 *   `e` from `s` to `t`: add a state `z` with edges $t -> z$ and
 *   $z -> z$, then count walks $s -> z$ of length $e + 1$; for
 *   (min, +) put 0 on the diagonal instead. Sum
 *   $I + M + dots + M^(e-1)$: the top-right block of
 *   $mat(M, I; 0, I)^e$. Many queries with different `e`:
 *   precompute $M^(2^j)$ and multiply a vector by them in
 *   $O(n^2 log e)$ per query. Order-$n$ recurrences with
 *   $n >= 100$: `linear-rec` is $O(n^2 log e)$.
 * Source: Standard.
 * Status: stress-tested against repeated multiplication and walk
 *   counting
 */
#pragma once
#include "../numtheory/modpow.hpp"

struct Matrix {
  static constexpr ll ZERO = 0, ONE = 1; // semiring 0 and 1
  vector<vector<ll>> a;
  Matrix(int n, bool id = false) : a(n, vector<ll>(n, ZERO)) {
    if (id)
      rep(i, 0, n) a[i][i] = ONE;
  }
  Matrix operator*(const Matrix& o) const {
    int n = sz(a);
    Matrix r(n);
    rep(i, 0, n)
      rep(k, 0, n)
        if (ll x = a[i][k]; x != ZERO)
          rep(j, 0, n)
            r.a[i][j] = (r.a[i][j] + x * o.a[k][j]) % MOD;
    return r;
  }
  Matrix pow(ll e) const {
    Matrix r(sz(a), true), b = *this;
    for (; e; e >>= 1, b = b * b)
      if (e & 1) r = r * b;
    return r;
  }
};
