/**
 * Title: Gaussian elimination (real)
 * Description: Gauss-Jordan elimination on doubles with partial
 *   pivoting. `gauss_real(a, m)` works like the mod p version: it
 *   reduces `a` using pivots in its first `m` columns and returns the
 *   pivot columns and, for square `a`, the determinant. `solve_real`
 *   classifies $A x = b$ as having no, one or infinitely many
 *   solutions and returns one of them. Only `gauss_real` changes its
 *   argument.
 * Precondition: `a` is `n x m` with any `n`, `m`; size `x` to `m`
 *   before `solve_real`. Entries of moderate size: `eps` is an
 *   absolute threshold below which a pivot counts as zero.
 * Time: $O(n m min(n, m))$; $n = m = 500$ takes about 0.07 s.
 * Usage:
 *   vector<double> x(m);
 *   int k = solve_real(a, b, x); // 0 none, 1 unique, 2 many
 *   double d = det_real(a);
 *   // Last: gauss_real works in place, a becomes its RREF.
 *   vi piv = gauss_real(a, m).first; // rank sz(piv)
 * Tweak: Replace `double` by `long double` for more precision.
 *   Inverse: append the identity as in the mod p version.
 *   Ill-conditioned or badly scaled input: divide each row by its
 *   largest entry first, or pass a larger `eps` (all three
 *   functions take one).
 * Example:
 *   // Expected steps to absorption: E_v = 1 + sum_w P_vw E_w,
 *   // E_t = 0, i.e. (I - P) E = 1 on the transient states.
 *   vector<vector<double>> a(n, vector<double>(n));
 *   vector<double> b(n, 1), e(n);
 *   rep(v, 0, n) {
 *     a[v][v] = 1;
 *     if (v == t) { b[v] = 0; continue; }
 *     for (auto [w, pr] : out[v]) a[v][w] -= pr;
 *   }
 *   solve_real(a, b, e);
 * Notes: The answer type decides the method: a probability or
 *   expectation as a real number needs this, "mod p" answers need
 *   the mod p version (exact). `eps` near $10^(-9)$ suits entries
 *   around 1 and $n <= 500$; entries near $10^6$ need a larger one
 *   or rescaling. Partial pivoting keeps errors small for
 *   well-conditioned systems, but a nearly singular system can
 *   still flip between "unique" and "many": if the model says the
 *   solution is unique, trust that. Cues: Markov chains with cycles
 *   (expected time, hitting probabilities), electrical networks
 *   (potentials from Kirchhoff's laws), fitting a polynomial through
 *   points.
 * Source: Standard; kactl SolveLinear (CC0).
 * Status: stress-tested against exact rational elimination
 */
#pragma once

pair<vi, double> gauss_real(vector<vector<double>>& a, int m,
                            double eps = 1e-9) {
  int n = sz(a), r = 0;
  vi piv;
  double det = 1;
  for (int c = 0; c < m && r < n; c++) {
    int p = r;
    rep(i, r + 1, n)
      if (abs(a[i][c]) > abs(a[p][c])) p = i;
    if (abs(a[p][c]) <= eps) continue;
    if (p != r) swap(a[p], a[r]), det = -det;
    det *= a[r][c];
    double inv = 1 / a[r][c];
    rep(j, c, sz(a[r])) a[r][j] *= inv;
    rep(i, 0, n)
      if (i != r && a[i][c] != 0) {
        double f = a[i][c];
        rep(j, c, sz(a[i])) a[i][j] -= f * a[r][j];
      }
    piv.push_back(c), r++;
  }
  return {piv, r == n ? det : 0};
}

double det_real(vector<vector<double>> a, double eps = 1e-9) {
  return gauss_real(a, sz(a), eps).second;
}

int solve_real(vector<vector<double>> a,
               const vector<double>& b, vector<double>& x,
               double eps = 1e-9) {
  int n = sz(a), m = sz(x);
  rep(i, 0, n) a[i].push_back(b[i]);
  vi piv = gauss_real(a, m, eps).first;
  int r = sz(piv);
  rep(i, r, n)
    if (abs(a[i][m]) > eps) return 0;
  x.assign(m, 0);
  rep(i, 0, r) x[piv[i]] = a[i][m];
  return r < m ? 2 : 1;
}
