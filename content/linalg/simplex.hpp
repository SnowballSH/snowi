/**
 * Title: Simplex (linear programming)
 * Description: Maximizes $c^T x$ subject to $A x <= b$, $x >= 0$.
 *   Returns the optimum and writes an optimal `x`, or returns `-inf`
 *   if no `x` is feasible and `inf` if the objective is unbounded.
 * Precondition: `A` is `m x n`, `b` has `m` entries (any signs),
 *   `c` has `n`. Data of moderate size: comparisons use `eps`.
 * Time: Exponential in the worst case, fast in practice: a random
 *   LP with 200 variables and 400 constraints takes about 0.15 s.
 * Usage:
 *   LPSolver lp(A, b, c);
 *   vector<double> x;
 *   double best = lp.solve(x); // -inf: infeasible, inf: unbounded
 * Tweak: `T = long double` on hard instances.
 * Example:
 *   // max 3x + 2y with x + y <= 4, x + 3y <= 7, x <= 3:
 *   LPSolver lp({{1, 1}, {1, 3}, {1, 0}}, {4, 7, 3}, {3, 2});
 *   vector<double> x;
 *   double best = lp.solve(x); // 11, x = {3, 1}
 *   // Dual prices: y_i = D[m][j] where N[j] = n + i, else 0.
 * Notes: Duality: $max {c^T x : A x <= b, x >= 0} =
 *   min {b^T y : A^T y >= c, y >= 0}$; the optimal dual `y` is in
 *   the objective row (see Example). The answer is real: an integer
 *   program needs more (branch and bound, or a known integral
 *   structure such as flows and bipartite matching, where the LP
 *   optimum is already integral). Steps use Dantzig's rule (most
 *   negative reduced cost) and switch to Bland's rule (smallest
 *   index) after a degenerate pivot, which cannot cycle. An
 *   infeasible start (some $b_i < 0$) runs a phase 1 with an
 *   artificial variable. Cues: "maximize a linear function under
 *   linear constraints" with few variables or constraints, LP
 *   duality to prove a bound or to turn min into max.
 * Source: kactl LPSolver (CC0), after Stanford ACM notebook.
 * Status: stress-tested against vertex enumeration and LP duality
 */
#pragma once

struct LPSolver {
  typedef double T;
  static constexpr T eps = 1e-9,
                     inf = numeric_limits<T>::infinity();
  int m, n;
  vi N, B; // nonbasic (columns) and basic (rows) variable ids
  vector<vector<T>> D;
  LPSolver(const vector<vector<T>>& A, const vector<T>& b,
           const vector<T>& c)
      : m(sz(b)), n(sz(c)), N(n + 1), B(m),
        D(m + 2, vector<T>(n + 2)) {
    rep(i, 0, m) {
      rep(j, 0, n) D[i][j] = A[i][j];
      B[i] = n + i, D[i][n] = -1, D[i][n + 1] = b[i];
    }
    rep(j, 0, n) N[j] = j, D[m][j] = -c[j];
    N[n] = -1, D[m + 1][n] = 1; // -1: the phase 1 variable
  }
  void pivot(int r, int s) {
    T inv = 1 / D[r][s];
    D[r][s] = 1;
    for (T& v : D[r]) v *= inv;
    rep(i, 0, m + 2)
      if (i != r && D[i][s] != 0) {
        T f = D[i][s];
        D[i][s] = 0;
        rep(j, 0, n + 2) D[i][j] -= f * D[r][j];
      }
    swap(B[r], N[s]);
  }
  bool simplex(int row) { // false if unbounded
    auto& z = D[row];
    for (bool bland = false;;) {
      int s = -1, r = -1;
      rep(j, 0, n + 1) {
        if (z[j] >= -eps || (row == m && N[j] == -1)) continue;
        if (s < 0 || (bland ? N[j] < N[s] : z[j] < z[s])) s = j;
      }
      if (s < 0) return true;
      T best = 0;
      rep(i, 0, m)
        if (D[i][s] > eps) {
          T q = D[i][n + 1] / D[i][s];
          if (r < 0 || pair(q, B[i]) < pair(best, B[r]))
            r = i, best = q;
        }
      if (r < 0) return false;
      bland = best <= eps; // degenerate step: Bland's rule next
      pivot(r, s);
    }
  }
  T solve(vector<T>& x) {
    int r = 0;
    rep(i, 1, m)
      if (D[i][n + 1] < D[r][n + 1]) r = i;
    if (D[r][n + 1] < -eps) {
      pivot(r, n);
      if (!simplex(m + 1) || D[m + 1][n + 1] < -eps)
        return -inf;
      rep(i, 0, m)
        if (B[i] == -1) { // still basic: swap it out
          int s = 0;
          rep(j, 1, n + 1)
            if (abs(D[i][j]) > abs(D[i][s])) s = j;
          pivot(i, s);
        }
    }
    bool ok = simplex(m);
    x.assign(n, 0);
    rep(i, 0, m)
      if (B[i] < n) x[B[i]] = D[i][n + 1];
    return ok ? D[m][n + 1] : inf;
  }
};
