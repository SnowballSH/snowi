"""
Title: Simplex (linear programming)
Description: Maximizes $c^T x$ subject to $A x <= b$, $x >= 0$.
  `solve()` returns the optimum and an optimal `x`, or `(-inf, None)`
  if no `x` is feasible and `(inf, None)` if the objective is
  unbounded.
Precondition: `A` is `m x n`, `b` has `m` entries (any signs), `c`
  has `n`. Data of moderate size: comparisons use `eps`.
Time: Exponential in the worst case, fast in practice: a random
  LP with 100 variables and 200 constraints takes about 0.2 s in
  PyPy.
Usage:
  best, x = LPSolver(A, b, c).solve()  # best = inf or -inf: no x
Tweak: Exact answers on small inputs: pass `Fraction` data and
  `eps=0`. Minimize: negate `c` and the answer. Equality: two
  inequalities. Free variable: $x = x' - x''$ with both $>= 0$.
Example:
  >>> # max 3x + 2y with x + y <= 4, x + 3y <= 7, x <= 3
  >>> lp = LPSolver([[1, 1], [1, 3], [1, 0]], [4, 7, 3], [3, 2])
  >>> best, x = lp.solve()
  >>> round(best, 9), [round(v, 9) for v in x]
  (11.0, [3.0, 1.0])
  >>> # dual prices: y_i = D[m][j] where N[j] = n + i, else 0
  >>> m, n = 3, 2
  >>> y = [0.0] * m
  >>> for j, var in enumerate(lp.N):
  ...   if var >= n:
  ...     y[var - n] = lp.D[m][j]
  >>> [round(v, 9) for v in y]  # b . y = 11 too
  [2.0, 0.0, 1.0]
  >>> LPSolver([[1], [-1]], [1, -2], [1]).solve()  # 2 <= x <= 1
  (-inf, None)
  >>> LPSolver([[-1]], [0], [1]).solve()  # max x, x >= 0
  (inf, None)
Notes: Duality: $max {c^T x : A x <= b, x >= 0} =
  min {b^T y : A^T y >= c, y >= 0}$; the optimal dual `y` is in
  the objective row (see Example). The answer is real: an integer
  program needs more (branch and bound, or a known integral
  structure such as flows and bipartite matching, where the LP
  optimum is already integral). Steps use Dantzig's rule (most
  negative reduced cost) and switch to Bland's rule (smallest
  index) after a degenerate pivot, which cannot cycle. An
  infeasible start (some $b_i < 0$) runs a phase 1 with an
  artificial variable. Cues: "maximize a linear function under
  linear constraints" with few variables or constraints, LP
  duality to prove a bound or to turn min into max.
Source: kactl LPSolver (CC0), after Stanford ACM notebook.
Status: stress-tested against vertex enumeration and LP duality
"""


class LPSolver:
  def __init__(self, A, b, c, eps=1e-9):
    m, n = len(b), len(c)
    self.m, self.n, self.eps = m, n, eps
    self.N = list(range(n)) + [-1]  # -1: the phase 1 variable
    self.B = list(range(n, n + m))
    # 0 * b[i] - 1: the number type of b (Fraction stays exact)
    D = [list(A[i]) + [0 * b[i] - 1, b[i]] for i in range(m)]
    D.append([-v for v in c] + [0, 0])
    D.append([0] * n + [1, 0])
    self.D = D

  def pivot(self, r, s):
    D, w = self.D, self.n + 2
    row = D[r]
    inv = 1 / row[s]
    row[s] = 1
    for j in range(w):
      row[j] *= inv
    for i, di in enumerate(D):
      f = di[s]
      if f and i != r:
        di[s] = 0
        for j in range(w):
          di[j] -= f * row[j]
    self.B[r], self.N[s] = self.N[s], self.B[r]

  def simplex(self, row):  # False if unbounded
    m, n, eps = self.m, self.n, self.eps
    D, N, B = self.D, self.N, self.B
    z, bland = D[row], False
    while True:
      s = -1
      for j in range(n + 1):
        if z[j] >= -eps or (row == m and N[j] == -1):
          continue
        if s < 0 or (N[j] < N[s] if bland else z[j] < z[s]):
          s = j
      if s < 0:
        return True
      r, best = -1, 0
      for i in range(m):
        if D[i][s] > eps:
          q = D[i][n + 1] / D[i][s]
          if r < 0 or q < best or (q == best and B[i] < B[r]):
            r, best = i, q
      if r < 0:
        return False
      bland = best <= eps  # degenerate step: Bland's rule next
      self.pivot(r, s)

  def solve(self):
    m, n, eps, D, B = self.m, self.n, self.eps, self.D, self.B
    inf = float("inf")
    r = min(range(m), key=lambda i: D[i][n + 1], default=0)
    if D[r][n + 1] < -eps:
      self.pivot(r, n)
      if not self.simplex(m + 1) or D[m + 1][n + 1] < -eps:
        return -inf, None
      for i in range(m):
        if B[i] == -1:  # still basic: swap it out
          s = 0
          for j in range(1, n + 1):
            if abs(D[i][j]) > abs(D[i][s]):
              s = j
          self.pivot(i, s)
    if not self.simplex(m):
      return inf, None
    x = [0] * n
    for i in range(m):
      if B[i] < n:
        x[B[i]] = D[i][n + 1]
    return D[m][n + 1], x
