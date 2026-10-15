"""
Title: Gaussian elimination (real)
Description: Gauss-Jordan elimination on floats with partial
  pivoting. `gauss_real(a, m)` works like the mod p version: it
  reduces `a` using pivots in its first `m` columns and returns the
  pivot columns and, for square `a`, the determinant. `solve_real`
  classifies $A x = b$ as having no, one or infinitely many
  solutions and returns one of them. Only `gauss_real` changes its
  argument.
Precondition: `a` is `n x m` with any `n`, `m`, entries of moderate
  size: `eps` is an absolute threshold below which a pivot counts as
  zero.
Time: $O(n m min(n, m))$; in PyPy $n = m = 300$ takes about 0.25 s,
  $n = m = 500$ about 1 s.
Usage:
  k, x = solve_real(a, b, m)  # k: 0 none, 1 unique, 2 many
  d = det_real(a)
  piv = gauss_real(a, m)[0]  # last: in place, a becomes its RREF
Tweak: Inverse: append the identity as in the mod p version.
  Exact answers on small inputs: pass `Fraction` entries and
  `eps=0` to any of the three functions (free variables come back
  as `0.0`). Badly scaled input: divide each row by its largest
  entry first, or pass a larger `eps`.
Example:
  >>> solve_real([[2.0, 1.0], [1.0, 3.0]], [4.0, 7.0], 2)
  (1, [1.0, 2.0])
  >>> solve_real([[1.0, 1.0], [2.0, 2.0]], [1.0, 2.0], 2)
  (2, [1.0, 0.0])
  >>> solve_real([[1.0, 1.0], [2.0, 2.0]], [1.0, 3.0], 2)[0]
  0
  >>> det_real([[1.0, 2.0], [3.0, 4.0]])
  -2.0
  >>> # Random walk on 0..3, reflecting at 0: expected steps to
  >>> # reach 3. E_3 = 0, else E_v = 1 + average of E_w over the
  >>> # neighbours w. The system has cycles, so no plain DP.
  >>> n = 4
  >>> a = [[float(i == j) for j in range(n)] for i in range(n)]
  >>> b = [1.0] * (n - 1) + [0.0]
  >>> for v in range(n - 1):
  ...   nb = [1] if v == 0 else [v - 1, v + 1]
  ...   for w in nb:
  ...     a[v][w] -= 1 / len(nb)
  >>> [round(e, 9) for e in solve_real(a, b, n)[1]]
  [9.0, 8.0, 5.0, 0.0]
Notes: The answer type decides the method: a probability or
  expectation as a real number needs this, "mod p" answers need
  the mod p version (exact). `eps` near $10^(-9)$ suits entries
  around 1 and $n <= 500$; entries near $10^6$ need a larger one
  or rescaling. Partial pivoting keeps errors small for
  well-conditioned systems, but a nearly singular system can
  still flip between "unique" and "many": if the model says the
  solution is unique, trust that. Cues: Markov chains with cycles
  (expected time, hitting probabilities), electrical networks
  (potentials from Kirchhoff's laws), fitting a polynomial through
  points.
Source: Standard; kactl SolveLinear (CC0).
Status: stress-tested against exact rational elimination
"""


def gauss_real(a, m, eps=1e-9):
  n, r, det, piv = len(a), 0, 1, []
  for c in range(m):
    if r == n:
      break
    p = r
    for i in range(r + 1, n):
      if abs(a[i][c]) > abs(a[p][c]):
        p = i
    if abs(a[p][c]) <= eps:
      continue
    if p != r:
      a[p], a[r] = a[r], a[p]
      det = -det
    row, w = a[r], len(a[r])
    det *= row[c]
    inv = 1 / row[c]
    for j in range(c, w):
      row[j] *= inv
    for i in range(n):
      ai = a[i]
      f = ai[c]
      if f and i != r:
        for j in range(c, w):
          ai[j] -= f * row[j]
    piv.append(c)
    r += 1
  return piv, (det if r == n else 0)


def det_real(a, eps=1e-9):
  return gauss_real([row[:] for row in a], len(a), eps)[1]


def solve_real(a, b, m, eps=1e-9):
  rows = [row + [v] for row, v in zip(a, b)]
  piv = gauss_real(rows, m, eps)[0]
  r = len(piv)
  if any(abs(row[m]) > eps for row in rows[r:]):
    return 0, None
  x = [0.0] * m
  for i, c in enumerate(piv):
    x[c] = rows[i][m]
  return (2 if r < m else 1), x
