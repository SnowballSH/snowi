"""
Title: Gaussian elimination mod p
Description: `gauss(a, m)`: reduced row echelon form in place,
  pivots from the first `m` columns (later columns ride along);
  returns the pivot columns and, if `a` is square and `m = n`, the
  determinant.
Precondition: `MOD` prime, entries in $[0, "MOD")$. `a` is
  `n x m` with any `n`, `m`.
Time: $O(n m min(n, m))$. In PyPy `solve_linear` takes about
  0.3 s for $n = m = 300$ and 1.4 s for $n = m = 500$; `mat_inv`
  works on $n times 2n$ and takes 0.9 s and 4 s.
Usage:
  r, x = solve_linear(a, b, m)  # r = -1, x = None: no solution
  d = det(a)
  inv = mat_inv(a)  # None if singular
  piv = gauss(a, m)[0]  # last: in place, a becomes its RREF
Tweak: Over $FF_2$ store each row as one int (bit `j` = column
  `j`), test `a[i] >> c & 1` and eliminate with `a[i] ^= a[r]`:
  $2000 times 2000$ takes about 1.2 s.
Uses: numtheory/modpow
Example:
  >>> a = [[1, 2, 3], [2, 4, 6]]  # rank 1
  >>> piv, d = gauss(a, 3)
  >>> piv, a[0]
  ([0], [1, 2, 3])
  >>> # kernel: one vector per free column f
  >>> ker = []
  >>> for f in range(3):
  ...   if f not in piv:
  ...     v = [0] * 3
  ...     v[f] = 1
  ...     for i, c in enumerate(piv):
  ...       v[c] = -a[i][f] % MOD
  ...     ker.append(v)
  >>> ker == [[MOD - 2, 1, 0], [MOD - 3, 0, 1]]
  True
  >>> # x + y = 3, x - y = 1
  >>> solve_linear([[1, 1], [1, MOD - 1]], [3, 1], 2)
  (2, [2, 1])
  >>> solve_linear([[1, 1], [2, 2]], [1, 3], 2)  # inconsistent
  (-1, None)
  >>> det([[2, 1], [1, 3]]), mat_inv([[1, 1], [1, 1]])
  (5, None)
Notes: When consistent, $A x = b$ has $"MOD"^(m - r)$ solutions:
  the one returned (free variables 0) plus any kernel vector. Rank
  and solvability depend on the field: a matrix can be singular
  mod $p$ but not over the rationals. Cues: "count assignments
  satisfying linear conditions mod p" ($p^(m - r)$), "determinant
  / number of spanning trees mod p" (matrix-tree), random-weight
  rank tests (Schwartz-Zippel). A composite modulus has no
  inverses: for a determinant, clear each column with Euclid-style
  row subtractions in $O(n^3 log "MOD")$.
Source: Standard Gauss-Jordan; kactl MatrixInverse-mod (CC0).
Status: stress-tested against permutation determinants and minor
  ranks
"""


def gauss(a, m):
  n, r, det, piv, mod = len(a), 0, 1, [], MOD
  for c in range(m):
    if r == n:
      break
    p = r
    while p < n and not a[p][c]:
      p += 1
    if p == n:
      continue
    if p != r:
      a[p], a[r] = a[r], a[p]
      det = -det
    row, w = a[r], len(a[r])
    det = det * row[c] % mod
    inv = pow(row[c], mod - 2, mod)
    for j in range(c, w):
      row[j] = row[j] * inv % mod
    for i in range(n):
      ai = a[i]
      f = ai[c]
      if f and i != r:
        for j in range(c, w):
          ai[j] = (ai[j] - f * row[j]) % mod
    piv.append(c)
    r += 1
  return piv, (det if r == n else 0)


def det(a):
  return gauss([row[:] for row in a], len(a))[1]


def solve_linear(a, b, m):
  rows = [row + [v] for row, v in zip(a, b)]
  piv = gauss(rows, m)[0]
  r = len(piv)
  if any(row[m] for row in rows[r:]):
    return -1, None
  x = [0] * m
  for i, c in enumerate(piv):
    x[c] = rows[i][m]
  return r, x


def mat_inv(a):
  n = len(a)
  rows = [row + [0] * n for row in a]
  for i in range(n):
    rows[i][n + i] = 1
  if len(gauss(rows, n)[0]) < n:
    return None
  return [row[n:] for row in rows]
