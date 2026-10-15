"""
Title: Matrix power
Description: Matrices modulo `MOD` as lists of rows, with product
  and fast power. $(M^e)_(i j)$ counts the walks of exactly `e`
  steps from `i` to `j` when $M_(i j)$ counts the edges $i -> j$.
Precondition: Entries in $[0, "MOD")$, $e >= 0$. `mat_mul` takes
  any `n x m` and `m x p` matrices, `mat_pow` a square one.
Time: $O(n^3 log e)$; in PyPy $e < 2^60$ takes up to about 0.35 s
  for $n = 64$ and 1.2 s for $n = 100$ ($e = 10^18$, with fewer
  one bits: 0.25 s and 0.85 s).
Usage:
  m = [[0] * n for _ in range(n)]
  m[i][j] = 1
  walks = mat_pow(m, e)[s][t]
Tweak: (min, +) for the shortest walk with exactly `e` edges:
  start each `acc` row at `INF` instead of 0, replace `if x:` by
  `if x != INF:` (0 is a real weight here), update
  `acc[j] = min(acc[j], x + brow[j])`, and start `mat_pow` from 0
  on the diagonal and `INF` elsewhere. Entries $>= "INF"\/2$
  then mean no walk, as long as every walk weighs less than
  $"INF" \/ 4$ in absolute value. (max, +): `-INF` for `INF` and
  `max`; entries $<= -"INF"\/2$ mean no walk. Use a large int
  such as `1 << 62` for `INF`, not `float("inf")`.
Uses: numtheory/modpow
Example:
  >>> mat_pow([[1, 1], [1, 0]], 10)  # Fibonacci numbers
  [[89, 55], [55, 34]]
  >>> # walks of length 3 in the triangle 0-1-2
  >>> tri = [[0, 1, 1], [1, 0, 1], [1, 1, 0]]
  >>> mat_pow(tri, 3)[0]
  [2, 3, 3]
  >>> # f_n = f_(n-1) + 2 f_(n-3), f_0..2 = 1, 1, 1: the state
  >>> # (f_n, f_(n-1), f_(n-2)) is a column vector
  >>> m = [[1, 0, 2], [1, 0, 0], [0, 1, 0]]
  >>> [row[0] for row in mat_mul(mat_pow(m, 3), [[1], [1], [1]])]
  [7, 5, 3]
Notes: State machine DP repeated $10^18$ times: the state
  vector times the transition matrix. Walks of length at most
  `e` from `s` to `t`: add a state `z` with edges $t -> z$ and
  $z -> z$, then count walks $s -> z$ of length $e + 1$; for
  (min, +) put 0 on the diagonal instead. Sum
  $I + M + dots + M^(e-1)$: the top-right block of
  $mat(M, I; 0, I)^e$. Many queries with different `e`:
  precompute $M^(2^j)$ and multiply a vector by them in
  $O(n^2 log e)$ per query. Order-$n$ recurrences with
  $n >= 50$: `linear-rec` is $O(n^2 log e)$.
Source: Standard.
Status: stress-tested against repeated multiplication and walk
  counting
"""


def mat_mul(a, b):
  m, res, mod = len(b[0]) if b else 0, [], MOD
  for row in a:
    acc = [0] * m
    for x, brow in zip(row, b):
      if x:
        for j in range(m):
          acc[j] = (acc[j] + x * brow[j]) % mod
    res.append(acc)
  return res


def mat_pow(a, e):
  n = len(a)
  res = [[int(i == j) for j in range(n)] for i in range(n)]
  while e:
    if e & 1:
      res = mat_mul(res, a)
    a = mat_mul(a, a)
    e >>= 1
  return res
