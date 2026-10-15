"""
Title: Subset transforms
Description: For arrays indexed by bit masks:
  `bit_transform(a, "|")` replaces $a_S$ by the subset sum
  $sum_(T subset.eq S) a_T$ (zeta), `"&"` by the superset sum and
  `"^"` by the Walsh-Hadamard sum $sum_T (-1)^(|S inter T|) a_T$;
  `inv=True` undoes it (for `"|"`: Möbius, inclusion-exclusion).
  `bit_conv(a, b, op)` returns $c_S = sum a_i b_j$ over the pairs
  with `i op j == S`; `subset_conv` takes only disjoint pairs
  (`i | j == S`, `i & j == 0`).
Precondition: `a` and `b` have the same length $2^n$ (pad with
  zeros), entries in `[0, MOD)`; results are modulo `MOD`.
Time: $O(n 2^n)$ for the transforms and `bit_conv`, 0.1 to 0.2 s
  per transform at $n = 20$. $O(n^2 2^n)$ for `subset_conv`: 0.8 s
  at $n = 18$, but 3.5 s and 450 MB at $n = 20$ (C++: 1 s).
Usage:
  bit_transform(a, "|")  # a[S] = sum of a[T], T in S
  c = bit_conv(a, b, "^")  # c[i ^ j] += a[i] b[j]
  d = subset_conv(a, b)  # i | j = S with i & j = 0
Tweak: Exact integers without a modulus: replace the last three
  lines of `bit_transform` by
  `if inv and op == "^": a[:] = [x // n for x in a]` and drop
  `% MOD` in `bit_conv`. Max, min or OR over subsets instead of a
  sum: `a[j + s] = max(a[j + s], a[j])` in the `"|"` loop
  (`a[j] = max(a[j], a[j + s])` for `"&"`), delete the last three
  lines, no inverse. $k$-fold power: transform, raise each entry
  to the $k$-th power, invert.
Uses: numtheory/modpow
Example:
  >>> bit_conv([1, 2, 0, 3], [2, 0, 1, 1], "^")
  [5, 7, 3, 9]
  >>> bit_conv([1, 2, 0, 3], [2, 0, 1, 1], "|")
  [2, 4, 1, 17]
  >>> subset_conv([1, 2, 0, 3], [2, 0, 1, 1])  # disjoint only
  [2, 4, 1, 9]
  >>> cnt = [1, 1, 1, 1]  # pairs of 2-bit masks by their OR
  >>> bit_transform(cnt, "|")
  >>> cnt = [x * x for x in cnt]
  >>> bit_transform(cnt, "|", True); cnt
  [1, 3, 3, 9]
Notes: The loops leave values unreduced (at most $2^n$ times
  larger) and reduce once at the end. `subset_conv` keeps for each
  mask a polynomial in the rank $k$ (the size of the subset); a
  pair is disjoint exactly when the ranks add up to $|S|$. Ranks
  above $|S|$ are zero after the zeta transform, and only ranks
  from $|S|$ up are read back, so the loops skip the rest. Cues:
  "for every mask, sum over its submasks" (zeta); "exactly $S$"
  from "within $S$" (Möbius); pairs combined by OR, AND or XOR;
  "split the set into disjoint parts" or covers by non-overlapping
  pieces (subset convolution); the same random XOR move repeated
  $k$ times (power in the Walsh-Hadamard domain). OR is subsets,
  AND is supersets; check with $f = g = [S = "all"]$, whose AND
  convolution is $f$.
Source: kactl FastSubsetTransform (CC0); ranked zeta transform for
  subset convolution (folklore).
Status: stress-tested against $O(4^n)$ definitions
"""


def bit_transform(a, op, inv=False):
  n, s, sg = len(a), 1, -1 if inv else 1
  while s < n:
    for i in range(0, n, 2 * s):
      if op == "|":
        for j in range(i, i + s):
          a[j + s] += sg * a[j]
      elif op == "&":
        for j in range(i, i + s):
          a[j] += sg * a[j + s]
      else:
        for j in range(i, i + s):
          a[j], a[j + s] = a[j] + a[j + s], a[j] - a[j + s]
    s *= 2
  k = pow(n, MOD - 2, MOD) if inv and op == "^" else 1
  for i in range(n):
    a[i] = a[i] % MOD * k % MOD


def bit_conv(a, b, op):
  a, b = a[:], b[:]
  bit_transform(a, op)
  bit_transform(b, op)
  a = [x * y % MOD for x, y in zip(a, b)]
  bit_transform(a, op, True)
  return a


def subset_conv(a, b):
  n, W = len(a), len(a).bit_length()  # ranks 0..n-bits
  pc = [0] * n
  for S in range(1, n):
    pc[S] = pc[S >> 1] + (S & 1)
  A, B = [0] * (n * W), [0] * (n * W)  # A[S * W + k]: rank k
  for S in range(n):
    A[S * W + pc[S]], B[S * W + pc[S]] = a[S], b[S]
  s = 1
  while s < n:  # zeta per rank; ranks above |S| stay 0
    for i in range(0, n, 2 * s):
      for j in range(i, i + s):
        x, y = j * W, (j + s) * W
        for k in range(pc[j] + 1):
          A[y + k] += A[x + k]
          B[y + k] += B[x + k]
    s *= 2
  for i in range(n * W):
    A[i] %= MOD
    B[i] %= MOD
  for S in range(n):  # product of rank polynomials
    p, x = pc[S], S * W
    for k in range(min(2 * p, W - 1), p - 1, -1):
      t = 0
      for i in range(k - p, p + 1):
        t = (t + A[x + i] * B[x + k - i]) % MOD
      A[x + k] = t
  s = 1
  while s < n:  # Mobius on the ranks that are read: k >= |S|
    for i in range(0, n, 2 * s):
      for j in range(i, i + s):
        x, y = j * W, (j + s) * W
        for k in range(pc[j + s], W):
          A[y + k] -= A[x + k]
    s *= 2
  return [A[S * W + pc[S]] % MOD for S in range(n)]
