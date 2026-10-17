"""
Title: Floor sum
Description: `floor_sum(n, m, a, b)` is
  $sum_(i=0)^(n-1) floor((a i + b) / m)$, for any signs of `a` and
  `b`, by a Euclid-like reduction.
Precondition: $n >= 0$, $m >= 1$.
Time: $O(log m)$; in PyPy about 1.5 µs per call for 32-bit
  inputs.
Usage:
  s = floor_sum(n, m, a, b)
Tweak: Sum over $i in [l, r)$: `floor_sum(r - l, m, a, a * l + b)`.
Example:
  >>> floor_sum(4, 3, 2, 1)  # 0 + 1 + 1 + 2
  4
  >>> floor_sum(5, 4, -3, 2)  # 0 - 1 - 1 - 2 - 3
  -7
  >>> n, m, a, b, c = 10, 7, 3, 5, 2  # i with (a i + b) % m < c
  >>> floor_sum(n, m, a, b) - floor_sum(n, m, a, b - c)
  3
  >>> sum(1 for i in range(n) if (a * i + b) % m < c)
  3
Notes: Counts lattice points under the line $y = (a x + b) / m$.
  Typical cues: a sum of `floor` or of `mod` over an arithmetic
  progression ($sum (a i + b) mod m$ is the plain sum minus `m`
  times the floor sum), counting `i` with $(a i + b) mod m < c$
  (indicator $[x mod m < c] = floor(x / m) - floor((x - c) / m)$),
  or comparing residues of two progressions, as in "Planting Trees"
  (UCup S4-X2 A). Each step swaps the roles of `a` and `m` like
  Euclid.
Source: AtCoder Library floor_sum (CC0).
Status: stress-tested against direct summation
"""


def floor_sum(n, m, a, b):
  res = a // m * (n * (n - 1) // 2) + b // m * n
  a, b = a % m, b % m
  while a * n + b >= m:
    n, b = divmod(a * n + b, m)
    a, m = m, a
    res += n * (n - 1) // 2 * (a // m) + n * (b // m)
    a, b = a % m, b % m
  return res
