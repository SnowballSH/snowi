"""
Title: Linear recurrences
Description: `berlekamp_massey(s)` returns the shortest `c` with
  $s_i = sum_(j = 1)^L c_j s_(i - j)$ for every $i >= L$, as
  `c[j - 1]` $= c_j$. `linear_rec(s, c, k)` returns $s_k$ from the
  first $L$ terms, via $x^k$ modulo the characteristic polynomial.
Precondition: Arithmetic modulo the prime `MOD`. Entries of `s`
  and `c` are below `MOD` in absolute value, so a hand-written
  `c = [2, -1]` is fine. A recurrence of order $L$ is found from
  $2L$ terms; give a margin, as fewer terms give a shorter wrong
  answer. $k >= 0$.
Time: $O(n^2)$ for $n$ terms; $O(L^2 log k)$ for `linear_rec`.
  About 0.7 s for $L = 1000$ and $k = 10^18$, 0.4 s for $n = 10^4$.
Usage:
  c = berlekamp_massey([1, 1, 3, 5, 11, 21])  # [1, 2]
  x = linear_rec([1, 1], c, 10**18)  # s_k, 0-indexed
Tweak: Squarefree composite modulus: rebind `MOD`, run modulo
  each prime factor and combine with CRT. A modulus with a square
  factor ($2^32$, $10^9$) needs a matrix power instead. For $L$ in
  the thousands with huge $k$ replace `mul` by NTT multiplication
  and division (or use Bostan-Mori).
Uses: numtheory/modpow
Example:
  >>> berlekamp_massey([i * i % MOD for i in range(8)])
  [3, 998244350, 1]
  >>> linear_rec([0, 1], [1, 1], 40)  # Fibonacci
  102334155
  >>> s = [0, 1, 1, 2, 3, 5, 8, 13]
  >>> linear_rec(s, berlekamp_massey(s), 10**18)
  23849548
Notes: Any sequence $s_k = u^T A^k v$ for a fixed $m times m$
  matrix has order at most $m$ (Cayley-Hamilton). The order of a
  sum is at most the sum of the orders, of a termwise product at
  most their product, and a prefix sum adds one; a polynomial of
  degree $d$ in $k$, or $k^d r^k$, has order $d + 1$. Cues: "a
  fixed transition repeated $k <= 10^18$ times", a DP over a small
  automaton, "the answer modulo a prime for huge $n$" where a
  matrix is awkward to write. Bound the hidden state size $m$,
  compute $2m + 10$ terms, run `berlekamp_massey` and check that
  $L$ is well below half the terms. Needs a field: no max-plus DP;
  for a composite modulus see Tweak.
  `mul` works modulo $x P(x)$, where
  $P(x) = x^L - sum_j c_j x^(L - j)$, so that $x$ itself needs no
  reduction even when $L = 1$.
Source: kactl BerlekampMassey and LinearRecurrence (CC0).
Status: stress-tested against Gaussian elimination and matrix power
"""


def berlekamp_massey(s):
  n, L, m, b = len(s), 0, 0, 1
  C, B = [1] + [0] * n, [1]
  for i in range(n):
    m += 1
    d = s[i] % MOD
    for j in range(1, L + 1):
      d = (d + C[j] * s[i - j]) % MOD
    if d == 0:
      continue
    T, coef = C[: L + 1], d * pow(b, MOD - 2, MOD) % MOD
    for j, x in enumerate(B):
      C[j + m] = (C[j + m] - coef * x) % MOD
    if 2 * L <= i:
      L, B, b, m = i + 1 - L, T, d, 0
  return [-x % MOD for x in C[1 : L + 1]]


def linear_rec(s, c, k):
  n = len(c)

  def mul(a, b):  # a b modulo x P(x)
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
      for j, y in enumerate(b):
        r[i + j] = (r[i + j] + x * y) % MOD
    for i in range(len(r) - 1, n, -1):
      x = r[i]
      for j, y in enumerate(c):
        r[i - 1 - j] = (r[i - 1 - j] + x * y) % MOD
    return r[: n + 1]

  r = [1]
  for bit in bin(k + 1)[2:]:  # r ends as x (x^k mod P)
    r = mul(r, r)
    if bit == "1":
      r = mul(r, [0, 1])
  return sum(x * y for x, y in zip(r[1:], s)) % MOD
