"""
Title: String hashing
Description: Polynomial hash of `s[l, r)` in $O(1)$:
  $sum_(i = l)^(r - 1) s_i B^(r - 1 - i) mod M$ with the prime
  $M = 2^61 - 1$ and a random base $B$ shared by all instances,
  so substrings of different strings compare too. Different
  strings of length $n$ collide with probability at most
  $n / M approx n dot 4.3 dot 10^(-19)$.
Precondition: `s` is a `str`, `bytes` or list of ints in $[0, M)$.
  Compare hashes of equal lengths only.
Time: $O(n)$ to build, $O(1)$ per `get`. PyPy, $n = 10^6$: build
  0.15 s, $10^6$ calls to `get` 0.3 s.
Usage:
  hs = StrHash(s)
  eq = hs.get(a, a + ln) == hs.get(b, b + ln)
Tweak: Hash of `x` followed by `y`, where `y` has length
  `k <= n`: `(hs.mul(x, hs.pw[k]) + y) % hs.M`. Keys of mixed
  lengths: `(hash, len)`.
Example:
  >>> s = "abracadabra"; hs = StrHash(s)
  >>> hs.get(0, 4) == hs.get(7, 11), hs.get(0, 4) == hs.get(1, 5)
  (True, False)
  >>> def lcp(a, b):  # of s[a:] and s[b:], by binary search
  ...   lo, hi = 0, len(s) - max(a, b)
  ...   while lo < hi:
  ...     mid = (lo + hi + 1) // 2
  ...     if hs.get(a, a + mid) == hs.get(b, b + mid):
  ...       lo = mid
  ...     else:
  ...       hi = mid - 1
  ...   return lo
  >>> lcp(0, 7), lcp(3, 5), lcp(1, 2)
  (4, 1, 0)
  >>> rev, n = StrHash(s[::-1]), len(s)
  >>> def is_pal(l, r):  # s[l, r)
  ...   return hs.get(l, r) == rev.get(n - r, n - l)
  >>> is_pal(3, 6), is_pal(0, 4)  # "aca", "abra"
  (True, False)
Notes: Why the bound holds: for two different strings the
  difference of the polynomials has at most $n - 1$ roots among
  the $M$ bases, so $q$ comparisons fail with probability at most
  $q n / M$. Measured alternatives: two 31-bit primes are as fast
  but only guarantee $n^2 / 2^62$ per comparison; `x * B % M` on
  big ints is 1.5 to 4 times slower. Never hash modulo $2^64$:
  a Thue-Morse string of length 2048 and its complement collide
  for every odd base. Fixed bases are hackable, and a single
  modulus near $10^9$ fails by the birthday bound once about
  $10^5$ hashes are compared all against all. Under point updates
  keep the hashes in a Fenwick or segment tree. LCP in
  $O(log n)$ (Example) sorts suffixes or rotations in
  $O(n log^2 n)$, a cheap substitute for a suffix array. 2D
  patterns: hash each row, then hash the column of row hashes
  with a second base. Cues: "count distinct substrings of length
  `k`" (put the hashes in a set), many "is `s[l, r)` equal to
  `t[x, y)`" queries, huge implicitly defined strings (a block
  repeated $q$ times hashes to a geometric sum).
Source: Karp, Rabin (1987); $2^61 - 1$ multiply-reduce as in
  common competitive programming libraries.
Status: stress-tested against exact strings and big-int
  evaluation
"""

import random


class StrHash:
  M = (1 << 61) - 1
  B = random.randrange(M)

  @staticmethod
  def mul(a, b):  # a * b % M, every value stays below 2^63
    M = StrHash.M
    ah, al = a >> 31, a & 0x7FFFFFFF
    bh, bl = b >> 31, b & 0x7FFFFFFF
    m = al * bh + ah * bl
    x = (ah * bh << 1) + (m >> 30) + ((m & 0x3FFFFFFF) << 31)
    x += al * bl
    x = (x >> 61) + (x & M)
    return x - M if x >= M else x

  def __init__(self, s):
    if isinstance(s, str):
      s = [ord(c) for c in s]
    n, M, B, mul = len(s), self.M, self.B, self.mul
    self.h = h = [0] * (n + 1)
    self.pw = pw = [1] * (n + 1)
    for i in range(n):
      x = mul(h[i], B) + s[i]
      h[i + 1] = x - M if x >= M else x
      pw[i + 1] = mul(pw[i], B)

  def get(self, l, r):  # s[l, r)
    x = self.h[r] - self.mul(self.h[l], self.pw[r - l])
    return x + self.M if x < 0 else x
