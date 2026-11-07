"""
Title: Manacher
Description: All maximal palindromes. A palindrome `s[l, r)` has
  center `c = l + r - 1` in $[0, 2n - 1)$: even `c` is the
  character `c/2`, odd `c` the gap after character `(c-1)/2`.
  `plen[c]` is the length of the longest palindrome with center
  `c`.
Time: $O(n)$. PyPy: about 0.2 s for $n = 10^6$.
Usage:
  plen = manacher(s)  # length 2n - 1
  pal = plen[l + r - 1] >= r - l  # s[l, r), l < r
  c = max(range(len(plen)), key=plen.__getitem__)  # n >= 1
  l = (c + 1 - plen[c]) // 2  # longest: s[l : l + plen[c]]
Example:
  >>> s = "abaaba"; plen = manacher(s)
  >>> plen
  [1, 0, 3, 0, 1, 6, 1, 0, 3, 0, 1]
  >>> [plen[l + r - 1] >= r - l for l, r in ((1, 4), (0, 3))]
  [False, True]
  >>> sum((x + 1) // 2 for x in plen)  # palindromes, with repeats
  11
  >>> max(r for r in range(1, len(s) + 1) if plen[r - 1] >= r)
  6
Notes: The palindromes with center `c` have lengths `p, p - 2,
  ...` down to 1 or 2, where `p = plen[c]`. They start at every
  position of $[(c + 1 - p)/2, floor(c/2)]$ and end at every
  position of $[ceil(c/2), (c + p - 1)/2]$, so the number of
  palindromes starting or ending at each position comes from a
  difference array (pairs of disjoint palindromes, cuts into
  palindromes). A string has at most `n` distinct palindromic
  substrings: each position ends at most one new one, the longest
  palindromic suffix so far; to count or list them use the
  eertree. Cues: many "is `s[l, r)` a palindrome" queries, "split
  into two palindromes" (a prefix and a suffix test), "longest
  palindromic substring".
Source: Manacher (1975); one-array form after the classic
  `#`-interleaved version.
Status: stress-tested against brute force
"""


def manacher(s):
  n = len(s)
  plen = [0] * max(2 * n - 1, 0)
  lo, hi = 0, -1  # rightmost palindrome found: s[lo..hi]
  for c in range(2 * n - 1):
    k = 1 - (c & 1)
    if c <= 2 * hi:  # start from the mirror image
      k = min(plen[2 * (lo + hi) - c], 2 * hi - c + 1)
    i, j = (c - k + 1) >> 1, (c + k - 1) >> 1
    while i > 0 and j + 1 < n and s[i - 1] == s[j + 1]:
      i -= 1
      j += 1
    plen[c] = j - i + 1
    if j > hi:
      lo, hi = i, j
  return plen
