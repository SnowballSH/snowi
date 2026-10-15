"""
Title: Z-function
Description: `z[i]` is the length of the longest common prefix of
  `s` and its suffix `s[i, n)`. By convention `z[0] = n`.
Precondition: `s` is a `str`, `bytes` or list.
Time: $O(n)$. PyPy: about 0.1 s for $n = 10^6$.
Usage:
  z = z_function(s)
  # pat occurs at i in s iff zt[m + 1 + i] == m, m = len(pat):
  zt = z_function(pat + "#" + s)  # "#" in neither
Tweak: For lists use a separator value that occurs nowhere, like
  -1: `z_function(pat + [-1] + s)`. Many implementations use
  `z[0] = 0`; check before reusing a formula.
Example:
  >>> z_function("abaabaab")
  [8, 0, 1, 5, 0, 1, 2, 0]
  >>> s, pat = "abababa", "aba"
  >>> m, zt = len(pat), z_function(pat + "#" + s)
  >>> [i for i in range(len(s)) if zt[m + 1 + i] == m]
  [0, 2, 4]
  >>> s = "abaabaab"; n = len(s); z = z_function(s)
  >>> [p for p in range(1, n + 1) if p == n or z[p] == n - p]
  [3, 6, 8]
  >>> t = "abx"  # lcp(t, s[i, n)) for every i:
  >>> zt = z_function(t + "#" + s)
  >>> zt[len(t) + 1 :]
  [2, 0, 1, 2, 0, 1, 2, 0]
Notes: `z[i]` compares one suffix with the whole string; the
  prefix function `pi` compares prefixes with their borders. Both
  carry the same information: `p < n` is a period iff
  `z[p] = n - p`, iff `n - p` is on the border chain of `pi`;
  `p = n` is always a period and has no `z[p]`. To get `pi` from
  `z`: for each `i >= 1` set `pi[i + z[i] - 1] = max(.., z[i])`,
  then from right to left `pi[i] = max(pi[i], pi[i + 1] - 1)`.
  Reach for Z when you need the LCP of `s`, or of a pattern `t`,
  with every suffix at once: matching, the longest prefix of `t`
  starting at each position, comparing `s[i, n)` with `t` in
  $O(1)$ (LCP, then one character), and on reversed strings the
  longest common suffix. Distinct substrings in $O(n^2)$: after
  appending `c` to `s`, the new substrings are the prefixes of
  `u`, the reverse of `s + c`, longer than $max_(i >= 1) z_u [i]$.
Source: Gusfield, Algorithms on Strings, Trees and Sequences
  (1997); kactl Zfunc (CC0).
Status: stress-tested against brute force
"""


def z_function(s):
  n = len(s)
  z = [0] * n
  if n:
    z[0] = n
  l = r = 0
  for i in range(1, n):
    k = min(r - i, z[i - l]) if i < r else 0
    while i + k < n and s[k] == s[i + k]:
      k += 1
    z[i] = k
    if i + k > r:
      l, r = i, i + k
  return z
