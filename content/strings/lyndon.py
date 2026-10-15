"""
Title: Lyndon factorization
Description: A Lyndon word is strictly smaller than all its
  proper suffixes, like `aab` or `abb`. `duval(s)` splits `s`
  into Lyndon words $w_1 >= w_2 >= dots >= w_k$ (the split is
  unique) and returns their start indices; the last factor is the
  smallest suffix of `s`. `min_rotation(s)` returns the smallest
  `i` for which the rotation `s[i, n) + s[0, i)` is
  lexicographically smallest.
Precondition: `s` is a `str`, `bytes` or list.
Time: $O(n)$. PyPy: about 0.1 s for $n = 10^6$.
Usage:
  st = duval(s)  # start of each factor
  suf = st[-1]  # smallest suffix: s[suf:], n >= 1
  st.append(len(s))  # now factor j is s[st[j]:st[j + 1]]
  r = min_rotation(s)
  t = s[r:] + s[:r]
Tweak: Largest rotation: the smallest one with the order flipped
  (`-x` for numbers, `chr(219 - ord(c))` for lowercase letters).
Example:
  >>> s = "abbabbaba"; st = duval(s)
  >>> st, [s[a:b] for a, b in zip(st, st[1:] + [len(s)])]
  ([0, 3, 6, 8], ['abb', 'abb', 'ab', 'a'])
  >>> min_rotation("baba"), min_rotation("cabcab")
  (1, 1)
  >>> len(duval("aabab")) == 1  # is it a Lyndon word?
  True
  >>> def canon(w):  # canonical form of a cyclic string
  ...   r = min_rotation(w); return w[r:] + w[:r]
  >>> canon("ballad"), canon("adball") == canon("ballad")
  ('adball', True)
Notes: Duval keeps `s[i, j)` as $w^m w'$ with $w$ Lyndon of
  length `j - k` and $w'$ a proper prefix of $w$: `s[j] > s[k]`
  merges everything into one Lyndon word, `s[j] == s[k]` extends
  the period, and `s[j] < s[k]` or the end emits the $m$ copies
  of $w$. Facts: a Lyndon word has no non-empty border, so it
  is primitive; if $u < v$ are Lyndon, $u v$ is Lyndon; `s`
  is Lyndon iff it is strictly smaller than all its other
  rotations. The minimal rotation is the canonical form of a
  cyclic string: two are equal up to rotation iff their minimal
  rotations are, and its start is unique iff `s` is primitive.
  Pitfall: the smallest suffix is not always the best start when
  more text follows: `a` is the smallest suffix of `aba`, yet
  `abac < ac`. Cues: "smallest rotation", "canonical form of a
  necklace", "lexicographically smallest result of a cut",
  "count substrings smaller than all their suffixes".
Source: Duval (1983); UCup Packet 19.
Status: stress-tested against brute force
"""


def duval(s):
  n, i, res = len(s), 0, []
  while i < n:
    j, k = i + 1, i
    while j < n and s[k] <= s[j]:
      k = i if s[k] < s[j] else k + 1
      j += 1
    while i <= k:
      res.append(i)
      i += j - k
  return res


def min_rotation(s):
  n, i, ans = len(s), 0, 0
  t = s + s
  while i < n:
    ans, j, k = i, i + 1, i
    while j < 2 * n and t[k] <= t[j]:
      k = i if t[k] < t[j] else k + 1
      j += 1
    while i <= k:
      i += j - k
  return ans
