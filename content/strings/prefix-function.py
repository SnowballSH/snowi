"""
Title: Prefix function (KMP)
Description: `pi[i]` is the length of the longest proper border
  of `s[0, i]`.
  `kmp_match(s, pat)` lists the start of every occurrence of `pat`
  in `s`, overlapping ones included.
Precondition: `pat` is non-empty. `s` and `pat` have the same type
  (`str`, `bytes` or list): a `str` character never equals a byte.
Time: $O(n)$; `kmp_match` $O(|s| + |"pat"|)$. PyPy: about 0.1 s for
  $n = 10^6$.
Usage:
  pi = prefix_function(s)
  per = n - pi[-1]  # smallest period, n >= 1
  occ = kmp_match(s, pat)  # start positions
Example:
  >>> prefix_function("abaabaab")
  [0, 0, 1, 1, 2, 3, 4, 5]
  >>> kmp_match("abababa", "aba"), kmp_match(b"xyz", b"y")
  ([0, 2, 4], [1])
  >>> s = "abaabaab"; n = len(s); pi = prefix_function(s)
  >>> n - pi[-1]  # period 3 does not divide 8: s is primitive
  3
  >>> k, borders = pi[-1], []  # all borders: walk the chain
  >>> while k:
  ...   borders.append(k); k = pi[k - 1]
  >>> borders
  [5, 2]
  >>> cnt = [0] * (n + 1)  # occurrences of each prefix s[0, L)
  >>> for x in pi: cnt[x] += 1
  >>> for L in range(n, 0, -1): cnt[pi[L - 1]] += cnt[L]
  >>> [c + 1 for c in cnt[1:]]
  [5, 3, 2, 2, 2, 1, 1, 1]
  >>> A = 2  # prefix automaton over "ab": nxt[k * A + c]
  >>> nxt = [0] * ((n + 1) * A)
  >>> for k in range(n + 1):
  ...   for c in range(A):
  ...     if k < n and ord(s[k]) - 97 == c:
  ...       nxt[k * A + c] = k + 1
  ...     elif k:
  ...       nxt[k * A + c] = nxt[pi[k - 1] * A + c]
  >>> nxt[5 * A + 0], nxt[5 * A + 1], nxt[8 * A + 0]
  (6, 0, 6)
Notes: The borders of `s[0, L)` form the chain `L -> pi[L - 1] ->
  pi[pi[L - 1] - 1] -> ... -> 0`; the step is `k = pi[k - 1]`,
  not `pi[k]`. They are the ancestors of `L` in the failure tree
  whose nodes are `0..n` and where the parent of `L` is
  `pi[L - 1]`; parents are smaller, so tree DP runs as a plain
  loop over `L` (as in the counting example) and no recursion is
  needed. Occurrences of every prefix of `s` in another text `t`:
  run `kmp_match`'s loop, add one to `cnt[k]` after each character
  (before falling back when `k = n`), then push counts to parents
  as above, without the final `+1`. The prefix automaton turns
  "the text must (not) contain `s`" into a DP over the matched
  length: count strings of length `m` avoiding `s` in
  $O(m n sigma)$, or with a matrix power for huge `m`. Some
  statements count the whole string as a border and some do not;
  check the samples. Cues: "prefix equals suffix", "shortest
  string whose repetition gives `s`", "the same piece at both
  ends".
Source: Knuth, Morris, Pratt (1977); kactl KMP (CC0).
Status: stress-tested against brute force
"""


def prefix_function(s):
  pi, k = [0] * len(s), 0
  for i in range(1, len(s)):
    c = s[i]
    while k and s[k] != c:
      k = pi[k - 1]
    if s[k] == c:
      k += 1
    pi[i] = k
  return pi


def kmp_match(s, pat):
  pi, res, k, m = prefix_function(pat), [], 0, len(pat)
  for i in range(len(s)):
    c = s[i]
    while k and pat[k] != c:
      k = pi[k - 1]
    if pat[k] == c:
      k += 1
      if k == m:
        res.append(i - m + 1)
        k = pi[k - 1]
  return res
