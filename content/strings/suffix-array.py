"""
Title: Suffix array
Description: `sa` lists the start positions of the suffixes of `s`
  in sorted order and `rank` is its inverse. `lcp[i]` is the length
  of the longest common prefix of the suffixes starting at `sa[i]`
  and `sa[i + 1]`, for `i < n - 1`, so the LCP of the suffixes
  ranked $i < j$ is `min(lcp[i:j])`.
Precondition: `s` is a `str`, `bytes` or list of hashable, mutually
  comparable values (ints, tuples, not lists); they are compressed
  first.
Time: $O(n log n)$. PyPy: about 0.3 s for $n = 2 dot 10^5$.
Usage:
  S = SuffixArray(s)
  # S.sa[r]: start of the suffix of rank r
  d = n * (n + 1) // 2 - sum(S.lcp)  # distinct substrings
Tweak: Sorted cyclic shifts: build on `s + s` and keep the starts `<
  n`.
Example:
  >>> S = SuffixArray("banana")
  >>> S.sa, S.rank, S.lcp
  ([5, 3, 1, 0, 4, 2], [3, 2, 5, 1, 4, 0], [1, 3, 0, 0, 2])
  >>> a, b = 1, 3  # LCP of two suffixes: min of lcp in between
  >>> i, j = sorted((S.rank[a], S.rank[b]))
  >>> min(S.lcp[i:j])  # O(1) with a sparse table over S.lcp
  3
  >>> def first_rank(pred):  # first rank r with not pred(sa[r])
  ...   lo, hi = 0, len(S.sa)
  ...   while lo < hi:
  ...     mid = (lo + hi) // 2
  ...     if pred(S.sa[mid]): lo = mid + 1
  ...     else: hi = mid
  ...   return lo
  >>> s, pat, m = "banana", "an", 2  # occurrences of pat:
  >>> lo = first_rank(lambda i: s[i : i + m] < pat)
  >>> hi = first_rank(lambda i: s[i : i + m] <= pat)
  >>> sorted(S.sa[lo:hi])
  [1, 3]
  >>> n = len(s)
  >>> n * (n + 1) // 2 - sum(S.lcp)  # distinct substrings
  15
Notes: The LCP of two suffixes is the minimum of `lcp` between their
  ranks (Example). The suffixes sharing a prefix of length `>= k`
  with a given one form a block of consecutive ranks: find it by
  binary search on a sparse table, or process all `k` at once with a
  monotonic stack over `lcp` (the LCP array walks the suffix tree).
  Distinct substrings: $n(n+1)/2 - sum_i "lcp"_i$; the suffix of
  rank `i` brings $n - "sa"_i - "lcp"_(i-1)$ new ones, which also
  finds the `k`-th smallest. The longest repeated substring has
  length $max_i "lcp"_i$. Several strings: join them with separators
  that are distinct and occur nowhere (lists of ints, letters
  shifted above them); the longest common substring of two strings
  is the maximum `lcp[i]` over adjacent ranks from different
  strings, and of `k` strings a sliding window over ranks that
  covers all `k`. Cues: "number of distinct substrings", "longest
  substring occurring twice", many LCP or comparison queries between
  substrings. For one pattern or one comparison, hashing or Z is
  shorter; for a string that grows online, use the suffix automaton.
Source: Prefix doubling with radix sort after kactl SuffixArray
  (CC0); LCP by Kasai et al. (2001).
Status: stress-tested against brute force
"""


class SuffixArray:
  def __init__(self, s):
    vals = sorted(set(s))
    code = {c: i + 1 for i, c in enumerate(vals)}
    t = [code[c] for c in s] + [0]  # sentinel 0 ends the string
    n, x = len(t), t[:]  # sorts the cyclic shifts of t
    sa, y = list(range(n)), [0] * n
    lim, j, p = len(vals) + 1, 0, 0
    while p < n:
      p = j
      for i in range(j):
        y[i] = n - j + i
      for i in sa:
        if i >= j:
          y[p] = i - j
          p += 1
      ws = [0] * max(n, lim)
      for v in x:
        ws[v] += 1
      for i in range(1, lim):
        ws[i] += ws[i - 1]
      for i in range(n - 1, -1, -1):
        v = y[i]
        ws[x[v]] -= 1
        sa[ws[x[v]]] = v
      x, y = y, x
      p, x[sa[0]] = 1, 0
      for i in range(1, n):
        a, b = sa[i - 1], sa[i]
        if y[a] != y[b] or y[a + j] != y[b + j]:
          p += 1
        x[b] = p - 1
      lim, j = p, max(1, 2 * j)
    n -= 1
    sa = sa[1:]
    rank, lcp, k = [0] * n, [0] * max(n - 1, 0), 0
    for i in range(n):
      rank[sa[i]] = i
    for i in range(n):
      r = rank[i]
      if r + 1 < n:
        j = sa[r + 1]
        while t[i + k] == t[j + k]:  # stops at the sentinel
          k += 1
        lcp[r] = k
      if k:
        k -= 1
    self.sa, self.rank, self.lcp = sa, rank, lcp
