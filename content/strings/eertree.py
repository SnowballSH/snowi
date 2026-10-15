"""
Title: Eertree
Description: Palindromic tree, built letter by letter. Its nodes are
  the distinct non-empty palindromic substrings of `s` plus two
  roots: node 0 of length $-1$ and node 1, the empty string. `nxt[26
  * v + c]` is the palindrome `c + v + c` (from node 0: just `c`), 0
  when absent, and `link[v]` the longest proper palindromic suffix
  of `v`. `suf[i]` is the node of the longest palindromic suffix of
  the prefix ending at `s[i]`. After `count()`, `cnt[v]` is the
  number of occurrences of `v` in `s`.
Precondition: Letters `'a'` to `'z'`.
Time: $O(n)$ amortized; at most $n + 2$ nodes. PyPy: about 0.3 s for
  $n = 2 dot 10^5$.
Usage:
  et = Eertree(s)  # or et.add(c), c in [0, 26)
  distinct = len(et.len) - 2  # palindromic substrings
  longest = et.len[et.suf[i]]  # ending at s[i]
  et.count()  # then et.cnt[v]: occurrences of node v >= 2
Tweak: Other alphabets: change every `26` and the offset 97.
Example:
  >>> et = Eertree("abacaba")
  >>> m = len(et.len)
  >>> m - 2, [et.len[v] for v in et.suf]
  (7, [1, 1, 3, 1, 3, 5, 7])
  >>> et.count()  # (length, occurrences) of each palindrome
  >>> sorted((et.len[v], et.cnt[v]) for v in range(2, m))
  [(1, 1), (1, 2), (1, 4), (3, 1), (3, 2), (5, 1), (7, 1)]
  >>> dep = [0] * m  # number of palindromic suffixes of each node
  >>> for v in range(2, m):
  ...   dep[v] = dep[et.link[v]] + 1
  >>> sum(dep[v] for v in et.suf)  # palindromes with repeats
  12
  >>> max(et.cnt[v] * et.len[v] for v in range(2, m))
  7
Notes: A string of length `n` has at most `n` distinct non-empty
  palindromic substrings, since each letter adds at most one new
  one: its longest palindromic suffix. Nodes are created in order of
  first occurrence, and `link[v] < v`, so loops over `v` in
  increasing order are top-down along links and decreasing order is
  bottom-up (as in `count()`). The lengths of the palindromic
  suffixes of any string form $O(log n)$ arithmetic progressions;
  with `diff[v] = len[v] - len[link[v]]` and the series link (the
  first node up the link path with a different `diff`) this gives
  $O(n log n)$ dynamic programming over partitions into palindromes
  (minimum number of pieces, number of ways). For "is `s[l, r)` a
  palindrome" queries use Manacher instead. Cues: "distinct
  palindromic substrings", "occurrences of each palindrome" (the
  number of occurrences times the length, maximized: Example),
  palindromes ending at each position, palindromic factorization.
Source: Rubinchik, Shur, "EERTREE" (2015).
Status: stress-tested against brute force
"""

from array import array


class Eertree:
  def __init__(self, s=""):
    self.len, self.link = [-1, 0], [0, 0]
    self.suf, self.cnt = [], []
    self.s = [-1]  # -1, then the letters
    self.nxt = array("i", bytes(4 * 26 * 2))
    for ch in s:
      self.add(ord(ch) - 97)

  def add(self, c):
    ln, link, nxt = self.len, self.link, self.nxt
    s, suf = self.s, self.suf
    s.append(c)
    i = len(s) - 1
    v = suf[-1] if suf else 1
    while s[i - ln[v] - 1] != c:
      v = link[v]
    w = nxt[26 * v + c]
    if not w:
      u = link[v]
      while s[i - ln[u] - 1] != c:
        u = link[u]
      link.append(nxt[26 * u + c] if ln[v] >= 0 else 1)
      w = len(ln)
      ln.append(ln[v] + 2)
      nxt.frombytes(bytes(4 * 26))
      nxt[26 * v + c] = w
    suf.append(w)

  def count(self):  # cnt[v]: occurrences of palindrome v
    link, cnt = self.link, [0] * len(self.len)
    for v in self.suf:
      cnt[v] += 1
    for v in range(len(cnt) - 1, 1, -1):
      cnt[link[v]] += cnt[v]
    self.cnt = cnt
