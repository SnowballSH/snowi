"""
Title: Suffix automaton
Description: Online automaton of the substrings of `s`: the paths
  from state 0 spell exactly the substrings. Each state `v` is a
  class of substrings with the same end positions: the suffixes of
  its longest member with lengths `len[link[v]] + 1` to `len[v]`.
  After `count()`, `cnt[v]` is how often each string of `v` occurs
  in `s`. Each `extend` adds `len[last] - len[link[last]]` new
  distinct substrings. Edges are flat: `nxt[26 * v + c]`, 0 when
  absent.
Precondition: Letters `'a'` to `'z'`.
Time: $O(n)$ amortized, at most $2n$ states and $3n$ edges. PyPy:
  about 0.5 s and 100 MB for $n = 2 dot 10^5$.
Usage:
  sam = SuffixAutomaton(s)  # or sam.extend(c), c in [0, 26)
  v = sam.find(t)  # -1: t is not a substring of s
  sam.count()  # once, after the last extend
  occ = sam.cnt[v] if v >= 0 else 0  # t non-empty
  d = sam.distinct()  # distinct non-empty, O(n)
Tweak: Other alphabets: change every `26` and the offset 97. `nxt`
  is an `array` of 4-byte ints: a list takes twice the memory and is
  slower.
Example:
  >>> sam = SuffixAutomaton("abcbc")
  >>> sam.distinct(), sam.find("cb") > 0, sam.find("ca")
  (12, True, -1)
  >>> sam.count()
  >>> [sam.cnt[sam.find(t)] for t in ("bc", "c", "abcbc")]
  [2, 2, 1]
  >>> on, d, ds = SuffixAutomaton(), 0, []
  >>> for ch in "abcbc":  # distinct substrings of each prefix
  ...   on.extend(ord(ch) - 97)
  ...   d += on.len[on.last] - on.len[on.link[on.last]]
  ...   ds.append(d)
  >>> ds
  [1, 3, 6, 9, 12]
  >>> t = "xbcbz"  # longest common substring with s = "abcbc"
  >>> v = l = best = 0
  >>> for ch in t:
  ...   c = ord(ch) - 97
  ...   while v and not sam.nxt[26 * v + c]:
  ...     v = sam.link[v]; l = sam.len[v]
  ...   if sam.nxt[26 * v + c]:
  ...     v = sam.nxt[26 * v + c]; l += 1
  ...   best = max(best, l)
  >>> best
  3
Notes: The links form a tree, the suffix tree of the reversed
  string: moving to `link[v]` drops letters from the front. Sorting
  the states by `len` orders both the links and the edges
  topologically, as in `count()`. Distinct substrings starting from
  each state: `ways[v] = 1 + sum ways[nxt[26 * v + c]]` in
  decreasing `len`; walk down with it for the `k`-th smallest
  substring. First occurrence: record `len[cur] - 1` when a state is
  created and copy it into clones; `t` then starts at `first[v] -
  |t| + 1`. The end positions of `v` are those of the non-clone
  states in its link subtree (small-to-large merging for position
  queries). Matching statistics (Example) take $O(|t|)$ for any text
  `t`. Several strings: join them with distinct separators (a larger
  alphabet). Cues: distinct substrings of a growing string (sum the
  increments of `extend`; `distinct()` after every letter is
  $O(n^2)$), occurrences of many query strings, longest common
  substring, `k`-th distinct substring. Prefer the suffix array when
  memory is tight or when LCPs of arbitrary suffixes are needed.
Source: Blumer et al. (1985); construction as in cp-algorithms.
Status: stress-tested against brute force
"""

from array import array


class SuffixAutomaton:
  def __init__(self, s=""):
    self.len, self.link, self.cnt = [0], [-1], [0]
    self.nxt, self.last = array("i", bytes(4 * 26)), 0
    for ch in s:
      self.extend(ord(ch) - 97)

  def extend(self, c):
    ln, link, cnt, nxt = self.len, self.link, self.cnt, self.nxt
    p = self.last
    cur = self.last = len(ln)
    ln.append(ln[p] + 1)
    link.append(0)
    cnt.append(1)
    nxt.frombytes(bytes(4 * 26))
    while p != -1 and not nxt[26 * p + c]:
      nxt[26 * p + c] = cur
      p = link[p]
    if p == -1:
      return
    q = nxt[26 * p + c]
    if ln[p] + 1 == ln[q]:
      link[cur] = q
      return
    cl = len(ln)
    ln.append(ln[p] + 1)
    link.append(link[q])
    cnt.append(0)
    nxt.extend(nxt[26 * q : 26 * q + 26])
    while p != -1 and nxt[26 * p + c] == q:
      nxt[26 * p + c] = cl
      p = link[p]
    link[q] = link[cur] = cl

  def count(self):  # cnt[v] = number of end positions of v
    ln, link, cnt = self.len, self.link, self.cnt
    m = len(ln)
    b, order = [0] * m, [0] * m  # order: states by len
    for l in ln:
      b[l] += 1
    for i in range(1, m):
      b[i] += b[i - 1]
    for v in range(m - 1, -1, -1):
      b[ln[v]] -= 1
      order[b[ln[v]]] = v
    for i in range(m - 1, 0, -1):
      v = order[i]
      cnt[link[v]] += cnt[v]

  def find(self, t):
    v, nxt = 0, self.nxt
    for ch in t:
      v = nxt[26 * v + ord(ch) - 97]
      if not v:
        return -1
    return v

  def distinct(self):
    ln, link = self.len, self.link
    return sum(ln[v] - ln[link[v]] for v in range(1, len(ln)))
