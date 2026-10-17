"""
Title: Aho-Corasick
Description: Automaton of a set of patterns on their trie.
  `nxt[26 * v + c]` is the move on letter `c`, `link[v]` the longest
  proper suffix of `v` in the trie, `words[v]` the number of
  patterns ending at the current position in state `v`.
Precondition: Non-empty patterns of letters `'a'` to `'z'`;
  duplicates are fine.
Time: $O(26 L)$ to build, $L$ the total pattern length; $O(|t| + L)$
  for `count(t)`. PyPy, 26 letters: about 0.3 s and 25 MB for $L =
  |t| = 2 dot 10^5$, 2 s and 110 MB for $L = |t| = 10^6$.
Usage:
  ac = AhoCorasick(pats)  # list of str
  occ = ac.count(t)  # occ[i]: occurrences of pats[i]
  v = total = 0  # matches of all patterns in t:
  for ch in t:
    v = ac.nxt[26 * v + ord(ch) - 97]
    total += ac.words[v]
Example:
  >>> ac = AhoCorasick(["he", "she", "his", "hers", "he"])
  >>> ac.count("ushers")
  [1, 1, 0, 1, 1]
  >>> m = len(ac.link)  # every (pattern, end position) match
  >>> pid, out = [-1] * m, [-1] * m
  >>> for i, u in enumerate(ac.node): pid[u] = i
  >>> for u in ac.q:  # out[u]: nearest pattern end on link path
  ...   out[u] = u if pid[u] >= 0 else out[ac.link[u]]
  >>> v, res = 0, []
  >>> for i, ch in enumerate("ushers"):
  ...   v = ac.nxt[26 * v + ord(ch) - 97]
  ...   u = out[v]
  ...   while u >= 0:
  ...     res.append((pid[u], i))
  ...     u = out[ac.link[u]]
  >>> res  # duplicate patterns report one index
  [(1, 3), (4, 3), (3, 5)]
Notes: The links form the failure tree: `words[v]` sums the pattern
  ends on the path from `v` to the root (BFS order `q`), and `count`
  sums the hits over subtrees (reverse BFS order). Listing every
  match (Example) can output $Theta(|t| sqrt(L))$ pairs; count when
  you can. Dynamic programming over the automaton handles "strings
  of length `m` that contain no pattern" (states with `words[v] = 0`
  only; matrix power for huge `m`), "maximum total score of
  matches", and "shortest string containing all patterns" (BFS over
  state and a mask of patterns seen; with `words` replaced by masks
  OR-ed along links). Queries "how often does pattern `i` occur in
  text `j`" offline: mark the text states and sum over the subtree
  of the pattern node in the failure tree (Euler tour and a Fenwick
  tree). Cues: many patterns against one text, "forbidden
  substrings", counting strings that avoid or contain given words.
  One pattern: use the prefix function.
Source: Aho, Corasick (1975); goto automaton as in kactl AhoCorasick
  (CC0).
Status: stress-tested against brute force
"""

from array import array


class AhoCorasick:
  def __init__(self, pats):
    nxt = array("i", bytes(4 * 26))  # trie, then all moves
    node = []
    for p in pats:
      v = 0
      for ch in p:
        i = 26 * v + ord(ch) - 97
        if not nxt[i]:
          nxt[i] = len(nxt) // 26
          nxt.frombytes(bytes(4 * 26))
        v = nxt[i]
      node.append(v)
    m = len(nxt) // 26
    link, words, q = [0] * m, [0] * m, [0]
    for v in node:
      words[v] += 1
    for v in q:  # BFS: q grows while we iterate
      lv = link[v]
      if v:
        words[v] += words[lv]
      for c in range(26):
        u, w = nxt[26 * v + c], nxt[26 * lv + c] if v else 0
        if u:
          link[u] = w
          q.append(u)
        else:
          nxt[26 * v + c] = w
    self.nxt, self.link, self.words = nxt, link, words
    self.node, self.q = node, q

  def count(self, t):
    nxt, link, q = self.nxt, self.link, self.q
    hits, v = [0] * len(link), 0
    for ch in t:
      v = nxt[26 * v + ord(ch) - 97]
      hits[v] += 1
    for i in range(len(q) - 1, 0, -1):
      hits[link[q[i]]] += hits[q[i]]
    return [hits[u] for u in self.node]
