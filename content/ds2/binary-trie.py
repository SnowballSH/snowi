"""
Title: Binary trie
Description: Multiset of integers in $[0, 2^B)$ stored as a trie of
  their bits, highest bit first, with a count per node.
  `kth(k, x)` is the `k`-th smallest (0-indexed) of the values
  $y xor x$ over stored `y`; so `kth(k)` is the `k`-th smallest
  value, and `kth(0, x)` and `kth(size() - 1, x)` are the minimum
  and maximum xor with `x`. `count_less(x, k)` counts stored `y`
  with $(x xor y) < k$.
Precondition: $0 <= x, y < 2^B$; `add(y, -1)` needs `y` stored;
  `kth` needs `0 <= k < size()`; `count_less` takes any
  $k >= 0$.
Time: $O(B)$ per operation, about 0.5 s for $2 dot 10^5$ inserts
  with $B = 30$. Memory: 3 list slots (24 bytes) per node, up to
  $n B$ nodes.
Usage:
  tr = BinaryTrie(30)  # B = 30 bits
  tr.add(y); tr.add(y, -1)  # insert, erase one copy
  lo, hi = tr.min_xor(x), tr.max_xor(x)
  v = tr.kth(k); c = tr.count_less(x, k)
  copies = tr.count_less(y, 1)  # count of y
Tweak: Values up to $2^60$: `BinaryTrie(60)`. For many tries (one
  per subtree or prefix), share the lists and keep one root per
  trie, copying the path on `add` for a persistent version.
Example:
  >>> tr = BinaryTrie(3)
  >>> for y in [1, 4, 5, 5]: tr.add(y)
  >>> tr.max_xor(2), tr.min_xor(4), tr.kth(2), tr.size()
  (7, 0, 5, 4)
  >>> tr.count_less(1, 2), tr.count_less(5, 1)  # count of 5
  (1, 2)
  >>> tr.add(5, -2); [tr.kth(k, 3) for k in range(tr.size())]
  [2, 7]
  >>> a, p, best = [8, 1, 2, 12], 0, 0  # max xor of a subarray
  >>> tr = BinaryTrie(4); tr.add(0)
  >>> for v in a:
  ...   p ^= v
  ...   best = max(best, tr.max_xor(p))
  ...   tr.add(p)
  >>> best
  15
Notes: Think of it for "maximum xor of a pair", "maximum xor of a
  subarray" (prefix xors, as in the Example), "k-th smallest
  pairwise xor" (binary search on the answer with `count_less`),
  and online "minimum xor with any stored value". A walk follows
  the bit of `x` to keep a 0 in the answer and crosses to the
  other child for a 1; `kth` crosses when the 0-side holds at most
  `k` values. Erasing only decrements counts, so empty nodes stay;
  a count of 0 means missing. When a query may only use the `y`
  from an index range `[l, r)`, use a persistent trie (one root
  per prefix, subtract the counts of two versions) or sort the
  queries offline. "Max xor of any subset" is a different problem:
  use a linear basis over $FF_2$.
Source: Folklore.
Status: stress-tested against brute force
"""


class BinaryTrie:
  def __init__(self, bits=30):
    self.B = bits
    self.ch = [0] * 4  # children of v: ch[2v], ch[2v + 1]
    self.cnt = [0, 0]  # node 0 is null, node 1 the root

  def add(self, y, c=1):  # c copies; c < 0 erases
    ch, cnt = self.ch, self.cnt
    v = 1
    cnt[v] += c
    for i in range(self.B - 1, -1, -1):
      j = 2 * v + (y >> i & 1)
      if not ch[j]:
        ch[j] = len(cnt)
        cnt.append(0)
        ch += (0, 0)
      v = ch[j]
      cnt[v] += c

  def size(self):
    return self.cnt[1]

  def kth(self, k, x=0):  # k-th smallest y ^ x
    ch, cnt = self.ch, self.cnt
    v, res = 1, 0
    for i in range(self.B - 1, -1, -1):
      j = 2 * v + (x >> i & 1)
      c = cnt[ch[j]]
      if c <= k:
        k -= c
        j ^= 1
        res |= 1 << i
      v = ch[j]
    return res

  def min_xor(self, x):
    return self.kth(0, x)

  def max_xor(self, x):
    return self.kth(self.size() - 1, x)

  def count_less(self, x, k):  # number of y with x ^ y < k
    if k >> self.B:
      return self.size()
    ch, cnt = self.ch, self.cnt
    v, res = 1, 0
    for i in range(self.B - 1, -1, -1):
      j = 2 * v + (x >> i & 1)
      if k >> i & 1:
        res += cnt[ch[j]]
        j ^= 1
      v = ch[j]
      if not v:
        break
    return res
