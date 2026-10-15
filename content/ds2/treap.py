"""
Title: Implicit treap
Description: A sequence with insert and erase at any position,
  range reverse, range sum and access by index, by splitting and
  merging a treap keyed by position. `split(x, k)` cuts tree `x`
  into its first `k` elements and the rest; `merge(a, b)`
  concatenates. Node 0 is the empty tree; node `x` has children
  `ch[2x]` and `ch[2x + 1]`, other fields in parallel lists.
Precondition: 0-indexed positions, half-open ranges `[l, r)`,
  `0 <= i <= size()` for `insert`, `i < size()` otherwise.
Time: $O(log n)$ expected per operation. No recursion: `split`
  and `merge` walk down, then fix the path bottom-up. At
  $n = 10^5$: $10^5$ inserts 0.6 s, $10^5$ reversals 1.5 s,
  $10^5$ range sums 0.4 s (cache misses dominate).
Usage:
  tr = Treap()
  tr.insert(i, v); tr.erase(i); tr.get(i)
  tr.prefix(k)  # sum of the first k elements
  tr.reverse(l, r); s = tr.query(l, r)  # sum of [l, r)
  out = tr.to_list()  # in order
  a, b = tr.split(tr.root, k)  # move the first k
  tr.root = tr.merge(b, a)  # elements to the end
Tweak: Other lazy tags (range add, range assign) go next to
  `flip`: apply them to `val` and `sm` (add `x * cnt` to `sm`)
  when set on a node, pass them to the children in `_push`. An
  ordered treap (a set) splits by key instead: compare `val` with
  the key where `split` compares `k` with the left size.
Example:
  >>> tr = Treap()
  >>> for i, v in enumerate([1, 2, 3, 4, 5]): tr.insert(i, v)
  >>> tr.reverse(1, 4); tr.to_list()
  [1, 4, 3, 2, 5]
  >>> tr.erase(0); tr.insert(2, 10); tr.to_list()
  [4, 3, 10, 2, 5]
  >>> tr.query(1, 4), tr.prefix(2), tr.get(2), tr.size()
  (15, 7, 10, 5)
  >>> a, bc = tr.split(tr.root, 1)  # cut [1, 3), paste at the end
  >>> b, c = tr.split(bc, 2)
  >>> tr.root = tr.merge(tr.merge(a, c), b); tr.to_list()
  [4, 2, 5, 3, 10]
Notes: Think of it when a sequence is edited in the middle: block
  moves, reversals, "insert after the k-th", rope-like strings.
  Node fields hold aggregates of the subtree, so a lazy tag must
  already be applied to the node that carries it (here `flip`
  swaps children lazily, which does not change `sm`). Reads
  (`get`, `prefix`, `query`) only walk down, which is several
  times cheaper than split and merge. Erased nodes are not
  reused. For "k-th element" or "insert in sorted order" only,
  `ds2/ordered-set` is simpler and faster in PyPy.
Source: Standard implicit treap (cp-algorithms); kactl Treap.
Status: stress-tested against a list
"""

from random import getrandbits


class Treap:
  def __init__(self):
    self.ch = [0, 0]  # children of x: ch[2 * x], ch[2 * x + 1]
    self.cnt, self.pri, self.val, self.sm = [0], [0], [0], [0]
    self.flip = [0]
    self.root = 0

  def _make(self, v):
    self.ch += (0, 0)
    self.cnt.append(1)
    self.pri.append(getrandbits(30))
    self.val.append(v)
    self.sm.append(v)
    self.flip.append(0)
    return len(self.cnt) - 1

  def _pull(self, x):
    ch, cnt, sm = self.ch, self.cnt, self.sm
    l, r = ch[2 * x], ch[2 * x + 1]
    cnt[x] = cnt[l] + 1 + cnt[r]
    sm[x] = sm[l] + self.val[x] + sm[r]

  def _push(self, x):
    flip = self.flip
    if flip[x]:
      ch = self.ch
      l, r = ch[2 * x + 1], ch[2 * x]
      ch[2 * x], ch[2 * x + 1] = l, r
      flip[l] ^= 1
      flip[r] ^= 1
      flip[x] = 0

  def split(self, x, k):  # (first k, rest)
    ch, cnt = self.ch, self.cnt
    path = []  # +x: x goes left, -x: x goes right
    while x:
      self._push(x)
      c = cnt[ch[2 * x]]
      if c < k:
        k -= c + 1
        path.append(x)
        x = ch[2 * x + 1]
      else:
        path.append(-x)
        x = ch[2 * x]
    a = b = 0
    for x in reversed(path):
      if x > 0:
        ch[2 * x + 1] = a
        a = x
      else:
        x = -x
        ch[2 * x] = b
        b = x
      self._pull(x)
    return a, b

  def merge(self, a, b):
    ch, pri = self.ch, self.pri
    path = []  # +a: a is the root, -b: b is the root
    while a and b:
      if pri[a] > pri[b]:
        self._push(a)
        path.append(a)
        a = ch[2 * a + 1]
      else:
        self._push(b)
        path.append(-b)
        b = ch[2 * b]
    x = a or b
    for y in reversed(path):
      if y > 0:
        ch[2 * y + 1] = x
      else:
        y = -y
        ch[2 * y] = x
      self._pull(y)
      x = y
    return x

  def insert(self, i, v):
    a, b = self.split(self.root, i)
    self.root = self.merge(self.merge(a, self._make(v)), b)

  def erase(self, i):
    a, bc = self.split(self.root, i)
    _, c = self.split(bc, 1)
    self.root = self.merge(a, c)

  def reverse(self, l, r):
    a, bc = self.split(self.root, l)
    b, c = self.split(bc, r - l)
    self.flip[b] ^= 1
    self.root = self.merge(self.merge(a, b), c)

  def prefix(self, k):  # sum of the first k elements
    ch, cnt, sm, x, s = self.ch, self.cnt, self.sm, self.root, 0
    while x:
      self._push(x)
      l = ch[2 * x]
      if k <= cnt[l]:
        x = l
      else:
        k -= cnt[l] + 1
        s += sm[l] + self.val[x]
        x = ch[2 * x + 1]
    return s

  def query(self, l, r):
    return self.prefix(r) - self.prefix(l)

  def get(self, i):
    ch, cnt, x = self.ch, self.cnt, self.root
    while True:
      self._push(x)
      left = cnt[ch[2 * x]]
      if i == left:
        return self.val[x]
      if i < left:
        x = ch[2 * x]
      else:
        i -= left + 1
        x = ch[2 * x + 1]

  def size(self):
    return self.cnt[self.root]

  def to_list(self):
    ch, val = self.ch, self.val
    out, stack, x = [], [], self.root
    while stack or x:
      while x:
        self._push(x)
        stack.append(x)
        x = ch[2 * x]
      x = stack.pop()
      out.append(val[x])
      x = ch[2 * x + 1]
    return out
