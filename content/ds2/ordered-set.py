"""
Title: Ordered set
Description: Sorted multiset with indexing, for PyPy (which has no
  balanced tree): elements are kept in sorted buckets of about
  $sqrt(n)$ elements. `s[k]` is the `k`-th smallest (0-indexed,
  negative counts from the end); `bisect_left(x)` counts elements
  `< x`, `bisect_right(x)` those `<= x`.
Precondition: Elements are comparable with `<` (ints, tuples,
  strings). `s[k]` and `pop(k)` need $-n <= k < n$.
Time: $O(sqrt(n))$ per operation; at $n = 2 dot 10^5$, about
  $8 dot 10^5$ operations per second.
Usage:
  s = SortedList(a)  # or SortedList()
  s.add(x); s.discard(x)  # one copy; remove(x) raises if absent
  x in s; s.count(x); len(s); list(s)
  k = s.bisect_left(x)  # elements < x
  s[k]; s[-1]; s.pop(0)  # k-th smallest, max, pop the min
Tweak: Set instead of multiset: `if x not in s: s.add(x)`. Keys
  with a custom order: store tuples `(key, x)`. Speed: `RATIO`
  trades scanning buckets against moving elements in a bucket.
Example:
  >>> s = SortedList([5, 1, 3, 3])
  >>> s.add(4); s.discard(1), s.discard(7)
  (True, False)
  >>> list(s), len(s)
  ([3, 3, 4, 5], 4)
  >>> s.count(3), s.bisect_left(4), s.bisect_right(4), 2 in s
  (2, 2, 3, False)
  >>> s[0], s[-1], s.pop(1), list(s)
  (3, 5, 3, [3, 4, 5])
  >>> i = s.bisect_left(4)  # successor and predecessor of 4
  >>> s[i] if i < len(s) else None, s[i - 1] if i else None
  (4, 3)
Notes: Think of it for C++ `set`/`multiset` jobs in Python: the
  sliding window median, the k-th smallest of a changing set,
  "nearest stored value to `x`" (`bisect_left` then index),
  counting smaller elements to the left (inversions) when values
  are not small. If values are small ints, or all values are
  known in advance, a Fenwick tree over (compressed) values with
  `lower_bound` for the k-th is faster. A heap with lazy deletion
  (a counter of deleted values) suffices when only the min or max
  is needed. A bucket splits when it outgrows `RATIO` times the
  number of buckets; empty buckets are dropped. Hash maps: dicts
  and sets are fine at ICPC, where nobody crafts tests against
  your code. Python hashes an int to its value modulo $2^61 - 1$,
  so on hackable judges (Codeforces) keys can be chosen to
  collide: store `x ^ R` for one random `R` (e.g.
  `random.getrandbits(60)`) instead of `x`.
Source: Bucket idea from tatyam's SortedMultiset (Python).
Status: stress-tested against a sorted list
"""

from bisect import bisect_left, bisect_right, insort
from itertools import chain


class SortedList:
  RATIO = 32

  def __init__(self, values=()):
    a = sorted(values)
    n = self.size = len(a)
    k = int((n * self.RATIO / 2) ** 0.5) + 1  # bucket size
    self.a = [a[i : i + k] for i in range(0, n, k)]

  def _bucket(self, x):  # first bucket whose last element >= x
    a = self.a
    for i in range(len(a) - 1):
      if x <= a[i][-1]:
        return i
    return len(a) - 1

  def add(self, x):
    if not self.a:
      self.a.append([x])
    else:
      i = self._bucket(x)
      b = self.a[i]
      insort(b, x)
      if len(b) > self.RATIO * len(self.a):
        self.a[i : i + 1] = [b[: len(b) >> 1], b[len(b) >> 1 :]]
    self.size += 1

  def _pop(self, i, j):
    b = self.a[i]
    x = b.pop(j)
    self.size -= 1
    if not b:
      del self.a[i]
    return x

  def discard(self, x):  # removes one copy; False if absent
    if not self.a:
      return False
    i = self._bucket(x)
    b = self.a[i]
    j = bisect_left(b, x)
    if j == len(b) or b[j] != x:
      return False
    self._pop(i, j)
    return True

  def remove(self, x):
    if not self.discard(x):
      raise ValueError(x)

  def _kth(self, k):
    if k < 0:
      k += self.size
    for i, b in enumerate(self.a):
      if k < len(b):
        return i, k
      k -= len(b)
    raise IndexError(k)

  def __getitem__(self, k):
    i, j = self._kth(k)
    return self.a[i][j]

  def pop(self, k=-1):
    return self._pop(*self._kth(k))

  def bisect_left(self, x):
    c = 0
    for b in self.a:
      if x <= b[-1]:
        return c + bisect_left(b, x)
      c += len(b)
    return c

  def bisect_right(self, x):
    c = 0
    for b in self.a:
      if x < b[-1]:
        return c + bisect_right(b, x)
      c += len(b)
    return c

  def count(self, x):
    return self.bisect_right(x) - self.bisect_left(x)

  def __contains__(self, x):
    if not self.a:
      return False
    b = self.a[self._bucket(x)]
    j = bisect_left(b, x)
    return j < len(b) and b[j] == x

  def __len__(self):
    return self.size

  def __iter__(self):
    return chain.from_iterable(self.a)
