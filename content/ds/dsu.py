"""
Title: Union-find (DSU)
Description: Disjoint sets with union by size and path halving.
Precondition: Elements are `0..n-1`.
Time: $O(alpha(n))$ amortized per operation.
Usage:
  d = DSU(n)
  d.join(a, b)  # False if already in the same set
  d.same(a, b); d.find(a); d.size(a)
Example:
  >>> d = DSU(5)
  >>> d.join(0, 1), d.join(3, 4), d.join(1, 0), d.join(1, 4)
  (True, True, False, True)
  >>> d.same(0, 3), d.same(0, 2), d.size(4)
  (True, False, 4)
Notes: Think DSU when edges only get added and the questions are
  about connectivity, component sizes or counts: Kruskal, offline
  "delete edges" read backwards as additions, merging equal cells
  of a grid. For deletions in arbitrary order use the rollback
  version with offline divide and conquer over time. For parity or
  offset constraints ("x and y differ"), store the distance to the
  parent (weighted DSU), or use $2n$ nodes for "same/different".
Source: kactl UnionFind (CC0).
Status: stress-tested against brute force
"""


class DSU:
  def __init__(self, n):
    self.parent = list(range(n))
    self.sizes = [1] * n

  def find(self, x):
    p = self.parent
    while p[x] != x:
      p[x] = p[p[x]]
      x = p[x]
    return x

  def same(self, a, b):
    return self.find(a) == self.find(b)

  def size(self, x):
    return self.sizes[self.find(x)]

  def join(self, a, b):
    a, b = self.find(a), self.find(b)
    if a == b:
      return False
    sizes = self.sizes
    if sizes[a] < sizes[b]:
      a, b = b, a
    self.parent[b] = a
    sizes[a] += sizes[b]
    return True
