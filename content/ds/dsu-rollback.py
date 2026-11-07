"""
Title: Union-find with rollback
Description: Union by size without path compression, so every
  `join` can be undone. `time()` returns a checkpoint;
  `rollback(t)` undoes all joins made after checkpoint `t`.
Precondition: Elements are `0..n-1`. Rollbacks are nested: `t`
  comes from an earlier `time()` not yet rolled past.
Time: $O(log n)$ per `find`/`join`, $O(1)$ amortized per undone
  join.
Usage:
  d = RollbackDSU(n)
  t = d.time()
  d.join(a, b); d.find(a); d.size(a)
  d.rollback(t)  # a and b are separate again
Tweak: Components: `n - time()`. Other state (a parity bit,
  per-set sums) must be recorded in `join` and restored in
  `rollback`, in reverse order.
Example:
  >>> d = RollbackDSU(4)
  >>> d.join(0, 1); t = d.time(); d.join(2, 3); d.join(1, 3)
  True
  True
  True
  >>> d.size(0), t, d.time()
  (4, 1, 3)
  >>> d.rollback(t); d.find(2) == d.find(0), d.size(0), d.size(3)
  (False, 2, 1)
Notes: Use it whenever a search or divide and conquer must undo
  unions: offline dynamic connectivity (segment tree over time),
  Mo's algorithm with rollback, or trying an edge set and then
  restoring. No path compression: it would make undo expensive,
  and union by size alone keeps depth at most $log_2 n$.
Source: kactl RollbackUF (CC0).
Status: stress-tested against brute force
"""


class RollbackDSU:
  def __init__(self, n):
    self.parent = list(range(n))
    self.sizes = [1] * n
    self.history = []  # roots b joined below another root

  def find(self, x):
    p = self.parent
    while p[x] != x:
      x = p[x]
    return x

  def size(self, x):
    return self.sizes[self.find(x)]

  def time(self):
    return len(self.history)

  def rollback(self, t):
    p, sizes, history = self.parent, self.sizes, self.history
    while len(history) > t:
      b = history.pop()
      sizes[p[b]] -= sizes[b]
      p[b] = b

  def join(self, a, b):
    a, b = self.find(a), self.find(b)
    if a == b:
      return False
    sizes = self.sizes
    if sizes[a] < sizes[b]:
      a, b = b, a
    self.parent[b] = a
    sizes[a] += sizes[b]
    self.history.append(b)
    return True
