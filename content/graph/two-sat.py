"""
Title: 2-SAT
Description: Decides whether boolean variables `0..n-1` can satisfy
  all clauses $(a or b)$, and finds such values. A literal is `x`
  for "x is true" or `~x` for "x is false". `solve()` returns `None`
  if impossible, else a list of `n` bools.
Precondition: Literals use variables below `n` (or ones made by
  `add_var`).
Time: $O(n + m)$ for $m$ clauses; `at_most_one` of $k$ literals
  adds $k - 2$ variables and $3k - 5$ clauses.
Usage:
  ts = TwoSat(n)
  ts.add_clause(a, ~b)  # a or not b
  ts.implies(a, b); ts.add_clause(a, a)  # force a
  ts.at_most_one([a, ~b, c])
  val = ts.solve()  # None if unsatisfiable
Example:
  >>> ts = TwoSat(3)
  >>> ts.add_clause(0, 1); ts.add_clause(~0, 2)
  >>> ts.add_clause(~1, ~2); ts.add_clause(~2, ~0)
  >>> ts.solve()  # the only solution: 0 false, 1 true, 2 false
  [False, True, False]
  >>> ts.add_clause(0, ~1); ts.solve() is None
  True
  >>> ts = TwoSat(4); ts.at_most_one([0, 1, ~2, 3])
  >>> ts.add_clause(1, 1); ts.add_clause(2, 3)
  >>> ts.solve()[:4]
  [False, True, True, False]
Notes: Think 2-SAT when every object has exactly two options and
  every constraint involves at most two objects ("if A picks X then
  B must not pick Y"). "At least one of three" is not 2-SAT: look
  for a reason it is automatic or guess one element. A global
  constraint (one clique, one special vertex) often becomes pairwise
  after fixing one element: try each. Minimize a threshold: binary
  search plus 2-SAT. Equal/different constraints only: a parity DSU
  also counts solutions ($2^"components"$). Counting 2-SAT solutions
  in general is #P-hard. Forced literal `x`: `~x` reaches `x`. "At
  most one per group" with pairwise clauses is quadratic: use
  `at_most_one`.
Uses: graph/scc
Source: kactl 2-SAT (CC0); UCup S4-S5 Training Packet 18.
Status: stress-tested against exhaustive search
"""


class TwoSat:
  def __init__(self, n=0):
    self.n = n
    self.g = [[] for _ in range(2 * n)]

  def add_var(self):
    self.g += [], []
    self.n += 1
    return self.n - 1

  def add_clause(self, a, b):
    a = 2 * a if a >= 0 else 2 * ~a + 1
    b = 2 * b if b >= 0 else 2 * ~b + 1
    self.g[a ^ 1].append(b)
    self.g[b ^ 1].append(a)

  def implies(self, a, b):
    self.add_clause(~a, b)

  def at_most_one(self, lits):
    if len(lits) < 2:
      return
    pre = lits[0]  # "some of lits[0..i-1] is true"
    for i in range(1, len(lits)):
      self.add_clause(~pre, ~lits[i])
      if i + 1 < len(lits):
        nxt = self.add_var()
        self.implies(pre, nxt)
        self.implies(lits[i], nxt)
        pre = nxt

  def solve(self):
    comp, _ = scc(self.g)
    val = []
    for x in range(self.n):
      t, f = comp[2 * x], comp[2 * x + 1]
      if t == f:
        return None
      val.append(t < f)  # x comes later in topological order
    return val
