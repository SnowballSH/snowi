"""
Title: Slope trick
Description: A convex piecewise-linear function $f$ with
  integer slopes, stored as its minimum `min_f` and its
  breakpoints: `L` (max-heap) left of the minimum, where the
  slope drops by 1 per point going left, and `R` (min-heap)
  right of it. The minimum is reached on `[top_l(), top_r()]`
  (infinite means unbounded). Shifts are lazy.
Precondition: `slide` needs `a <= b`. Keep values below
  $2^63$ for speed.
Time: $O(log n)$ per hinge, $O(1)$ for the rest (clearing is
  amortised), $O(n)$ for `eval`.
Usage:
  st = SlopeTrick()  # f = 0
  st.add_abs(a)  # f += |x - a|; also add_x_minus_a,
  # add_a_minus_x (one-sided hinges), add_const
  st.prefix_min()  # f(x) = min of f(y) over y <= x
  st.slide(a, b)  # f(x) = min of f(y), x - b <= y <= x - a
  lo, hi, v = st.top_l(), st.top_r(), st.min_f
Tweak: Weight $c dot |x - a|$: call `c` times, or store
  `(point, count)` pairs and split a pair when only part of it
  moves. Start pinned at $x_0$ ($f = 0$ there, $infinity$
  elsewhere): push `K` copies of $x_0$ into both heaps, `K`
  above any slope added later. Domain $x >= 0$: add `K` times
  `add_a_minus_x(0)`, or show the optimum never leaves it.
Example:
  >>> st = SlopeTrick()  # fewest +-1 steps: non-decreasing
  >>> for x in [5, 2, 4, 1]:
  ...   st.prefix_min()
  ...   st.add_abs(x)
  >>> st.min_f
  6
  >>> st = SlopeTrick()  # sum of distances: a median
  >>> for x in [4, 1, 6]: st.add_abs(x)
  >>> st.min_f, st.top_l(), st.top_r(), st.eval(0)
  (5, 4, 4, 11)
  >>> st.slide(-2, 3)  # x may move 2 left or 3 right
  >>> st.top_l(), st.top_r(), st.eval(0)
  (2, 7, 7)
Notes: Think of it when a DP has one numeric state `x` (a
  value, position or time) and costs like $|x - a|$ or
  $max(0, x - a)$, with transitions "min over $y <= x$", "min
  over a window around $x$", "shift", or "add two functions".
  Signals: make an array monotone with fewest $plus.minus 1$
  steps (record `top_l()` per step and clamp backwards for an
  optimal array; strictly increasing: use `a[i] - i`;
  neighbours within `d`: `slide(-d, d)` in place of
  `prefix_min()`, rising by 0 to `d`: `slide(0, d)`), move at
  speed at most 1 and pay distances at given times
  (`slide(-t, t)`), equalise root-to-leaf path lengths on a tree
  (merge children's heaps small-to-large). Adding a right hinge
  pays `top_l - a` if `a` is left of the minimum, then moves the
  new largest left point to `R`; the left hinge mirrors it.
  `min` of two different convex functions or a concave term
  break convexity. Track `min_f`; `eval` is for the final answer
  only. A concave function of a count works the same way with
  roles swapped: store its marginal gains in a heap (buy low,
  sell high: for each price `p`, if the cheapest kept price
  `q < p`, gain `p - q`, pop `q`, push `p`; always push `p`).
Source: Packet 13 (UCup); "Slope trick" by maspy.
Status: stress-tested against brute force on a grid
"""

from heapq import heappush, heappushpop


class SlopeTrick:
  def __init__(self):
    self.min_f = self.add_l = self.add_r = 0
    self.L = []  # add_l - point: a max-heap of the points
    self.R = []  # point - add_r

  def top_l(self):
    return self.add_l - self.L[0] if self.L else -float("inf")

  def top_r(self):
    return self.R[0] + self.add_r if self.R else float("inf")

  def add_const(self, c):
    self.min_f += c

  def add_x_minus_a(self, a):  # f += max(0, x - a)
    L, add_l = self.L, self.add_l
    if L and add_l - L[0] > a:
      self.min_f += add_l - L[0] - a
    top = add_l - heappushpop(L, add_l - a)
    heappush(self.R, top - self.add_r)

  def add_a_minus_x(self, a):  # f += max(0, a - x)
    R, add_r = self.R, self.add_r
    if R and R[0] + add_r < a:
      self.min_f += a - R[0] - add_r
    top = heappushpop(R, a - add_r) + add_r
    heappush(self.L, self.add_l - top)

  def add_abs(self, a):
    self.add_x_minus_a(a)
    self.add_a_minus_x(a)

  def prefix_min(self):
    self.R = []

  def suffix_min(self):
    self.L = []

  def slide(self, a, b):
    self.add_l += a
    self.add_r += b

  def eval(self, x):
    res, add_l, add_r = self.min_f, self.add_l, self.add_r
    for v in self.L:
      res += max(0, add_l - v - x)
    for v in self.R:
      res += max(0, x - v - add_r)
    return res
