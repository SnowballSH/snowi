"""
Title: Aliens trick (WQS binary search)
Description: Minimum cost $F(k)$ using exactly `k` items when
  $F$ is convex in `k`, from a solver of the relaxed problem
  with no count limit where each item costs an extra `lam`.
  Finds the smallest `lam` whose solver count is at most `k`;
  then $F(k) = v - "lam" dot k$.
Precondition: `solve(lam)` returns `(v, c)`: the minimum `v` of
  cost + `lam` * count, and the smallest count `c` among
  solutions with value `v`. `k` is a feasible count, and every
  slope $F(j + 1) - F(j)$ lies in $[-"hi", -"lo"]$.
Time: $O(log("hi" - "lo"))$ calls to `solve`.
Usage:
  best = wqs(k, -C, C, solve)  # |slopes| <= C
Tweak: Ties toward more items: search the largest `lam` with
  count `>= k`. At most `k` items: `lo = 0`. Smallest count in the
  solver: `(v, c)` tuples; with states in lists, one int
  $v dot 2^20 + c$ ($|v| < 2^42$).
Example:
  >>> def solve(lam, a=(5, 1, 1, 5, 2)):  # max sum, non-adjacent
  ...   skip, take = (0, 0), (1 << 62, 0)  # (value, count)
  ...   for x in a:
  ...     t = (skip[0] + lam - x, skip[1] + 1)
  ...     skip, take = min(skip, take), t
  ...   return min(skip, take)
  >>> [-wqs(k, -14, 14, solve) for k in range(4)]
  [0, 5, 10, 8]
Notes: Why the tie rule: with ties toward fewer items, the
  reported count is the smallest optimal count, which never
  increases with `lam`, and it first drops to `k` or below at
  $"lam" = F(k) - F(k + 1)$, where `k` and `k + 1` are both
  optimal. So `k` is optimal there and $F(k) = v - "lam" k$ even
  if the solver reports fewer items; always use `k`, not `c`.
  Where $F$ is straight, many counts tie at one `lam`; a solver
  that breaks ties inconsistently (one comparison on values
  alone) makes the search overshoot to a `lam` where `k` is not
  optimal; comparing `(v, c)` (Tweak, Example) avoids it.
  Slopes can be as large as $F$ itself (in the Example,
  `a = 0, M, 0, ..., M, 0` with `k` near its maximum), so packed
  values may pass $2^63$: still correct, but big ints make
  `solve` about 6 times slower, while `(v, c)` tuples stay machine
  ints. Convexity sources: min-cost flow
  (`k` units), matroids (spanning tree with `k` white edges),
  partitions with a quadrangle-inequality cost (`dp/dc-opt`),
  exchange arguments.
  If unsure, print $F$ from a brute force and check
  $2 F(k) <= F(k - 1) + F(k + 1)$. Real-valued slopes: 60 rounds
  of bisection on a float `lam`. All $F(k)$ at once: query `lam`
  at chord slopes of known points; each call finds a new vertex
  or proves a segment straight.
Source: IOI 2016 "Aliens"; Packet 07 (UCup).
Status: stress-tested against exhaustive search
"""


def wqs(k, lo, hi, solve):
  while lo < hi:
    lam = (lo + hi) // 2
    if solve(lam)[1] <= k:
      hi = lam
    else:
      lo = lam + 1
  return solve(lo)[0] - lo * k
