#import "/typst/notes.typ": *

#full[Write the transition, then match its shape below. Check a
claimed property (monotone `opt`, convexity in $k$) by brute force
first.]

=== Which speed-up applies
#fact[Quadrangle inequality][$w(a, c) + w(b, d) <= w(a, d) + w(b, c)$
  for $a <= b <= c <= d$: crossing intervals cost at most nested
  ones. Holds for $f(x_r - y_l)$ with $f$ convex and $x$, $y$
  non-decreasing (so $f(S_r - S_l)$, or distances between two
  sorted point sets), and for sums over the pairs inside `[l, r)`
  of weights $>= 0$; it survives sums and adding $g(l) + h(r)$.
  Everywhere (any two rows and columns) it is a *Monge* matrix.]
#fact[Monotone `opt`][With the inequality, the smallest optimal
  $j$ in $min_j "old"[j] + w(j, i)$ never decreases in $i$: one
  layer in $O(n log n)$ (`dp/dc-opt`), or $O(n)$ by SMAWK. If
  `dp` depends on itself ($"dp"[i] = min_(j < i) "dp"[j] + w(j, i)$),
  use CDQ on $i$ or a deque of candidates with binary search.]
#fact[Knuth][$"dp"[l][r] = w(l, r) + min_k "dp"[l][k] + "dp"[k][r]$
  with the inequality and $w$ monotone under inclusion: $O(n^2)$
  (`dp/knuth`).]
#fact[Lines][$"dp"[i] = min_j ("dp"[j] + a_j b_i + c_j)$ plus terms
  in $i$ alone is a minimum of lines $a_j x + ("dp"[j] + c_j)$ at
  $x = b_i$; add line $j$ once $"dp"[j]$ is known. A deque if
  slopes and queries are sorted, else `ds/li-chao` or
  `ds/line-container`. No inequality needed.]
#fact[Convolution][$c_k = min_(i + j = k) (a_i + b_j)$. Both convex:
  $c_0 = a_0 + b_0$, then the differences of $a$ and of $b$ merged
  in sorted order, $O(n + m)$, and $c$ is convex. Only $b$ convex:
  the smallest optimal $i$ never decreases in $k$ (`dp/dc-opt`,
  SMAWK). Concave sequences under (max, +): the same.]
#fact[Exactly $k$][The best value is convex (minimum) or concave
  (maximum) in $k$: drop $k$, charge $lambda$ per item, binary
  search $lambda$ (`dp/wqs`). Convex: min-cost flow in the amount
  ($k$ disjoint subarrays, $k$ buy-sell trades), matroids,
  partitions with a quadrangle-inequality cost.]
#fact[Unit steps][$min sum_i f_i (b_i)$ over integers $b_i >= 0$
  with $sum_i b_i = m$ and every $f_i$ convex: take the $m$
  cheapest unit steps $f_i (x + 1) - f_i (x)$ over all $i$ (each
  $f_i$ yields them in increasing order). Huge $m$: binary search
  the price of the last step. Many $m$: the optimum is
  $max_lambda (lambda m - sum_i max_x (lambda x - f_i (x)))$, with
  $lambda$ over the step prices only.]
#fact[Ratios][The largest $sum a \/ sum b$ (all $b > 0$) over a
  family of sets is $>= lambda$ iff $max sum (a - lambda b) >= 0$:
  binary search $lambda$, or Dinkelbach: from any feasible ratio,
  set $lambda$ to the ratio of a set maximising
  $sum (a - lambda b)$ until it stops growing (a few rounds).]
#fact[Convex in a value][A DP over one numeric state $x$ that stays
  convex piecewise linear in $x$ (costs $|x - a|$, prefix minima,
  windows): `dp/slope-trick`.]
#fact[Block maximum][$"dp"[j] = min_(i < j) ("dp"[i] + max(a_i,
  dots, a_(j - 1)))$ in $O(n)$: a stack of groups of $i$ that share
  the current block maximum $M$, each with its least $"dp"[i]$ and
  the least $M + "dp"[i]$ over itself and the groups below it. For
  each $j$ push $i = j - 1$, merge into it the top groups whose
  maximum is $<= a_(j - 1)$; the top then holds $"dp"[j]$.]

=== Knapsack
#fact[Subset sums][Reachable sums up to $S$ in $O(n S / 64)$.#py[
  A big int is a fast bitset: about $10^9$ bit operations per
  0.3 s.]]
#lang-code(
  cpp: "bitset<S + 1> b; // b[s]: sum s is reachable
b[0] = 1;
for (int x : a) b |= b << x;",
  py: "b, mask = 1, (1 << (S + 1)) - 1  # bit s: sum s reachable
for x in a:
  b = (b | b << x) & mask",
)
- Sizes summing to $S$ take at most $sqrt(2 S)$ distinct values.
  If size $x$ occurs $c >= 3$ times, keep one copy if $c$ is odd,
  two if even, and add $floor((c - 1) / 2)$ copies of $2 x$; go
  from small $x$ up. $O(sqrt(S))$ items remain:
  $O(S sqrt(S) / 64)$.
- Bounded knapsack (count $c$, weight $w$, value $v$): split $c$
  into $1, 2, 4, dots$ while the next power fits in what is left,
  plus the remainder ($10 = 1 + 2 + 4 + 3$), as 0/1 items:
  $O(W sum log c_i)$. $O(n W)$: for each residue $r$ of $w$,
  $"new"[r + t w] = t v + max_(t - c <= s <= t) ("dp"[r + s w] - s v)$,
  a sliding window maximum (monotone deque).

=== Bitmasks and digits
- All submasks `s` of all masks `m`, $3^n$ in total: start at
  `s = m`, step `s = (s - 1) & m`, stop after `s = 0`.
- Sum over submasks in $O(n 2^n)$: for each bit `b`, for each
  mask `m` with bit `b`: `f[m] += f[m ^ 1 << b]`. Supersets,
  inverses, subset convolution: `poly/subset-transforms`.
- Broken profile (tilings of $n times m$, $m$ small): cells in
  row-major order, the profile holds the state of the last $m$
  cells (covered or not; base 3 when a cell needs three states):
  $O(n m s^m)$ for $s$ states per cell.
- Digit DP: count $x$ in $[0, N]$ whose digits drive an automaton
  (`nxt`, `Q` states) to `accept`; for $[L, R]$ take
  `count(R) - count(L - 1)`. Leading zeros are read as digits: add
  a "started" flag to the state when they matter.
#lang-code(
  cpp: "ll count(const string& N) { // x in [0, N]
  vector<ll> f(Q), g; // f[q]: prefixes already below N
  int t = 0;          // state of the prefix equal to N
  for (char ch : N) {
    g.assign(Q, 0);
    rep(q, 0, Q) rep(d, 0, 10) g[nxt(q, d)] += f[q];
    rep(d, 0, ch - '0') g[nxt(t, d)]++;
    f = g, t = nxt(t, ch - '0');
  }
  ll res = accept(t);
  rep(q, 0, Q) if (accept(q)) res += f[q];
  return res;
}",
  py: "def count(N):  # x in [0, N], N a decimal string
  f, t = [0] * Q, 0  # f: prefixes below N; t: the equal one
  for ch in map(int, N):
    g = [0] * Q
    for q in range(Q):
      for d in range(10):
        g[nxt(q, d)] += f[q]
    for d in range(ch):
      g[nxt(t, d)] += 1
    f, t = g, nxt(t, ch)
  return sum(f[q] for q in range(Q) if accept(q)) + accept(t)",
)

=== Orders and trees
- Exchange argument: if we choose the order (jobs, stacking),
  find when swapping neighbours `a, b` helps: `a` goes first iff
  `cost(a, b) <= cost(b, a)`, e.g. $t_a w_b < t_b w_a$ for
  weighted completion time. Sort by it, then run the DP or
  knapsack in that order. Cross-multiply instead of dividing and
  check that the order is transitive.
- Tree knapsack in $O(n^2)$, or $O(n K)$: _Tree techniques_.

#full[
=== More patterns
- Insertion (component) DP for permutations: insert values in
  increasing order and track the number of blocks $j$ that the
  final order will glue together. A new value opens a block
  ($j + 1$ places), sticks to an end of a block ($2 j$ places)
  or joins two neighbouring blocks ($j - 1$ places). Counts
  permutations with given up/down patterns, or costs that depend
  on whether a neighbour is larger.
- Sorting first fixes the order of a DP: sort by deadline and
  knapsack over time ("take a job only before its deadline"), by
  $w + s$ for stacking boxes of weight $w$ and strength $s$, by
  increasing $v$ when only the largest triggered $v$ counts
  ($f <- p v + (1 - p) f$).
- Huge `n` with a fixed linear transition: matrix power
  (`linalg/matrix`), over (min, +) or (max, +) too. Linear
  recurrences: `poly/linear-rec`.
- Monotone `opt` for a single array ($"dp"$ depends on itself):
  keep a deque of (candidate $j$, first $i$ where $j$ is best);
  a new $j$ pops candidates from the back while it beats them at
  their first $i$, then binary searches where it starts to win.
  $O(n log n)$.
- Keeping only the best value per state is exact when every
  transition is non-decreasing in that value. Counting optimal
  plans as (best, count) pairs needs each transition strictly
  increasing in the value (only best predecessors reach the best)
  or independent of it (every plan ties: keep the total count
  too).
- Counting inputs that a deterministic check accepts (a greedy
  with a few counters, a small DP): the check's state is the DP
  state ("DP of DP"). Prove the greedy exact first.
- Counting distinct results, not ways: build each result through
  one canonical decomposition (the leftmost cut, the shortest last
  block) and forbid or subtract the choices that rebuild a result
  already counted.
- A product $product_v c_v$ inside the sum counts the ways to mark
  one item in every box ($c_v$ items in box $v$): count (object,
  marking) pairs instead, with the number of boxes still waiting
  for their mark in the state.
- Lexicographically smallest optimum, or rebuilding without
  parent pointers: run the DP backwards (best completion, or "can
  finish", from each state), then walk forwards taking the
  smallest transition that keeps the optimum.
- Matrix chain order for dimensions $w_0, dots, w_n$ is the
  cheapest triangulation of the polygon with vertex weights $w_i$,
  a triangle costing the product of its weights. Hu-Shing: some
  optimal triangulation joins the smallest weight to the second
  and third smallest (by sides or diagonals).
- Huffman in $O(n)$ after sorting: merged weights come out
  non-decreasing, so keep two FIFO queues (sorted leaves, merged
  nodes) and take the two smallest fronts. A group of $c$ equal
  minimum weights $v$ merges at once into $floor(c\/2)$ nodes of
  weight $2 v$, keeping the odd one.
- L1 isotonic regression, $min sum |a_i - b_i|$ over
  non-decreasing $b$ (`dp/slope-trick`), equals the largest
  $sum (a_i - a_j)$ over sets of disjoint pairs $i < j$.
- All $F(k)$ of a concave $F$ by the Aliens trick: non-negative
  integer slopes with total at most $S$ take $O(sqrt(S))$ distinct
  values. Bisect $lambda$ ranges, stop when both ends report the
  same count, and fill $F$ linearly between consecutive vertices:
  $O(sqrt(S) log S)$ solver calls.
- Weighted squared distances: for $u_j <= x$,
  $min_j m_j (x - u_j)^2 = (min_j sqrt(m_j) (x - u_j))^2$, a
  minimum of lines of slope $sqrt(m_j)$ (`ds/li-chao`, adding line
  $j$ once $x$ passes $u_j$); mirror for $u_j >= x$. Recompute
  the winner's value exactly.
- Huge counts of a few kinds of items: an optimum is mostly copies
  of one or two patterns, adjacent on the upper convex hull of
  value against composition. Fix all but a margin of a few copies
  greedily, solve the rest with an exact DP, and test the margin
  against a brute force.
- Min-plus power of a cost on weights $0..n$, $c$ items: some
  $floor(c\/2)$ of $c$ weights with total $W$ sum to within $n$ of
  $W floor(c\/2) \/ c$ (slide a window over them cyclically), so
  $F_c$ comes from $F_(floor(c\/2))$ and $F_(ceil(c\/2))$ on
  windows of width $O(n)$: $O(n^2 log c)$.
- Meet in the middle on coordinates: to count pairs where an
  earlier vector dominates a later one in every coordinate, when
  each half of the coordinates takes few distinct patterns, keep
  `cnt[h2][h1]`, the items whose second half is `h2` and whose
  first half dominates `h1`. An item asks the sum of
  `cnt[h2'][h1]` (its own `h1`) over the `h2'` that dominate its
  second half, then adds 1 to `cnt[h2][h1']` (its own `h2`) for
  every `h1'` that its first half dominates.
- LIS of a concatenation $A B$: the maximum over $x$ of the number
  of `tail` values of $A$ below $x$ plus the number of `head`
  values of $B$ that are $>= x$ (`head[k]`: the largest first
  value of an increasing subsequence of length $k + 1$). With
  distinct values the pairs (longest increasing, longest
  decreasing subsequence ending at $i$) are pairwise distinct, so
  more than $r s$ values contain an increasing subsequence of
  length $r + 1$ or a decreasing one of length $s + 1$
  (Erdős-Szekeres).
]
