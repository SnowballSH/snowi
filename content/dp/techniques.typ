#import "/typst/notes.typ": *

Write the transition, then match its shape below. Check a claimed
property (monotone `opt`, convexity in $k$) by brute force first.

=== Which speed-up applies
#fact[Quadrangle inequality][$w(a, c) + w(b, d) <= w(a, d) + w(b, c)$
  for $a <= b <= c <= d$: crossing intervals cost at most nested
  ones. Holds for $f(S_r - S_l)$ with $f$ convex and $S$
  non-decreasing, and for sums over the pairs inside $[l, r)$ of
  weights $>= 0$; it survives sums and adding $g(l) + h(r)$.
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
#fact[Exactly $k$][The best value is convex (minimum) or concave
  (maximum) in $k$: drop $k$, charge $lambda$ per item, binary
  search $lambda$ (`dp/wqs`). Convex: min-cost flow in the amount,
  matroids, partitions with a quadrangle-inequality cost.]
#fact[Convex in a value][A DP over one numeric state $x$ that stays
  convex piecewise linear in $x$ (costs $|x - a|$, prefix minima,
  windows): `dp/slope-trick`.]
#fact[Window][$"dp"[i] = c_i + min_(i - L <= j < i) "dp"[j]$: keep a
  deque of indices with increasing `dp`; pop the back while it is
  not better than the new one, the front once out of range.
  $O(n)$.]

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
  a sliding window maximum (Window, mirrored).
- 0/1 items: loop capacity downwards; unbounded: upwards. Huge
  weights, small values: `dp[v]` = least weight for value $v$.

=== Bitmasks and digits
- All submasks `s` of all masks `m`, $3^n$ in total: start at
  `s = m`, step `s = (s - 1) & m`, stop after `s = 0`.
- Sum over submasks in $O(n 2^n)$: for each bit `b`, for each
  mask `m` with bit `b`: `f[m] += f[m ^ 1 << b]`. Supersets,
  inverses, subset convolution: `poly/subset-transforms`.
- Paths through all vertices: $"dp"["mask"][v]$, masks in
  increasing order, $O(2^n n^2)$. Assigning $n <= 20$ items:
  `dp[mask]`, the next item is `popcount(mask)`.
- Broken profile (tilings of $n times m$, $m$ small): cells in
  row-major order, the mask holds the next $m$ cells (covered or
  not): $O(n m 2^m)$.
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
- Tree knapsack: merging child tables of sizes $s_1$ and $s_2$
  costs $s_1 s_2$, which sums to $O(n^2)$ (each pair of vertices
  meets once, at its LCA). Cap tables at $K$ entries for $O(n K)$.
  Size the merged table $min(s_1 + s_2, K + 1)$ and loop only over
  existing entries.

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
  $w + s$ for stacking boxes of weight $w$ and strength $s$.
- Rolling arrays: keep two layers (or one, with the right loop
  direction). With `-INF` for unreachable states, never add to it
  (`LLONG_MIN + x` overflows); use `LLONG_MIN / 4` or skip.
- Huge `n` with a fixed linear transition: matrix power
  (`linalg/matrix`), over (min, +) or (max, +) too. Linear
  recurrences: `poly/linear-rec`.
- Probability and expectation DPs usually run backwards from the
  final states; see the probability note.
- Monotone `opt` for a single array ($"dp"$ depends on itself):
  keep a deque of (candidate $j$, first $i$ where $j$ is best);
  a new $j$ pops candidates from the back while it beats them at
  their first $i$, then binary searches where it starts to win.
  $O(n log n)$.

=== Example: digit sums
Count $x <= N$ with digit sum divisible by $D$: states $q$ = sum
mod $D$, `nxt(q, d) = (q + d) % D`, `accept(q) = q == 0`. Count
in $[L, R]$ as `count(R) - count(L - 1)` (handle $L = 0$).
]
