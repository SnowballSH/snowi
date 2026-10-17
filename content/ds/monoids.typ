#import "/typst/notes.typ": *

#full[
A segment tree works for any *monoid*: an associative `op` with an
identity `e`. It need not be commutative, since the tree always
combines left before right. All the thinking goes into the node.

=== Recipe
+ Start from the answer for one range.
+ Cut the range into a left part `a` and a right part `b`. Ask what
  goes wrong when the answer crosses the cut, and add exactly the
  fields that case needs.
+ Repeat for the new fields until they merge on their own.
+ Define every field by its meaning on the segment, write `e`
  explicitly, and stress-test `op` on short arrays against brute
  force.

A range whose answer comes from a left-to-right scan with a small
state is a composition of functions (one per element): store the
composed function if it stays in a small family.
]

=== Standard nodes
#full[
#fact[Max subarray sum][Store `(sum, pre, suf, best)`; the best
  crossing subarray is `a.suf + b.pre`. Identity
  `(0, -INF, -INF, -INF)`. If the empty subarray counts, leaves
  are `(v, m, m, m)` with `m = max(v, 0)` and `e` is all zeros.]
#lang-code(
  cpp: "struct T { ll sum, pre, suf, best; };
T op(T a, T b) {
  return {a.sum + b.sum, max(a.pre, a.sum + b.pre),
          max(b.suf, b.sum + a.suf),
          max({a.best, b.best, a.suf + b.pre})};
}",
  py: "def op(a, b):  # (sum, pre, suf, best)
  return (a[0] + b[0], max(a[1], a[0] + b[1]),
          max(b[2], b[0] + a[2]), max(a[3], b[3], a[2] + b[1]))",
)
#fact[Min and its count][Store `(mn, cnt)`: keep the smaller
  `mn`, add counts on ties. Under range add, the number of zeros is `cnt`
  when `mn = 0` and values stay $>= 0$: area of a union of
  rectangles, positions covered by no interval.]
]
#fact[Affine maps][Store $x -> b x + c$ as `(b, c)`. "`f`, then
  `g`" is `(g.b f.b, g.b f.c + g.c)`, identity `(1, 0)`.]
#fact[Clamped shift][Store $x -> min(h, max(ell, x + a))$ as
  `(a, lo, hi)`: any sequence of "add", "chmin", "chmax".
  "`f`, then `g`" is `(f.a + g.a, g(f.lo), g(f.hi))`, identity
  `(0, -INF, INF)`.]
#fact[Dynamic DP][A DP with $k$ states and transition
  $"dp"_i [t] = max_s ("dp"_(i-1) [s] + M_i [s][t])$ is a product of
  $k times k$ matrices over (max, +). Identity: 0 on the diagonal,
  $-infinity$ elsewhere. Rows are the state before, columns the
  state after. Independent set on a path: $M(w) = mat(0, w; 0,
  -infinity)$. On a cycle, take the trace: $max_s P[s][s]$.]
#fact[First index][On sums of values $>= 0$,
  `r = max_right(l, s <= k)` makes `[l, r + 1)` the shortest range
  from `l` with sum `> k` (none if `r = n`). Over compressed values
  with leaves `(count, sum)`, `min_left(V, c <= t)` finds where the
  `t` largest start: their sum is that range's sum plus the missing
  copies of the value just below.]
#fact[Prefix maxima][Records of a range, or
  $sum_i max(a_l, dots, a_i)$, under point updates. Node `v` keeps
  its max and `R[v] = f(right, max(left))`, where `f(u, h)` is the
  answer for `u` entered with running maximum `h`. `f` follows one
  child: if `h >= max(left)`, the left part is all `h` (no
  records) and recurse right; else return `f(left, h) + R[v]`.
  So `f` is $O(log n)$, a pull-up calls it once ($O(log^2 n)$ per
  update), and a query carries `h` through its $O(log n)$ nodes
  from left to right. A running min and max together (sums of
  $min dot max$) work the same, with a one-sided `f` where only
  one cap still binds: $O(log^2 n)$ per call.]
#fact[Another colour][Best value among elements whose colour is
  not `c` (Borůvka: the cheapest edge leaving my component). Store
  `(v1, c1, v2)`: the best value, its colour, and the best value of
  any other colour. Merge: `v1` is the larger best; `v2` is the
  max of both `v2` and of each child's `v1` whose colour differs
  from the new `c1`. Answer `c1 != c ? v1 : v2`. Or sort the
  elements by colour and query prefix and suffix structures
  around the block of `c`.]

=== Lazy tags
#fact[Linear tags][If an update applies the same linear map to
  each element's vector $(x_1, dots, x_k, 1)$, the node stores the
  sum of each coordinate (the length for the 1), a tag is a
  $(k + 1) times (k + 1)$ matrix, and `composition` is the product
  "new times old". This covers add, assign, affine, adding
  $c + d i$ (keep $sum i$), power sums (on $(x, x^2, 1)$, adding
  $d$ maps $x^2$ to $x^2 + 2 d x + d^2$) and historic sums: with
  $(A, C)$, "add $C$ into $A$" each step makes $sum A$ the total of
  $C$ over time. Store only the entries that can be non-zero.]
#fact[Sorted arrays][If the array stays non-decreasing and an
  update applies a non-decreasing map made of "unchanged",
  "constant" and "shift" pieces (chmax with `c` on a suffix,
  $K -> K$ for $K < a$, $a - 1$ for $a <= K < 2 a$, $K - a$
  above), each piece is a contiguous range of positions: find the
  boundaries by descent, then range assign or add. The array stays
  sorted.]

=== When it does not fit
#py[
- The tag cannot update a summary without looking inside: range
  chmin with range sum. Use segment tree beats: store `mx`, the
  strict second max `smx` and `cnt` of `mx`; stop if `mx <= x`,
  tag if `smx < x < mx`, else recurse. $O((n + q) log n)$
  amortized, $O(log^2 n)$ with adds.
]
- Range `a_i %= x`, `a_i /= x`, `a_i = floor(sqrt(a_i))` on values
  $>= 0$: store `mx` and recurse only into nodes that change:
  `mx >= x` for `%=`, `mx >= 1` for `/=` (skip `x = 1`), `mx >= 2`
  for sqrt. A changed element at least halves: $O(log V)$ changes.
#full[
- k-th smallest in a range, count distinct: merge sort tree,
  wavelet tree, persistent tree, or offline sweeps.
]

#full[
=== Pitfalls
- Non-commutative `op`: the query keeps separate left and right
  accumulators, `b = op(t[--r], b)`. Swapping them passes tests
  with max or sum and fails with matrices.
- $-infinity$ arithmetic: `LLONG_MIN + x` overflows. Use
  `-INF = LLONG_MIN / 4` and clamp with `max(-INF, ...)`, or `int`
  with `min(INF, ...)`.
- The identity really must be neutral: the all-zero matrix is not
  the (max, +) identity, and padding leaves hold `e`.
- In C++ a member function named `max` hides `std::max` inside the
  struct; name queries `qmax`, `qsum`.
- With many test cases, build each tree with that test's `n`.
- Composed maps grow fast (affine, $2^k$ coefficients): reduce
  modulo a prime or saturate at a threshold.
- A negate tag swaps the extremes: `(sum, mn, mx)` becomes
  `(-sum, -mx, -mn)`.
- Assign and add in one tag: keep `(assign or none, add)`; an
  assign wipes the earlier adds, and a later add only grows the
  add part.

=== More patterns
- "Count positions with value `> 0`" under range add when the
  global minimum is known to be 0: compare with the node minimum
  instead, a predicate that survives lazy adds.
- "Best total with exactly $k$ items" on a range, $g(k)$ concave:
  a node stores the whole sequence $g$ and merges by (max, +)
  convolution, which for concave sequences merges the sorted
  slopes. Or use the Aliens trick instead: charge $lambda$ per
  item, keep ordinary nodes (plus the item count) for that fixed
  $lambda$, and binary search $lambda$.
- Components of the graph on `[l, r)` (cells of a $k times n$
  grid, edges to the nearest matching element): a node stores its
  component count and a label for each of the few vertices that
  can have edges leaving the range (the boundary columns, the
  first and last occurrence of each bit). Merge with a small DSU
  on the children's labels along the crossing edges; count = sum
  minus successful joins; then relabel the new boundary.
- Sign of $W = sum_(i = l)^r 2^(i - l) a_i$ with $|a_i| <= M$
  under point updates: evaluate $v <- 2 v + a_i$ from $i = r$ down
  to $l$, clamping $v$ to $[-C, C]$ with $C = M + 1$; once
  $|v| >= C$ the sign never changes. A node of length $k$ is the
  map $v -> "clamp"(2^k v + W)$ on $[-C, C]$. Store $W$ exactly
  while $2^k <= 2 C$; a longer node is a step, $C$ for $v > -A$
  and $-C$ for $v < -A$ with $A$ the integer nearest to $W \/ 2^k$:
  store $A$ and the value at $-A$.
- Unit jobs with deadlines (_Matroids_): a set fits iff for every
  $D$ at most $D$ of them are due by $D$. Keep $"slack"(D)$, that
  is $D$ minus the number of chosen jobs due by $D$, in a range
  add, range min tree; choosing a job due at $d$ adds $-1$ on
  `[d, n]`. A new job fits iff no tight $D >= d$ (slack 0);
  otherwise, with $D^*$ the first tight $D >= d$, the chosen jobs
  due by $D^*$ are the ones it can replace. After removing a
  chosen job, an outside job fits iff it is due after the last
  tight $D$.
- Cues: "after each update, print the answer", "first position
  where ...", "the sign of a huge sum" (saturating maps),
  "minimum deletions so that ..." on substrings (automaton DP as
  min-plus matrices), "components of a range" (boundary labels),
  "sum over all subarrays" under updates (historic sums).
]
