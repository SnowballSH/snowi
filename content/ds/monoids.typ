#import "/typst/notes.typ": *

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

=== Standard nodes
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
#fact[Runs][Store `(first, last, runs)`; merging subtracts one
  run when `a.last == b.first`. The identity needs a flag
  (`runs = 0`).]
#fact[Hash][Store `(h, p)` with $p = B^"len"$; merge to
  `(a.h * b.p + b.h, a.p * b.p)` mod a prime, identity `(0, 1)`.
  Keep a reversed hash too for palindrome checks under point
  updates.]
#fact[First index]["First `i >= l` with `a[i] >= x`" is
  `max_right(l, v < x)` on max (`n` if none). On sums of values
  $>= 0$, `r = max_right(l, s <= k)` makes `[l, r + 1)` the
  shortest range from `l` with sum `> k` (none if `r = n`).]

=== Lazy tags
A tag `f` needs three laws: it acts on a node summary in $O(1)$
(`mapping`), two tags fold into one (`composition(new, old)`,
new on the outside), and it distributes over `op`.
#fact[Add][On min or max: `s + f`. On sums the node also
  stores `len`: `sum + f * len`.]
#fact[Assign][Tag `f` or "none"; `composition` keeps the newer
  one. On sums: `f * len`. With adds too, the tag is
  `(assign or none, add)` and an assign wipes earlier adds.]
#fact[Affine][On `(sum, len)`: `(b sum + c len, len)`; compose as
  above.]
#fact[Negate][On `(sum, mn, mx)`: `(-sum, -mx, -mn)`.]
#fact[Flip bits][On `(ones, len)`: `(len - ones, len)`.]

=== When it does not fit
- The tag cannot update a summary without looking inside: range
  chmin with range sum. Use segment tree beats: store `mx`, the
  strict second max `smx` and `cnt` of `mx`; stop if `mx <= x`,
  tag if `smx < x < mx`, else recurse. $O((n + q) log n)$
  amortized, $O(log^2 n)$ with adds.
- Range `a_i %= x`, `a_i /= x`, `a_i = floor(sqrt(a_i))` on values
  $>= 0$: store `mx` and recurse only into nodes that change:
  `mx >= x` for `%=`, `mx >= 1` for `/=` (skip `x = 1`), `mx >= 2`
  for sqrt. A changed element at least halves: $O(log V)$ changes.
- k-th smallest in a range, count distinct: merge sort tree,
  wavelet tree, persistent tree, or offline sweeps.
- `op` is too heavy (sets, sorted lists): sqrt decomposition or Mo.

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

=== More patterns
- Visible buildings (count prefix maxima of a range): a node
  answers "how many prefix maxima above `h`" by recursing into one
  child, $O(log^2 n)$ per operation.
- "Count positions with value `> 0`" under range add when the
  global minimum is known to be 0: compare with the node minimum
  instead, a predicate that survives lazy adds.
- "Best total with exactly $k$ items" on a range, $g(k)$ concave:
  a node stores the whole sequence $g$ and merges by (max, +)
  convolution, which for concave sequences merges the sorted
  slopes. Or use the Aliens trick instead: charge $lambda$ per
  item, keep ordinary nodes (plus the item count) for that fixed
  $lambda$, and binary search $lambda$.
- Cues: "after each update, print the answer", "first position
  where ...", "the sign of a huge sum" (saturating maps),
  "minimum deletions so that ..." on substrings (automaton DP as
  min-plus matrices).
]
