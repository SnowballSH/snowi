#import "/typst/notes.typ": *

Read the limits, decide what is asked, check the cost of the
candidate tools, and rewrite the problem until one fits.

=== Step 1: budget from the limits
#table(
  columns: (auto, auto, 1fr),
  [*$n$ up to*], [*time*], [*typical tools*],
  [10–12], [$n!$, $n! dot n$], [all orders, backtracking],
  [20–25], [$2^n n$], [bitmask DP; to 40: meet in the middle],
  [500], [$n^3$], [interval DP, Floyd-Warshall, Hungarian],
  [5000], [$n^2$], [all pairs, 2D DP],
  [$10^5$–$10^6$], [$n log n$, $n sqrt(n)$], [sort, trees, sweeps, Mo],
  [$10^(12)$], [$sqrt(n)$], [divisors, floor blocks],
  [$10^18$], [$log n$, $O(1)$], [formula, matrix power, binary search],
)
About $10^8$ simple operations per second in #box[C++], 3 to 10 times
fewer in PyPy. Count the queries and the sum of $n$ over tests. Small
values or a small $k$ are budgets too (index by value, $2^k$).

=== Step 2: the lane
#fact[A. Count][Ways, probabilities, expectations mod $p$. Sum
  of contributions of single elements or pairs? Complement or
  inclusion-exclusion? A canonical order that builds each object
  once? A huge parameter (polynomial in it)? All $k$ at once
  (convolution)? DP, `combo/binom`, _Counting formulas_,
  _Probability and expectation_, `poly/ntt`, `poly/lagrange`,
  `linalg/matrix`.]
#fact[B. Optimise][Monotone in the answer (binary search)? An
  exchange argument that fixes the order (sort, then DP)? The
  shape of the DP transition? Exactly $k$ items (WQS)? A graph:
  shortest path, MST, flow, cut? _DP techniques_, `dp/dc-opt`,
  `dp/wqs`, `ds/li-chao`, `graph/dijkstra`, `ds/dsu`,
  `flow/min-cost-flow`, `misc/ternary-search`.]
#fact[C. Decide or build][An invariant or colouring that rules out
  the bad cases? Does the brute force agree on all small inputs?
  Two options per object (2-SAT)? Matching and Hall, degrees and
  Euler? `graph/two-sat`, `flow/hopcroft-karp`, _Flow and cut
  modelling_, `graph/euler-path`, _Graph theory facts_.]
#fact[D. Many queries][Offline allowed (sort, sweep, reverse time,
  Mo)? Updates? What does a node store and how do two merge? Paths
  and subtrees: Euler tour, HLD. `ds/fenwick`, `ds/segtree`,
  `ds/lazy-segtree`, _Designing a segment tree node_,
  `ds/sparse-table`, `tree/lca`, `tree/hld`, `offline/mo`.]
#fact[E. Game][Impartial with independent parts (Grundy values,
  XOR)? Brute force small positions and find the pattern; try
  pairing and symmetry. _Game theory_.]

=== Step 3: check the fit
Count operations with the constant before coding; check the model
on the samples. Right model but too slow (a min cut on $10^6$
vertices, all $n^2$ pairs): keep it, cut its cost with step 4.

=== Step 4: rewrites
- Complement: "at least one" is all minus none; "some has $P$" is
  "not all fail $P$" ("a 0 in $[l, r]$" is "not all ones").
- Contribution: how often each element or pair counts (linearity);
  a DP can carry (count, sum) pairs per state.
- Think in runs and intervals: maximal equal blocks, gaps, cuts.
  A reachable or feasible set that is an interval is two numbers,
  and winning states that form a prefix are one threshold (removing
  a block from a scan: _DP techniques_).
- Split: rules that separate by coordinate or by bit are solved
  one coordinate or bit at a time (counts multiply, costs add).
- Permutations: cycles; fewest swaps to sort is $n - "cycles"$.
- Parities and invariants: what never changes under the moves?
- Sort and sweep by value, time or coordinate, with events.
- Binary search the answer: "least maximum" becomes yes or no.
- Fix one dimension (the maximum, the first element), solve the rest.
- Offline: reorder queries; reverse time (deletions become inserts).
- Backwards: run the process from the end, where choices are forced.
  To remove one block so a scan still succeeds, compare the state
  after each prefix with the largest state that succeeds from each
  suffix (_DP techniques_).
- Small-to-large: merge the smaller set into the larger.
- Pattern: print brute force answers for small inputs, guess, prove.
- Put the update on the edge: store differences (between neighbours,
  child minus parent), so a range update touches two places and a
  suffix or subtree update one. Adding a linear function on a
  range: second differences, prefix-summed twice.
- Reduce to a graph: states as vertices, moves as edges.
- Polynomial in a parameter: sample $d + 1$ values, `poly/lagrange`.

#full[
=== More on the budget
- $n <= 10^7$: $O(n)$ only (linear sieve, two pointers); an `int`
  array of $10^7$ is 40 MB.
- $n approx 10^(11)$ with one query: sublinear prefix sums, such as
  $O(n^(3\/4) \/ log n)$ prime counting (`numtheory/prime-count`)
  or $O(n^(2\/3))$ with Du's sieve.
- $n = 2 dot 10^5$ in 2 s: $n log^2 n approx 6 dot 10^7$ is safe,
  $n sqrt(n) approx 9 dot 10^7$ is fine in #box[C++] with a small
  constant and risky in PyPy.
- Two sizes $n m <= 10^7$ suggest an $O(n m)$ DP; one dimension
  $<= 20$ suggests going exponential in it.
- Memory: 256 MB holds about $6 dot 10^7$ `int`s in #box[C++] and
  about $3 dot 10^7$ small ints in PyPy lists (8 bytes each).
- The output is a hint: mod 998244353 suggests convolutions; a
  real number with error $10^(-6)$ suggests probability, geometry
  or a binary search on reals; "any valid answer" means a
  construction.
- Online queries (decoded with the last answer) rule out the
  offline tools of lane D. Interactive: compare the query limit
  with $log_2$ of the number of candidates (`misc/interactive`).

=== More questions per lane
- A: a sum over a structure (subsets, paths, trees) of a product:
  DP over that structure. Symmetric objects: Burnside. Labelled
  components: exp and log of generating functions (_Generating
  functions_). Gcd and divisors: Möbius (_Number theory facts_).
- B: the best bottleneck (largest minimum edge): sort the edges
  and use DSU, or binary search. Convex in one variable: ternary
  search or `dp/slope-trick`. A ratio: binary search $lambda$.
  Feasibility rule with exchange: _Matroids_. Few distinct values:
  enumerate the value.
- C: "fewest operations": BFS over the states for small sizes,
  then find the formula. "Print any": _Constructive problems and
  invariants_.
- D: no updates: prefix sums, `ds/sparse-table`, offline sweeps.
  Point updates: `ds/fenwick`, `ds/segtree`. Range updates:
  `ds/lazy-segtree`. Order statistics in a range:
  `ds2/wavelet-matrix`. Edges that come and go:
  `offline/dynamic-connectivity`. See _Offline techniques_.
- E: different moves for the two players: plain win/lose DP, no
  Grundy values. Cycles: retrograde analysis with draws.
  Randomness: expectations over the chance moves.

=== Worked lanes (UCup S5 Extra Stage 1)
- "At least one 0 in $[l, r]$", min cut correct but too slow
  (B): complement to "all ones" and think in runs; DP over the
  runs with a range-add, range-max segment tree.
- "Orient $K_n$ with out-degree minus in-degree divisible by $d$"
  (C): the out-degrees form a score list; Landau's test and the
  greedy build (_Constructive problems and invariants_).
- "Count arrays, values up to $m = 10^9$", only comparisons
  matter (A): the answer is a polynomial in $m$; DP for small $m$,
  then `poly/lagrange`.
- Components of "edge between $i < j$ iff $p_i > p_j$" (D): they
  are intervals, cut after $k$ iff $p_0, dots, p_k$ are the
  smallest $k + 1$ values; count cut points with a sweep and
  `ds/fenwick`.
- "Fewest swaps" across permutations (B): the mismatches are $n$
  minus the fixed points of one permutation $q$, and a swap splits
  or merges one cycle of $q$.
- "$x$ mod $a_i$ along a path" (D): only $O(log x)$ of the mods
  change $x$, since each one that does at least halves it.

=== When stuck
- Re-read for an unused constraint: distinct values, connected,
  $n$ even, the sum of $n$ bounded. An odd constraint is a hint.
- Solve the line, tree or small-$k$ version first.
- Look from the other side: the dual (cut instead of flow,
  potentials instead of paths), the reverse process, the
  complement graph, coordinates rotated by 45 degrees.
- Iterate over the answer value instead of the input, over
  divisors instead of multiples, or over the smaller side (split
  by size at $sqrt(n)$).
- Randomise: a random order, random weights, random sampling when
  witnesses are many (_Randomisation and hashing_).
]
