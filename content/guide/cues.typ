#import "/typst/notes.typ": *
#let xor = sym.plus.o
#let cues(..cells) = table(
  columns: (1.1fr, 1fr),
  stroke: none,
  inset: (x: 1.5pt, y: 1.3pt),
  fill: (_, y) => if calc.odd(y) { luma(236) },
  ..cells,
)

Left: what the statement says. Right: what to think of first, with
the snowi item. The notes named under each heading hold the details.

=== Counting
_Counting formulas_, _Generating functions_, _Polynomial
techniques_.
#cues(
  ["count ...", $n <= 5000$], [DP in a canonical build order],
  ["for every $k = 1, dots, n$"], [convolution: `poly/ntt`],
  [binomials mixing $i$ and $k$], [split into $i$ and $k - i$ parts],
  ["every value used", "no bad position"], [inclusion-exclusion],
  ["exactly $k$" fixed points or pairs], ["at least", binomial inversion],
  [easy if items could collide], [Stirling inversion],
  [labelled, "connected", "components"], [EGF: all $= exp$(connected)],
  [a root with children of the same kind], [Lagrange inversion],
  [steps $>= -1$, prefix never below 0], [cycle lemma, reflection],
  [product of many small polynomials], [divide and conquer + NTT],
  [$f_n$ from a convolution of earlier $f$], [online (CDQ) convolution],
  [huge $m$ ($10^9$), only comparisons], [polynomial: `poly/lagrange`],
  [fixed linear step, $n <= 10^18$], [`linalg/matrix`, `poly/linear-rec`],
  ["up to rotation", distinct colourings], [Burnside],
  [non-crossing lattice paths], [LGV determinant],
  [sum $equiv r$ mod a small $k$], [DP over residues],
  [ups and downs of a permutation], [insertion (component) DP],
  [digit conditions on $[L, R]$], [digit DP, $f(R) - f(L - 1)$],
  [convolution mod $10^9 + 7$], [`poly/fft-mod`],
)

=== Probability
_Probability and expectation_.
#cues(
  ["expected number of ..."], [sum of probabilities per item],
  [expected maximum or time], [$E[X] = sum_(k >= 1) P(X >= k)$],
  [random process until it stops], [DP backwards from the end],
  ["output $p q^(-1)$ modulo ..."], [count outcomes, divide once],
  [random removals, shuffled deck], [a uniform random permutation],
  ["may stop at any time"], [optimal stopping, a threshold],
  [expected time of the last event], [min-max inclusion-exclusion],
)

=== Number theory
_Number theory facts_.
#cues(
  [gcd, coprime, "lowest terms" in a sum], [Möbius over divisors],
  [$floor(n\/i)$, multiples up to $n$], [$O(sqrt(n))$ floor blocks],
  [$n approx 10^(11)$, one answer], [`numtheory/prime-count`],
  [coprime $x + y = s$], [$phi(s)$ choices of $x$],
  [binomials mod a small prime], [Lucas: digits in base $p$],
  [composite modulus], [CRT: `numtheory/euclid`],
  [$a^x equiv b$ modulo $m$], [`numtheory/discrete-log`],
  [$x^2 equiv a$ modulo a prime], [`numtheory/mod-sqrt`],
  [product is a perfect square], [exponent parities as XOR vectors],
  [$sum_i floor((a i + b) \/ m)$], [`numtheory/floor-sum`],
  [$f_(i + 1) = a_i f_i + f_(i - 1)$], [continued fractions, Euclid],
  ["is a square modulo $p$"], [Legendre symbol: one bit],
  [power towers, huge exponents], [exponent mod $phi(m)$, $+ phi(m)$],
  [repeated $x mod a_i$], [only $O(log x)$ steps change $x$],
)

=== Bits and subsets
_Bit tricks_, _Polynomial techniques_.
#cues(
  [sum over submasks, for every mask], [SOS: `poly/subset-transforms`],
  [pairs whose OR, AND or XOR is $x$], [transform, multiply, invert],
  [split a set into disjoint parts], [subset convolution],
  [same XOR move repeated $k$ times], [transform, power pointwise],
  [only how many bits flip matters], [collapse to Hamming weights],
  [wildcards `?` in bit patterns], [ternary transform, $O(3^n)$],
  [subsets of $n <= 40$ items], [meet in the middle],
)

=== XOR and linear algebra
_Linear algebra facts_.
#cues(
  ["XOR of a subset": max, count], [`linalg/xor-basis`, $2^(n - r)$],
  [XOR of a walk, edges reusable], [cycle basis + root XORs],
  [span or max XOR of $a_l, dots, a_r$], [prefix basis: `linalg/xor-basis`],
  [toggles, Lights Out], [system over $FF_2$: `linalg/gauss-mod`],
  [count solutions of parity rules], [$2^f$, $f$ free variables],
  [only "equal" and "different" rules], [parity DSU, 2-colouring],
  [count spanning trees], [matrix-tree determinant],
  [does a perfect matching exist?], [random Tutte matrix],
)

=== Strings
_String facts_.
#cues(
  [prefix equals suffix (borders)], [`strings/prefix-function`],
  ["made of copies", periods], [period $n - pi_(n - 1)$, Fine-Wilf],
  [two palindromic prefixes], [a period; reflect $x -> 2 y - x$],
  [smallest rotation, suffix or cut], [Lyndon: `strings/lyndon`],
  [many patterns in one text], [`strings/aho-corasick`],
  [distinct substrings, $k$-th substring], [`strings/suffix-automaton`],
  [compare substrings, LCP], [`strings/string-hash`, suffix array],
  [longest palindrome at each centre], [`strings/manacher`],
  [distinct palindromic substrings], [`strings/eertree`],
  [matching with wildcards], [FFT, one product per letter],
  [lengths $10^18$ from repeated prefixes], [Euclid-like shrinking],
)

=== Geometry
_Geometry facts_.
#cues(
  [max of $a_k x + b_k y$ per query], [hull, tangent binary search],
  [sums of two choices, collisions], [`geometry/minkowski`],
  [farthest pair, width], [rotating calipers],
  [best window as (length, sum)], [D&C, Minkowski sums of hulls],
  [mixtures with weights summing to 1], [convex hull of the options],
  [region inside many lines], [`geometry/half-plane`],
  [exact answers, coordinates $10^9$], [integer cross products],
  [one parameter added to everything], [each candidate is a line],
  [Manhattan distance], [rotate: $(x + y, x - y)$],
)

=== DP speed-ups
_DP techniques_.
#cues(
  [$(h_i - h_j)^2$ or $a_i b_j$ in a transition], [`ds/li-chao`, line container],
  ["exactly $k$ pieces", "at most $k$"], [`dp/wqs` if the cost is convex in $k$],
  [cost with the quadrangle inequality], [`dp/dc-opt`, `dp/knuth`],
  [costs $|x - a|$, "make it sorted"], [`dp/slope-trick`],
  [$sum_i b_i = m$, convex costs $f_i (b_i)$], [the $m$ cheapest unit steps],
  [maximise $sum a \/ sum b$], [binary search $lambda$],
  [best with $c$ units from two parts], [concave: merge the slopes],
  [minimum over a sliding window], [monotone deque],
  [longest chain, increasing subsequence], [`dp/lis`; Dilworth],
)

=== Ranges and updates
_Designing a segment tree node_, _Offline techniques_, _Sqrt
decomposition and small tricks_.
#cues(
  [range answer by a left-to-right scan], [node = composed function],
  ["after each update, print the DP"], [(min, +) matrices, `ds/segtree`],
  ["first position where ..."], [descend the segment tree],
  [range chmin or chmax], [`ds2/segtree-beats`],
  [distinct values in $[l, r]$], [offline by $r$, last occurrence],
  [$k$-th smallest in $[l, r]$], [`ds2/wavelet-matrix`],
  [$a_j <= a_i$, $b_j <= b_i$, $c_j <= c_i$], [`offline/cdq`],
  [edges added and removed], [`offline/dynamic-connectivity`],
  ["the first time when ..." per query], [`offline/parallel-binary-search`],
  [static ranges, $O(1)$ add and remove], [`offline/mo`],
  [a cyclic sequence of transitions], [trace of the matrix product],
  [next greater, span of a maximum], [monotone stack],
  [sum of the max of all subarrays], [contribution of each maximum],
  ["at least one 0 in $[l, r]$"], [complement: runs of ones],
  [order statistics, online], [`ds2/ordered-set`, `ds2/treap`],
  [$max_j (a_j xor x)$], [`ds2/binary-trie`],
)

=== Trees
_Tree techniques_.
#cues(
  [per subtree, a multiset statistic], [small-to-large],
  ["for every vertex as the root"], [`tree/rerooting`],
  [count paths with a property], [`tree/centroid`],
  [nearest marked vertex, updates], [centroid ancestors],
  [queries on $k$ vertices, $sum k$ bounded], [`tree/virtual-tree`],
  ["only edges of weight $<= x$"], [Kruskal reconstruction tree],
  [path updates and path queries], [`tree/hld`],
  [walk up, pay at each new maximum], [next greater ancestor],
  [link and cut edges online], [`tree/link-cut`],
  [subtree sums with updates], [Euler tour + `ds/fenwick`],
)

=== Graphs and flows
_Graph theory facts_, _Flow and cut modelling_, _Matroids_.
#cues(
  [two options per object, pair rules], [`graph/two-sat`],
  ["at most one per group" in 2-SAT], [prefix variables],
  [fewest deletions to separate], [min cut: `flow/dinic`],
  [taking $x$ forces $y$, max profit], [project selection],
  [assign with capacities], [bipartite flow, Hall],
  ["is it always possible?", $n = 10^6$], [cut condition, then greedy],
  [most items, no two conflict, bipartite], [$|V| - nu$ (Kőnig)],
  [cheapest assignment or transport], [`flow/hungarian`, min-cost flow],
  [cover a DAG with disjoint paths], [$n - nu$],
  [use every edge exactly once], [`graph/euler-path`],
  [no cycle / independent / matchable], [matroid greedy],
  [two such rules at once], [matroid intersection],
  ["differ by exactly 1 on each edge"], [bipartite: parity = colour],
  [state = vertex plus extra info], [BFS or Dijkstra on pairs],
  [edges $(u + k, v + k)$ for all $k < L$], [DSU per power of two],
  [fewest swaps to sort], [$n - "cycles"$; adjacent: inversions],
  [bridges, cut vertices], [`graph/bcc`],
  [reachability, cycles in a digraph], [`graph/scc`, condensation],
  [orient $K_n$, given out-degrees], [Landau, greedy build],
)

=== Games and interaction
_Game theory_.
#cues(
  [sum of independent games], [Grundy values, XOR],
  [huge numbers, simple moves], [brute force table, period],
  [token game, visited vertices deleted], [wins iff start in all max matchings],
  [symmetric position], [mirror or pairing strategy],
  [different moves per player], [win/lose DP, no Grundy],
  [players collect points], [DP on my score minus yours],
  [query limit near $log_2$(candidates)], [halving: `misc/interactive`],
  [hidden vertex in a tree], [centroid search],
)

=== Randomisation
_Randomisation and hashing_.
#cues(
  [all or none per value, even counts], [random XOR labels],
  [equal multisets, "is a permutation"], [sums of random labels],
  [if the answer is no, many witnesses], [sample random candidates],
  [count equal subarrays or paths], [prefix fingerprints],
  [greedy fast for a random order], [shuffle first],
)

=== Constructions
_Constructive problems and invariants_.
#cues(
  ["print any, or $-1$"], [necessary condition, brute force],
  [reach B from A by moves], [invariant: parity, sum mod $k$],
  [the input is only $n$ (and $k$)], [brute force $n <= 8$, pattern],
  [budget of $floor(k\/3)$ changes], [3 classes, take the smallest],
  [regular: "each in exactly 3"], [Hall, perfect matchings],
  [a process that repairs violations], [decreasing potential],
  [the answer always seems small], [prove $<= c$, then the rest],
  ["your program runs twice"], [encode a unique structure],
)
