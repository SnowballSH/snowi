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
techniques_, _DP techniques_.
#cues(
  ["count ...", $n <= 5000$], [DP in a canonical build order],
  ["number of distinct results"], [one canonical build per result],
  [count inputs that a greedy accepts], [DP over the greedy's counters],
  [$product_v c_v$ inside the sum], [mark one item in every box],
  ["for every $k = 1, dots, n$"], [convolution: `poly/ntt`],
  [binomials mixing $i$ and $k$], [split into $i$ and $k - i$ parts],
  [sums of products of binomials], [Vandermonde, hockey stick],
  ["every value used", "no bad position"], [inclusion-exclusion],
  ["exactly $k$" fixed points or pairs], ["at least", binomial inversion],
  [easy if items could collide], [Stirling inversion],
  [labelled, "connected", "components"], [EGF: all $= exp$(connected)],
  [a root with children of the same kind], [Lagrange inversion],
  [prefix never below 0, brackets, trees], [reflection, cycle lemma, Catalan],
  [$plus.minus 1$ walk inside a strip, Kadane], [reflection sum; sticky walk],
  [product of many small polynomials], [divide and conquer + NTT],
  [paths under a monotone staircase], [D&C over columns, NTT],
  [$f_n$ from a convolution of earlier $f$], [online (CDQ) convolution],
  [huge $m$ ($10^9$), only comparisons], [polynomial: `poly/lagrange`],
  [a DP dimension that only adds up], [evaluate at $n + 1$ points],
  [fixed linear step, $n <= 10^18$], [`linalg/matrix`, `poly/linear-rec`],
  ["up to rotation", distinct colourings], [Burnside],
  [non-crossing lattice paths], [LGV determinant],
  [sum $equiv r$ mod a small $k$], [residue DP; roots of unity],
  [ups and downs of a permutation], [insertion (component) DP],
  [bounds $p_i <= b_i$, "$max_i i a_i <= k$"], [sort the bounds: position bounds],
  [prefix maxima, records], [insert values from the largest],
  [sum of a value over all choices], [DP with (count, sum) per state],
  [no three lengths form a triangle], [sorted, they grow like Fibonacci],
  [LIS or nesting over all permutations], [RSK, hook length formula],
  [digit conditions on $[L, R]$], [digit DP, $f(R) - f(L - 1)$],
  [convolution mod $10^9 + 7$], [`poly/fft-mod`],
)

=== Probability
_Probability and expectation_.
#cues(
  ["expected number of ..."], [sum of probabilities per item],
  [expected maximum or time], [$E[X] = sum_(k >= 1) P(X >= k)$],
  [random process until it stops], [DP backwards from the end],
  [random walk with cycles], [equations; tridiagonal on a line],
  ["output $p q^(-1)$ modulo ..."], [count outcomes, divide once],
  [random removals, shuffled deck], [a uniform random permutation],
  ["may stop at any time"], [optimal stopping, a threshold],
  [DP over $n$ random steps too wide], [keep a band of a few $sqrt(n)$],
  [several chains, one advances per step], [index policy (Gittins)],
  [expected time of the last event], [min-max inclusion-exclusion],
  [random shift, expected overlap], [Fubini: $integral "overlap" = |P| |Q|$],
)

=== Number theory
_Number theory facts_.
#cues(
  [gcd, coprime, "lowest terms" in a sum], [Möbius over divisors],
  [$floor(n\/i)$, multiples up to $n$], [floor blocks; multiples: $n ln n$],
  [$n approx 10^(11)$, one answer], [`numtheory/prime-count`],
  [coprime $x + y = s$], [$phi(s)$ choices of $x$],
  [$a x + b y = c$, $a x equiv c$ mod $m$], [Bézout: `numtheory/euclid`],
  [equations in $phi$, $d$, $sigma$], [prime by prime (Totients)],
  [binomials mod a small prime], [Lucas: digits in base $p$],
  [composite modulus], [CRT: `numtheory/euclid`],
  [$a^x equiv b$ modulo $m$], [`numtheory/discrete-log`],
  [powers modulo $2^k$], [units $plus.minus 5^j$, orders powers of 2],
  [$x^2 equiv a$ modulo a prime], [`numtheory/mod-sqrt`],
  [product is a perfect square], [exponent parities as XOR vectors],
  [$sum_i floor((a i + b) \/ m)$], [`numtheory/floor-sum`],
  [$f_(i + 1) = a_i f_i + f_(i - 1)$, resistors], [continued fractions, Euclid],
  [subtract the smaller, many rounds], [subtractive Euclid: one division],
  ["is a square modulo $p$"], [Legendre symbol: one bit],
  [power towers, huge exponents], [exponent mod $phi(m)$, $+ phi(m)$],
  [repeated $x mod a_i$, values that halve], [only $O(log x)$ steps change $x$],
  [which sums $sum c_i a_i$ are reachable], [Dijkstra on residues mod $min a_i$],
  [the answer divides $n$ or a gcd], [all divisors: $<= 1344$ up to $10^9$],
  [more subsets or prefixes than residues], [pigeonhole: two of them collide],
  [evenly spread, $floor(i q \/ n)$ ones], [balanced words],
  [digit sums, multiples of $b^s - 1$], [Digit sums: at least $s (b - 1)$],
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
  [AND, OR, XOR rules on numbers], [solve each bit on its own],
  [OR, AND or gcd of subarrays], [$O(log V)$ distinct values per end],
  [$min_j (f_j + (x xor x_j))$, $x < 2^(18)$], [split the bits: $O(sqrt(V))$ each],
  [edges usable if their mask $supset.eq m$], [one DSU per mask],
)

=== XOR and linear algebra
_Linear algebra facts_.
#cues(
  ["XOR of a subset": max, count], [`linalg/xor-basis`, $2^(n - r)$],
  [XOR of a walk, edges reusable], [cycle basis + root XORs],
  [span or max XOR of $a_l, dots, a_r$], [prefix basis: `linalg/xor-basis`],
  [toggles, Lights Out], [chase rows, then eliminate over $FF_2$],
  [rectangle sums or flips mod 2], [2D prefix or difference parities],
  [count solutions of parity rules], [$2^f$, $f$ free variables],
  [only "equal" and "different" rules], [parity DSU, 2-colouring],
  [rank of unions, disjoint bases], [random coefficients mod $p$],
  [traces of $A^k$, closed walk sums], [Newton: power sums give $det$],
  [$det(I + u v^T)$, rank-one updates], [matrix determinant lemma],
  [count spanning trees], [matrix-tree determinant],
  [does a perfect matching exist?], [random Tutte matrix],
)

=== Strings
_String facts_.
#cues(
  [prefix equals suffix (borders)], [`strings/prefix-function`],
  ["made of copies", periods, rotations], [period $n - pi_(n - 1)$, Fine-Wilf],
  [two palindromic prefixes], [a period; reflect $x -> 2 y - x$],
  [smallest rotation, suffix or cut], [Lyndon: `strings/lyndon`],
  [many patterns in one text], [`strings/aho-corasick`],
  [many words sharing prefixes], [a trie, as in `strings/aho-corasick`],
  [distinct substrings, $k$-th substring], [`strings/suffix-automaton`],
  [compare substrings, LCP], [`strings/string-hash`, suffix array],
  [longest palindrome at each centre], [`strings/manacher`],
  [distinct palindromic substrings], [`strings/eertree`],
  [squares, cubes, repetitions], [runs: fewer than $n$ of them],
  [delete one letter, compare the results], [runs of equal letters],
  [is it a subsequence, embeddings], [greedy leftmost embedding],
  [brackets: balanced, matching], [prefix balance $>= 0$, ends at 0],
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
  [around one point, or over directions], [angular sort, sweep, two pointers],
  [a beam bouncing off mirrors], [unfold the reflections],
  [regions or components of a drawing], [Euler: $V - E + F = 1 + C$],
  [polygon sides glued in pairs], [surface words: orientable, $chi$],
  [sum of areas over many polygons], [shoelace: count each edge],
  [groups split by lines], [lines through two input points],
  [lattice points in a polygon], [Pick; or column by column],
  [all pairs closer than $d$], [grid of cells, scan the neighbours],
  [integral of a density over a region], [Green's theorem, quadrature],
  [shortest path around obstacles], [tangents and arcs, visibility graph],
  [order points with only left turns], [always possible: a spiral],
  [join boundary pairs, no crossings], [no two pairs alternate],
  [triangle centres, circles, ellipses], [the formulas in _Geometry facts_],
)

=== DP design and speed-ups
_DP techniques_.
#cues(
  [$(h_i - h_j)^2$, $a_i b_j$, costs crossing once], [`ds/li-chao`, line container],
  ["exactly $k$ pieces", "at most $k$"], [`dp/wqs` if the cost is convex in $k$],
  [cost with the quadrangle inequality], [`dp/dc-opt`, `dp/knuth`],
  [costs $|x - a|$, "make it sorted"], [`dp/slope-trick`],
  [$sum_i b_i = m$, convex costs $f_i (b_i)$], [the $m$ cheapest unit steps],
  [one constraint couples convex terms], [binary search its multiplier (KKT)],
  [maximise $sum a \/ sum b$], [binary search $lambda$],
  [best with $c$ units from two parts], [concave: merge the slopes],
  [minimum over a sliding window], [monotone deque],
  [DP pays the maximum of each block], [stack of groups: Block maximum],
  [longest chain, increasing subsequence], [`dp/lis`; Dilworth],
  [choose a subset and its order], [exchange argument, sort, then DP],
  [subset sums, $n S <= 10^(10)$], [bitset: $O(n S \/ 64)$],
  [items with a small weight or value], [knapsack indexed by the small one],
  [grid with a short side ($<= 10$)], [broken profile DP],
  [lexicographically smallest optimum], [DP backwards, walk forwards],
  [best value and how many plans give it], [(best, count) pairs per state],
  [huge counts of a few item kinds], [greedy bulk, exact DP on a margin],
  [values blow up, only the sign matters], [saturate at a bound],
)

=== Ranges and updates
_Designing a segment tree node_, _Offline techniques_, _Sqrt
decomposition and small tricks_.
#cues(
  [point updates, range queries], [`ds/segtree`: design the node],
  [range answer by a left-to-right scan], [node = composed function],
  ["after each update, print the DP"], [(min, +) matrices, `ds/segtree`],
  ["first position where ..."], [descend the segment tree],
  [range add, sums of squares, history], [linear tags (matrices)],
  [range adds of linear functions], [second differences, summed twice],
  [prefix maxima of a range, updates], [Prefix maxima: $O(log^2 n)$],
  [components of the graph on a range], [node keeps boundary labels],
  [range chmin or chmax], [`ds2/segtree-beats`],
  [sweep $r$, answers for all $l$ at once], [all left ends: lazy tree over $l$],
  [good subarrays inside $[L, R]$], [rectangles in the $(l, r)$ plane],
  [distinct values in $[l, r]$], [offline by $r$, last occurrence],
  [$k$-th smallest in $[l, r]$], [`ds2/wavelet-matrix`, persistent],
  [questions about an old version], [`ds2/persistent-segtree`],
  [$a_j <= a_i$, $b_j <= b_i$, $c_j <= c_i$], [`offline/cdq`],
  [edges added and removed], [`offline/dynamic-connectivity`],
  [only deletions], [reverse time: they become inserts],
  ["the first time when ..." per query], [`offline/parallel-binary-search`],
  [static ranges, $O(1)$ add and remove], [`offline/mo`; rollback variant],
  [static range minimum], [`ds/sparse-table`],
  [for every $i$, all items except $i$], [all but one: divide and conquer],
  [a node's part depends on the query], [requests per segment-tree node],
  [cost from $|x - y|$ or $max - min$], [one 0/1 problem per threshold],
  [a cyclic sequence of transitions], [trace of the matrix product],
  [next greater, span of a maximum], [monotone stack],
  [pairs $(l, r)$ ruled by their maximum], [split at it, loop the shorter side],
  [sum of the max of all subarrays], [contribution of each maximum],
  ["at least one 0 in $[l, r]$"], [complement: runs of ones],
  [sum of the $t$ largest, any $t$], [tree over values, descend],
  [order statistics, insert at the $k$-th], [`ds2/ordered-set`, `ds2/treap`],
  [$max_j (a_j xor x)$], [`ds2/binary-trie`],
  [no $O(log n)$ node, many updates], [sqrt blocks, periodic rebuild],
  [queries with a step or modulus $d$], [tables for small $d$, walk large],
  [paint ranges, each cell once], [skip pointers (DSU)],
)

=== Trees
_Tree techniques_.
#cues(
  [per subtree, a multiset statistic], [small-to-large],
  [tree DP over sizes, $n <= 5000$], [tree knapsack: $O(n^2)$ in total],
  [DP indexed by the depth], [long paths: $O(n)$],
  ["for every vertex as the root"], [`tree/rerooting`],
  [all children but one, modulo $p$], [prefix and suffix products],
  [count paths or distances], [`tree/centroid`, FFT of depths],
  [nearest marked vertex, updates], [centroid ancestors],
  [queries on $k$ vertices, $sum k$ bounded], [`tree/virtual-tree`],
  ["only edges of weight $<= x$"], [Kruskal reconstruction tree],
  [path updates and path queries], [`tree/hld`],
  [walk up, pay at each new maximum], [next greater ancestor],
  [compare sequences read upwards], [hashes of $2^k$ steps up],
  [link and cut edges online], [`tree/link-cut`],
  [subtree sums with updates], [Euler tour + `ds/fenwick`],
  [cover or match greedily on a tree], [deepest vertex first],
  [pairs of leaves routed along paths], [Leaf pairings: edge loads],
  [sum over pairs of vertices or pieces], [per edge: $s (n - s)$ pairs],
  [diameter, eccentricity, centre], [two BFS; root at the centre],
  [LCA under any root, three paths meet], [the median of three vertices],
  [parents before children, orders], [Tree orders: $n! \/ product "size"$],
  [trees with given degrees], [Prüfer: degrees $>= 1$, sum $2n - 2$],
  [XOR with the children, $t$ times], [$binom(t, d)$ odd iff $d subset.eq t$],
)

=== Graphs and flows
_Graph theory facts_, _Flow and cut modelling_, _Matroids_.
#cues(
  [two options per object, pair rules], [`graph/two-sat`],
  ["at most one per group" in 2-SAT], [prefix variables],
  [edges from a vertex to a range], [segment-tree graph],
  [each item takes one of its two ends], [edge per item; $E <= V$ per part],
  [fewest deletions to separate], [min cut: `flow/dinic`],
  [taking $x$ forces $y$, max profit], [project selection],
  [assign with capacities], [bipartite flow, Hall],
  [lower and upper bounds per edge], [flow with lower bounds],
  ["is it always possible?", $n = 10^6$], [cut condition, then greedy],
  [most items, no two conflict, bipartite], [$|V| - nu$ (Kőnig)],
  [cheapest assignment or transport], [`flow/hungarian`, min-cost flow],
  [tokens moved on a grid, unit steps], [transport dual: potentials],
  [trade-off between two objectives], [hull with a min-cut oracle],
  [cover a DAG with disjoint paths], [$n - nu$],
  [use every edge exactly once], [`graph/euler-path`],
  [degrees fixed modulo 2 or $k$], [parity fixing on a spanning tree],
  [no cycle / independent / matchable], [matroid greedy],
  [rank, forest, matching of $[l, r]$], [latest basis: sweep $r$],
  [two such rules at once], [matroid intersection],
  ["differ by exactly 1 on each edge"], [bipartite: parity = colour],
  [closed walks of length mod $k$], [SCC period: labels mod $k$],
  [state = vertex plus extra info], [BFS or Dijkstra on pairs],
  [nearest of many sources], [multi-source Dijkstra, owners],
  [times set by earlier events], [Dijkstra in time order],
  [shortest cycle, chordless cycle], [BFS from every vertex],
  [MST of an implicit dense graph], [Borůvka with a data structure],
  [MST with small integer weights], [$sum_v (c_v - 1)$ over thresholds],
  [edges $(u + k, v + k)$ for all $k < L$], [DSU per power of two],
  [fewest swaps to sort], [$n - "cycles"$; adjacent: inversions],
  [bridges, cut vertices, cactus], [`graph/bcc`, block-cut tree],
  [is a free grid cell a cut vertex?], [ring of its 8 neighbours],
  [series-parallel, no $K_4$ minor], [contract degree 2, merge parallels],
  [reachability, cycles in a digraph], [`graph/scc`, condensation],
  [reachability from every vertex, DAG], [bitsets in topological order],
  [smallest order under precedences], [`graph/topo-sort` with a heap],
  [orient $K_n$, given out-degrees], [Landau, greedy build],
)

=== Games and interaction
_Game theory_, `misc/interactive`.
#cues(
  [sum of independent games], [Grundy values, XOR],
  [huge numbers, simple moves], [brute force table, period],
  [every play has the same length], [only its parity matters],
  [token game, visited vertices deleted], [wins iff start in all max matchings],
  [symmetric position], [mirror or pairing strategy],
  [different moves per player], [win/lose DP, no Grundy],
  [players collect points], [DP on my score minus yours],
  [chance moves, cycles of states], [value iteration],
  [hidden position, adaptive opponent], [play against the candidate set],
  [query limit near $log_2$(candidates)], [halving: `misc/interactive`],
  [replies linear in hidden numbers], [a linear system (least squares)],
  [sort with comparison queries], [`stable_sort`, the judge compares],
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
  [$k$ objects that must all differ], [colour coding],
  [random input, typical behaviour], [concentration: $O(sqrt(n))$ spread],
)

=== Constructions
_Constructive problems and invariants_.
#cues(
  ["print any, or $-1$"], [necessary condition, brute force],
  [reach B from A by moves], [invariant: parity, sum mod $k$, gcd],
  [tokens that cannot pass each other], [their order never changes],
  [the input is only $n$ (and $k$)], [brute force $n <= 8$, pattern],
  [a construction for every $n$], [induction, periodic blocks],
  [budget of $floor(k\/3)$ changes], [3 classes, take the smallest],
  [grid cells, no two touching], [4 classes $(i mod 2, j mod 2)$],
  [exactly $k$ of something], [gadgets that add up; subset sum],
  [regular: "each in exactly 3"], [Hall, perfect matchings],
  [a process that repairs violations], [decreasing potential],
  [fair split, envy-free up to one item], [envy-cycle elimination],
  [the answer always seems small], [prove $<= c$, then the rest],
  ["your program runs twice"], [encode a unique structure],
  [3 values: equal stays, else the third], [it is $-(x + y)$ mod 3: linear],
)

=== Greedy, search and ad hoc
_Technique finder_ (step 4), _Constructive problems and
invariants_ (exchange arguments).
#cues(
  [largest $k$ that is still possible], [binary search $k$, greedy check],
  [$k$-th smallest of $n^2$ pair values], [binary search the value, count],
  [best order of jobs or items], [exchange two neighbours, sort],
  [pair two lists, least max or sum], [sort both, pair in or against order],
  [jobs with deadlines], [earliest deadline; heap of the taken],
  [disjoint intervals, stabbing points], [sort by right end, greedy],
  [longest window with a monotone rule], [two pointers],
  [best $k$ items, a max-versus-sum rule], [a window of the sorted order],
  [reachable or feasible values], [an interval: keep both ends],
  [only a few candidates can be optimal], [prove which, try them all],
  [$n <= 10$ orders, $<= 20$ yes/no choices], [try all: $n!$ or $2^n$],
  [one parameter decides the rest], [enumerate it, solve the rest],
  [split into a left and a right part], [prefix and suffix aggregates],
  [long stretches of equal values], [run-length: work per run],
  [operations merge neighbours], [each result is a block: DP over cuts],
  [cyclic order, most $sum |x_i - x_(i + 1)|$], [top half minus bottom half],
  [coverage by many intervals], [sweep the sorted endpoints],
  [conditions that split by axis], [each axis alone, then combine],
  [cyclic array, cover a circle], [double it; total minus largest gap],
  [hard forwards, simple from the end], [run the process backwards],
  [lengths summing to $n$], [fewer than $sqrt(2n)$ distinct],
  [two costs trade off through a size], [equate them, as for sqrt blocks],
  [merge two for cost $a + b$, repeatedly], [Huffman: two sorted queues],
  [long statement, explicit rules], [simulate exactly: _Troubleshooting_],
)
