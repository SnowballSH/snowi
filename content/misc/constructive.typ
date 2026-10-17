#import "/typst/notes.typ": *
#let xor = sym.plus.o
#full[

Workflow: find a necessary condition (an invariant, a colouring, a
count), brute force every small input to see whether it is also
sufficient and what the answers look like, then build the answer
and compare it with the brute force through a checker for all
small inputs, $n = 1, 2, 3$ included.

=== Invariants and monovariants
#fact[Invariant][A quantity $I$ that no move changes: $I(A) != I(B)$
  proves that $B$ is unreachable from $A$. It is only necessary;
  compare it with a BFS over all states of small instances to see
  whether it is the whole answer.]
- Parity of a sum or of the number of odd values; the sum mod $k$
  when a move adds the same amount to $k$ cells.
- A weighted sum $sum_i w_i a_i$ with weights chosen so that every
  move sums to 0: alternating signs for "subtract 1 from two
  neighbours".
- Combine constraints: add (or XOR) several of them so that every
  unknown but one cancels. If the rows, columns and diagonals of a
  $3 times 3$ square XOR to 0, then rows 1 and 3, column 2 and both
  diagonals hold the centre three times and every other cell an
  even number of times, so the centre is 0.
- Three values, where equal inputs give that value and different
  ones give the third: on $0, 1, 2$ this is $-(x + y) mod 3$, a
  linear map (on $1, 2, 3$ the third value is $x xor y$). Chains
  of it compute signed sums; applied to all neighbours of
  $a_0, dots, a_(n - 1)$ until one value is left, it gives
  $(-1)^(n - 1) sum_i binom(n - 1, i) a_i mod 3$ (Lucas).
]
- Permutation parity: a swap flips it, a 3-cycle keeps it,
  reversing a block of length $L$ flips it iff $floor(L\/2)$ is
  odd. 15 puzzle: the sign of the permutation of all cells (blank
  included) times $(-1)^d$, with $d$ the blank's taxicab distance
  from its home cell, never changes; on an $m times n$ board with
  $m, n >= 2$ it is the only obstruction.
#full[
- Swapping the tokens on the ends of any edge of a graph reaches
  exactly the arrangements that keep every token inside its
  connected component.
- The gcd of the values under $x -> x - y$; under mirror jumps
  $x -> 2y - x$ the gcd $g$ of all differences and each value mod
  $2g$.
- The sum of positions when a move shifts one token by $+d$ and
  another by $-d$. Relay jumps (the token to move jumps from $x$
  over the token at $y$ to $2 y - x$, and the token at $y$ moves
  next): the sum of positions minus twice the position of the
  token to move. The sum itself changes.
- The order of tokens that cannot pass each other (colliding ants
  pass through each other with their labels swapped).
- Degrees mod 2 (there is an even number of odd vertices), the
  number of components, the number of cells of each colour.
#fact[Monovariant][An integer $Phi >= 0$ that every move strictly
  decreases allows at most $Phi("start")$ moves. Typical $Phi$:
  inversions, $sum_i a_i^2$ under averaging moves, the number of
  bad edges, the sorted state in lexicographic order. With costs:
  if a move of cost $c$ lowers $Phi$ by at most $c$, every
  solution costs at least $Phi("start") - Phi("end")$, and one
  whose moves lower $Phi$ by exactly their cost is optimal.]
- Local repair: start anywhere; while a constraint is violated,
  fix it so that $Phi$ drops. Example: put each vertex in one of
  $k$ classes so that $v$ has at most $floor(deg(v)\/k)$
  neighbours in its own class. Move a vertex with more to the
  class least frequent among its neighbours; the number of edges
  inside classes drops by at least 1, so this ends within $m$
  moves ($k = 2$: a cut with at least $m\/2$ edges).
- Example: a move subtracts 1 from $a_i$ and $a_(i + 1)$; can all
  values reach 0? The alternating sum $sum_i (-1)^i a_i$ is
  invariant, so it must be 0, but that is not enough
  ($[1, 0, 0, 1]$). Only move 1 touches $a_1$, so it is used
  $x_1 = a_1$ times, then $x_i = a_i - x_(i - 1)$. The answer is
  yes iff every $x_i >= 0$ and $x_(n - 1) = a_n$. A forced greedy
  is the construction and the proof at once.
#fact[Discrete continuity][If every move changes an integer $f$ by
  at most 1, every value between $f(P)$ and $f(Q)$ occurs on the
  way from $P$ to $Q$. If such moves connect all candidates, the
  achievable values form an interval: find its two ends. Examples:
  replacing one element of a subsequence changes its LIS by at
  most 1, so the LIS values of the length-$m$ subsequences of a
  sequence form an interval; monotone lattice paths with the same
  ends are connected by flips of one corner cell.]

=== Parity and colouring
- Chessboard colouring: a domino covers one black and one white
  cell. Two opposite corners have the same colour iff $m + n$ is
  even, and then the board without them has no tiling. A
  knight changes colour on every move, so a closed knight's tour
  has even length.
- $1 times k$ tiles: colour cell $(i, j)$ with $(i + j) mod k$.
  Every tile covers each colour once, so a tiling needs all $k$
  colour counts equal.
- The 4 classes $(i mod 2, j mod 2)$ are each king-independent (no
  two cells of a class touch, even diagonally). Pigeonhole: of $N$
  objects some class gets at least $ceil(N\/4)$ and some at most
  $floor(N\/4)$. A budget of $floor(k\/3)$ changes means 3 classes:
  use the smallest.
- Bipartite graphs: every closed walk has even length, and
  $"dist"(u, v)$ is even iff $u$ and $v$ have the same colour. An
  odd cycle is the only obstruction to a 2-colouring. So when all
  vertices of a tree are paired up, the total distance has the
  parity of the number of black vertices.
- The $n times m$ grid ($n, m >= 2$) has a Hamiltonian cycle iff
  $n m$ is even (a cycle alternates colours). For even $n$: along
  row 1, snake over columns $2, dots, m$ of the other rows, then
  back up column 1. If both are odd, all cells but one corner
  form a cycle. A snake that follows a Hamiltonian cycle never
  runs into itself.
- Euler tours need even degrees (in = out when directed) and one
  connected component with edges; pair the odd vertices first
  (_Graph theory facts_).

=== Bounds
- Double counting: if each of $N$ windows (or constraints) must
  contain a bad element and every element lies in at most $c$ of
  them, at least $ceil(N\/c)$ elements are bad.
- Triangle-free sorted lengths grow like Fibonacci numbers:
  _Numbers, estimates and numerics_.
- Pigeonhole halving: if a range of values holds fewer elements
  than it has values, it misses one, and so does one of its
  halves. One counting pass per halving finds a missing value
  in $O(log n)$ passes and $O(1)$ memory (the counter can stop at
  the half's size).

=== Building the answer
- Brute force first: all objects for $n <= 8$, or a BFS over the
  states for "is it reachable?". Print the lexicographically
  smallest answers, as pictures for grids; many solutions means a
  simple rule exists.
- Induction: extend a solution for $n - 1$, $n - 2$ or $n\/2$.
  Strengthen the hypothesis until the step works ("the tour of a
  subtree starts at its root and ends next to it"); peel a layer
  off a grid; solve the base cases by brute force.
- Extremal choices: handle the largest, the leftmost or a leaf
  first. A forced choice is never wrong.
- Think backwards: what must hold just before the last move? Run
  the process in reverse: release times become deadlines, and the
  last slot can take anyone.
- Exact test, then greedy: with a test "can the rest still be
  completed?", fill positions left to right with the smallest
  choice that passes it. This never gets stuck and gives the
  lexicographically smallest answer. Example: a multiset with
  fixed outside neighbours $L$ and $R$ (when present) can be
  ordered with no two equal neighbours iff
  $c_x - ("len" - c_x) + [L = x] + [R = x] <= 1$ for every value
  $x$, with $c_x$ its count and len the multiset's size.
- Almost-solution, then repair: meet every local constraint (a
  perfect matching, degrees mod $k$), then fix the global one
  (merge the cycles into one, fix the total parity) by small
  surgeries. Two cycles that use opposite sides of a 4-cycle (a
  unit square in a grid) become one when those two edges are
  replaced by the other two sides.
- Periodic blocks: tile with a small gadget, then fix the borders
  and small $n$ separately; most wrong answers come from $n <= 3$
  or odd sizes.
- Exactly $k$: use gadgets whose contributions add up and take the
  largest that fits ($binom(m, 2)$, powers of two), a bitset
  subset sum over the gadget sizes, or independent gadgets whose
  counts multiply (7 gadgets with 4 states give $4^7$).
- Read the budget: $floor(k\/3)$ means 3 classes, $2 log_2 n$
  means two queries per halving (or doubling) step, $3 n$ means
  $O(1)$ moves per element.
- "Minimise $d$" with tiny answers: prove $d <= c$ by a
  construction, then characterise $d = 1$, $d = 2$, ... one by one.
- Two runs (encode, then decode after the judge shuffles or
  relabels): make the message the unique structure of its kind,
  or add parity bits and XOR checksums. To add one element to a
  set $S$ of $k < n\/2$ values from $1, dots, n$ so that it can be
  recovered from the result: read $1, dots, n$ as brackets (`(`
  for members, `)` for the rest), match them greedily left to
  right and add the last unmatched `)`. The decoder runs the same
  matching on the new set and outputs the bottom of the stack.
- Hat guessing: a guess that is a linear function of the visible
  colours is right iff a linear form of the whole configuration
  vanishes. Over $FF_q$, give $q + 1$ people pairwise
  non-proportional forms from one 2-dimensional space of forms,
  each nonzero on its owner's colour and zero on colours the owner
  cannot see. One of them vanishes on every configuration
  (Subspaces in _Linear algebra facts_), so someone guesses right.
  If everyone sees all other hats and may
  first announce a bit: use $2^(n - 1)$ colours as $n$-bit codes
  of even parity; person $i$ announces bit $i$ of the XOR of the
  codes it sees. Then everyone knows the other bits of its own
  code (bit $j$ of the total XOR is announcement $j$ XOR bit $j$
  of person $j$'s code) and gets bit $i$ from the parity.

]
=== Standard constructions
#full[
- Gray code: _Bit tricks_.
]
#fact[de Bruijn][$B(k, n)$ is a cyclic word of length $k^n$ that
  contains every length-$n$ word over $k$ letters exactly once:
  the Lyndon words whose length divides $n$, concatenated in
  lexicographic order (or an Euler circuit on the words of length
  $n - 1$). Append its first $n - 1$ letters for a linear word.
  $O(k^n)$.]
#lang-code(
  cpp: "vi deBruijn(int k, int n) { // letters 0..k-1, cyclic
  vi seq, w{-1};
  while (sz(w)) {
    w.back()++;
    int m = sz(w);
    if (n % m == 0) seq.insert(seq.end(), all(w));
    while (sz(w) < n) w.push_back(w[sz(w) - m]);
    while (sz(w) && w.back() == k - 1) w.pop_back();
  }
  return seq;
}",
  py: "def de_bruijn(k, n):  # letters 0..k-1, cyclic
  seq, w = [], [-1]
  while w:
    w[-1] += 1
    m = len(w)
    if n % m == 0:
      seq += w
    while len(w) < n:
      w.append(w[len(w) - m])
    while w and w[-1] == k - 1:
      w.pop()
  return seq",
)
#fact[Latin squares][$L_(i j) = (i + j) mod n$. For odd $n$ the
  squares $(i + j) mod n$ and $(i + 2 j) mod n$ are orthogonal
  (all $n^2$ pairs of entries differ); for a prime $p$ the squares
  $(a i + j) mod p$, $a = 1, dots, p - 1$, are mutually orthogonal.
  An orthogonal pair exists for every $n$ except 2 and 6.]
#fact[Magic square][Odd $n$, $0 <= i, j < n$, residues in
  $[0, n)$: $M_(i j) = n ((i + j + (n + 1)\/2) mod n) +
  ((i - j + (n - 1)\/2) mod n) + 1$ uses $1, dots, n^2$, and every
  row, column and both diagonals sum to $n (n^2 + 1)\/2$.]
#fact[Round robin][Even $n$: in round $r = 0, dots, n - 2$ pair
  player $n - 1$ with $r$, and $(r + i) mod (n - 1)$ with
  $(r - i) mod (n - 1)$ for $1 <= i < n\/2$. Every pair meets
  exactly once (the edges of $K_n$ split into perfect matchings).
  Odd $n$: add a dummy player, whose partner rests.]
#fact[Hamiltonian decomposition][$K_(2m)$ splits into the $m$
  zigzag paths $i, i + 1, i - 1, i + 2, i - 2, dots, i + m$
  (mod $2m$), $i = 0, dots, m - 1$. For $K_(2m + 1)$ join one new
  vertex to both ends of each path: $m$ Hamiltonian cycles.]
#full[
#fact[Tournament scores][Which out-degree lists occur: Landau
  (_Graph theory facts_). To build one, the vertex with the
  smallest score $s$ beats the $s$ next smallest and loses to the
  rest, whose scores drop by 1; repeat. It fails only if no
  tournament exists. $O(n^2 log n)$.]
#lang-code(
  cpp: "// beats[u][v]: u beats v. Empty if s is not a score list.
vector<vector<bool>> tournament(vi s) {
  int n = sz(s);
  vector<vector<bool>> beats(n, vector<bool>(n));
  vi a(n);
  iota(all(a), 0);
  while (sz(a)) {
    sort(all(a), [&](int x, int y) { return s[x] < s[y]; });
    int u = a[0], k = s[u];
    if (k < 0 || k >= sz(a)) return {};
    rep(i, 1, sz(a)) {
      if (i <= k) beats[u][a[i]] = 1;
      else beats[a[i]][u] = 1, s[a[i]]--;
    }
    a.erase(a.begin());
  }
  return beats;
}",
  py: "def tournament(s):  # beats[u][v]; None if s is no score list
  n, s = len(s), s[:]
  beats = [[0] * n for _ in range(n)]
  alive = list(range(n))
  while alive:
    alive.sort(key=s.__getitem__)
    u = alive.pop(0)
    k = s[u]
    if not 0 <= k <= len(alive):
      return None
    for i, v in enumerate(alive):
      if i < k:
        beats[u][v] = 1
      else:
        beats[v][u] = 1
        s[v] -= 1
  return beats",
)
]
#fact[Permutations][Exactly $k$ inversions, which needs
  $0 <= k <= binom(n, 2)$: position $i = 0, dots, n - 1$ takes the
  unused value with exactly $c_i$ smaller unused values, where
  $c_i = min(k', n - 1 - i)$ and $k'$ is what is left of $k$.
  Longest increasing $a$ and longest decreasing $b$ exist
  together iff $a + b - 1 <= n <= a b$: $a$ blocks of consecutive
  values, each block decreasing, one of size $b$ and none larger,
  blocks in increasing order. Always $"LIS" dot "LDS" >= n$ and
  $"LIS" + "LDS" <= n + 1$.]
#full[
- Appending an element raises at most one of LIS and LDS. Every
  prefix can be given its own LIS and LDS: grow a Young diagram
  one corner cell at a time and give cell (column $c$, row $r$)
  the rank of the key $(c, -r)$ among all cells; each prefix then
  has LIS = its number of columns and LDS = its number of rows.
- Swaps that only create inversions (swap $p_i < p_j$ with
  $i < j$) turn $p$ into $q$ iff for all $i$ and $j$ at least as
  many of $p_1, dots, p_i$ as of $q_1, dots, q_i$ are at most $j$.
  Build it by placing the values from the smallest: move $x$ right
  by swapping it with the smallest larger value between its
  position and its target.
- Positions and values: swapping positions $i$ and $j$ of $p$ swaps
  the values $i$ and $j$ in $p^(-1)$; pick the view in which the
  moves are local. Adjacent swaps at $i, i + 1, dots, j$ in this
  order rotate $a_i, dots, a_(j + 1)$ left by one.
#fact[Sorting networks][Odd-even transposition sort compares and
  swaps the pairs $(0, 1), (2, 3), dots$ and $(1, 2), (3, 4), dots$
  in alternate rounds; $n$ rounds sort any $n$ values, whichever
  round comes first. 0-1 principle: a fixed sequence of
  compare-swaps sorts every input iff it sorts every 0/1 input, so
  it is enough to reason about "large" (at least a threshold)
  versus "small" values.]
#fact[Prefix codes][A prefix-free binary code with word lengths
  $l_1 <= dots <= l_k$ exists iff $sum_i 2^(-l_i) <= 1$ (Kraft).
  Canonical code: each word is the previous one plus 1, shifted
  left to its own length. A protocol may finish with outcomes of
  different lengths as long as these weights sum to at most 1.]
#fact[Envy-cycle elimination][Items arrive one by one; valuations
  are monotone. Give the next item to someone nobody envies. If
  everyone is envied, follow "is envied by" edges to a cycle, let
  each person on it take the bundle it envies, and retry. The
  result is EF1: whoever envies a bundle stops once some single
  item is removed from it. Only the bundle that grew needs new
  values: $n$ valuation queries per item.]
- Degree sequences: Havel-Hakimi, Erdős-Gallai, Gale-Ryser
  (_Graph theory facts_). A tree with degrees $d_i >= 1$ exists iff
  $sum_i d_i = 2 n - 2$ ($n >= 2$): join a vertex of remaining
  degree 1 to one of remaining degree at least 2, repeat, then join
  the last two. A $k$-regular bipartite multigraph as $k$ rounds of
  perfect matchings, or a bipartite edge colouring with $Delta$
  colours: _Graph theory facts_.

=== Exchange arguments
Take an optimal solution and the first place where it differs from
the greedy one; swap it towards the greedy choice without making it
worse. For an order, compare neighbours $a, b$: $a$ goes first iff
$"cost"(a, b) <= "cost"(b, a)$. Check that this comparison is
transitive, and cross-multiply instead of dividing.
- Smallest concatenation of strings: sort by $a + b < b + a$ (this
  order is transitive).
- Pairing two lists: $sum_i a_i b_(pi(i))$ is largest with both
  sorted the same way and smallest with opposite orders. Opposite
  orders also maximise $min_i (a_i + b_(pi(i)))$ and minimise
  $max_i (a_i + b_(pi(i)))$.
- Uncrossing: if two chosen tuples can be replaced by their
  coordinatewise maximum and minimum without becoming invalid or
  worse (a product of positive coordinates, say), some optimum is
  a chain: match the largest with the largest. Example: disjoint
  index quadruples $i < j < k < l$ spelling K, U, P, C with the
  largest $sum i j k l$: scan from the right, a FIFO queue per
  letter.
- Cyclic orders: $sum_i max(0, x_(i + 1) - x_i)$ around a cycle is
  half of $sum_i |x_(i + 1) - x_i|$, and at most the sum of the
  $floor(m\/2)$ largest of the $m$ values minus the sum of the
  $floor(m\/2)$ smallest; alternate large and small values to reach
  it. On a closed tour, step costs $max(Delta a, Delta b) =
  Delta b + max(Delta (a - b), 0)$ reduce to this with
  $x = a - b$.
- Canonical forms: show that some optimum is sorted, contiguous
  or has at most one fractional item, and search only those.

]
=== Scheduling
#full[
- Minimise $sum_i w_i C_i$ (weighted completion times): sort by
  $t_i \/ w_i$ (Smith's rule).
- Minimise the maximum lateness: earliest deadline first. With
  release times: at every moment run the released job with the
  earliest deadline; for unit jobs (or with preemption) this meets
  all deadlines whenever any schedule does.
]
- Most jobs on time: by deadline, add each job; if the total
  time exceeds its deadline, drop the longest job so far (a
  max-heap; Moore-Hodgson).
- From the end: the total time $T$ is known, so choose the last
  job. To minimise $max_j f_j (C_j)$ with non-decreasing $f_j$,
  put last the job with the smallest $f_j (T)$ and repeat on the
  rest (Lawler). Least $sum_j C_j$ subject to deadlines
  $C_j <= d_j$: put last the longest job with $d_j >= T$.
- Two machines in sequence (times $a_i$, then $b_i$), least
  makespan: jobs with $a_i < b_i$ by increasing $a_i$, then the
  others by decreasing $b_i$ (Johnson). If job $j$ takes $W_j$ on
  each of $m$ machines in sequence (unlimited buffers), every
  order has makespan $sum_j W_j + (m - 1) max_j W_j$.
#full[
- Cache of size $k$ with all requests known: on a miss evict the
  item whose next request is farthest away (Belady); this
  minimises the misses.
- Unit jobs with deadlines and profits: _Matroids_.

=== Intervals and circles
- The largest number of pairwise disjoint intervals equals the
  fewest points that stab them all; both come from a scan by
  increasing right end.
- The fewest intervals whose indicator vectors sum to $v >= 0$
  ("add 1 on a range" operations): $sum_i max(0, v_i - v_(i - 1))$
  with $v_0 = 0$.
- Split intervals into $k$ groups to maximise the total length of
  the groups' unions: the optimum is the integral of
  $min(k, "depth"(t))$. Give the intervals, by left end, each to
  the group whose union ends first.
- Arcs on a circle, fewest stabbing points: put the first point at
  the end of some arc, then stab the rest greedily as on a line.
  The best over all arcs as the start is optimal; with next
  pointers and binary lifting each start costs $O(log n)$.
- A circle cut into at most $k$ arcs of sum at most $M$: run the
  line greedy from position 0. A greedy segment other than the
  last, with the element that did not fit, sums to more than $M$,
  so every feasible cutting starts an arc at one of its elements
  after the first. Try those of the shortest such segment as
  starts: $O(n)$ steps of the `nxt` pointers in total.
]
