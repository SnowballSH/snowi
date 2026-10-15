#import "/typst/notes.typ": *
#let xor = sym.plus.o

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
- Permutation parity: a swap flips it, a 3-cycle keeps it,
  reversing a block of length $L$ flips it iff $floor(L\/2)$ is
  odd. 15 puzzle: the sign of the permutation of all cells (blank
  included) times $(-1)^d$, with $d$ the blank's taxicab distance
  from its home cell, never changes; on an $m times n$ board with
  $m, n >= 2$ it is the only obstruction.
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
  bad edges, the sorted state in lexicographic order.]
- Local repair: start anywhere; while a constraint is violated,
  fix it so that $Phi$ drops. Example: put each vertex on one of
  two sides so that $v$ has at most $deg(v)\/2$ neighbours on its
  own side. Flip any vertex with more; the number of edges inside
  a side drops by at least 1, so this ends within $m$ flips.
- Example: a move subtracts 1 from $a_i$ and $a_(i + 1)$; can all
  values reach 0? The alternating sum $sum_i (-1)^i a_i$ is
  invariant, so it must be 0, but that is not enough
  ($[1, 0, 0, 1]$). Only move 1 touches $a_1$, so it is used
  $x_1 = a_1$ times, then $x_i = a_i - x_(i - 1)$. The answer is
  yes iff every $x_i >= 0$ and $x_(n - 1) = a_n$. A forced greedy
  is the construction and the proof at once.

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
  odd cycle is the only obstruction to a 2-colouring.
- Euler tours need even degrees (in = out when directed) and one
  connected component with edges; pair the odd vertices first
  (_Graph theory facts_).

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
- Almost-solution, then repair: meet every local constraint (a
  perfect matching, degrees mod $k$), then fix the global one
  (merge the cycles into one, fix the total parity) by small
  surgeries.
- Periodic blocks: tile with a small gadget, then fix the borders
  and small $n$ separately; most wrong answers come from $n <= 3$
  or odd sizes.
- Exactly $k$: use gadgets whose contributions add up and take the
  largest that fits ($binom(m, 2)$, powers of two), or a bitset
  subset sum over the gadget sizes.
- Read the budget: $floor(k\/3)$ means 3 classes, $2 log_2 n$
  means two queries per halving (or doubling) step, $3 n$ means
  $O(1)$ moves per element.
- "Minimise $d$" with tiny answers: prove $d <= c$ by a
  construction, then characterise $d = 1$, $d = 2$, ... one by one.
- Two runs (encode, then decode after the judge shuffles or
  relabels): make the message the unique structure of its kind,
  or add parity bits and XOR checksums.

=== Standard constructions
#fact[Gray code][$g_i = i xor floor(i\/2)$ (`i ^ i >> 1`):
  consecutive codes differ in one bit, cyclically too. Inverse:
  $i$ is the XOR of $floor(g\/2^j)$ over all $j >= 0$. Recursively: $0$
  followed by $G_(n - 1)$, then $1$ followed by $G_(n - 1)$
  reversed.]
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
#fact[Tournament scores][Landau: $s_1 <= dots <= s_n$ are the
  out-degrees of a tournament iff
  $sum_(i <= k) s_i >= binom(k, 2)$ for all $k$, with equality at
  $k = n$. Greedy: the vertex with the smallest score $s$ beats the
  $s$ next smallest and loses to the rest, whose scores drop by 1;
  repeat. It fails only if no tournament exists.
  $O(n^2 log n)$.]
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
#fact[Permutations][Exactly $k$ inversions, which needs
  $0 <= k <= binom(n, 2)$: position $i = 0, dots, n - 1$ takes the
  unused value with exactly $c_i$ smaller unused values, where
  $c_i = min(k', n - 1 - i)$ and $k'$ is what is left of $k$.
  Longest increasing $a$ and longest decreasing $b$ exist
  together iff $a + b - 1 <= n <= a b$: $a$ blocks of consecutive
  values, each block decreasing, one of size $b$ and none larger,
  blocks in increasing order. Always $"LIS" dot "LDS" >= n$ and
  $"LIS" + "LDS" <= n + 1$.]
- Degree sequences: Havel-Hakimi, Erdős-Gallai, Gale-Ryser
  (_Graph theory facts_). A $k$-regular bipartite multigraph
  splits into $k$ perfect matchings: $k$ disjoint rounds, or an
  edge colouring with $Delta$ colours (_Flow and cut modelling_).

=== Exchange arguments
Take an optimal solution and the first place where it differs from
the greedy one; swap it towards the greedy choice without making it
worse. For an order, compare neighbours $a, b$: $a$ goes first iff
$"cost"(a, b) <= "cost"(b, a)$. Check that this comparison is
transitive, and cross-multiply instead of dividing.
- Minimise $sum_i w_i C_i$ (weighted completion times): sort by
  $t_i \/ w_i$ (Smith's rule).
- Minimise the maximum lateness: earliest deadline first.
- Most jobs on time: by deadline, add each job; if the total
  time exceeds its deadline, drop the longest job so far (a
  max-heap; Moore-Hodgson).
- Two machines in sequence (times $a_i$, then $b_i$), least
  makespan: jobs with $a_i < b_i$ by increasing $a_i$, then the
  others by decreasing $b_i$ (Johnson).
- Smallest concatenation of strings: sort by $a + b < b + a$ (this
  order is transitive).
- Intervals: the largest number of pairwise disjoint ones equals
  the fewest points that stab them all; both come from a scan by
  increasing right end.
- Canonical forms: show that some optimum is sorted, contiguous
  or has at most one fractional item, and search only those.
- Unit jobs with deadlines and profits: _Matroids_.
