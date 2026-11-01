#import "/typst/notes.typ": *

$S(n, k)$ and $c(n, k)$ are Stirling numbers of the second and
(unsigned) first kind, $x^(underline(k)) = x (x - 1) dots
(x - k + 1)$ and $x^(overline(k)) = x (x + 1) dots (x + k - 1)$.

=== Binomials
- Upper indices:
  $sum_(k = 0)^N binom(k, a) binom(N - k, b) = binom(N + 1, a + b + 1)$
  ($k$ is the $(a + 1)$-th chosen element).
- Bounded stars and bars: $x_1 + dots + x_k = n$ with
  $0 <= x_i <= c$ has
  $sum_j (-1)^j binom(k, j) binom(n - j (c + 1) + k - 1, k - 1)$
  solutions; for a sum $<= n$ use $binom(n - j (c + 1) + k, k)$.
- Over the compositions of $n$ into $p$ positive parts, the
  product of $r$ fixed parts sums to $binom(n + r - 1, p + r - 1)$.
#full[
- $binom(n, m) binom(m, k) = binom(n, k) binom(n - k, m - k)$ (pick
  the inner set first).
- Vandermonde: $sum_i binom(m, i) binom(n, k - i) = binom(m + n, k)$;
  squares: $sum_k binom(n, k)^2 = binom(2n, n)$.
- Hockey stick: $sum_(i = r)^n binom(i, r) = binom(n + 1, r + 1)$,
  $sum_(i = 0)^k binom(n + i, i) = binom(n + k + 1, k)$.
]

=== Paths and Catalan numbers
#fact[Reflection][$plus.minus 1$ walks of $n$ steps from $a$ to $b$:
  $W(b - a)$ with $W(d) = binom(n, (n + d)\/2)$ (0 if parity
  differs). For $a, b >= 0$, those that never go below 0:
  $W(b - a) - W(a + b + 2)$ (reflect the start across $-1$). Inside
  $[0, w]$: $sum_(j in ZZ) (W(b - a + 2j(w + 2)) - W(-a - b - 2 +
  2j(w + 2)))$, $O(n\/w)$ nonzero terms. Grid paths to $(n, m)$,
  $m <= n$, never more ups than rights:
  $binom(n + m, m) - binom(n + m, m - 1)$.]
#full[#fact[Sticky walk][$x_0 = 0$, $x_t = max(0, x_(t - 1) + s_t)$
  for steps $s_t = plus.minus 1$ (Kadane: $max_t x_t$ is the largest
  subarray sum, empty allowed). The free walk $y$ from 0 that steps
  $s_t$ from $y >= 0$ and $-s_t$ from $y < 0$ is a bijection with
  $x = y$ or $x = -1 - y$, so $max_t x_t <= k$ iff $y$ stays in
  $[-k - 1, k]$: a strip, as above.]
#fact[Ballot][$A$ gets $a$ votes, $B$ gets $b < a$: $A$ stays
  strictly ahead during a random count with probability
  $(a - b)\/(a + b)$.]]
#fact[Catalan][1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796;
  #full[$C_n = binom(2n, n)\/(n + 1)$,
  $C_(n + 1) = sum_(i = 0)^n C_i C_(n - i)$, ]$C_n = binom(2n, n) -
  binom(2n, n + 1)$, $C_(n + 1) = (4n + 2)\/(n + 2) dot C_n$.#full[ Counts balanced
  bracket strings with $n$ pairs, Dyck paths, binary trees with $n$
  nodes, plane trees with $n + 1$ nodes, full binary trees with
  $n + 1$ leaves, triangulations of an $(n + 2)$-gon, non-crossing
  matchings of $2n$ points, non-crossing partitions of $n$,
  231-avoiding (stack-sortable) permutations.]]

=== Stirling, Bell, partitions
#fact[First kind][$c(n, k)$: permutations of $n$ with $k$ cycles,
  or with $k$ prefix maxima.
  $c(n + 1, k) = n c(n, k) + c(n, k - 1)$;
  $x^(overline(n)) = sum_k c(n, k) x^k$; $c(n, 1) = (n - 1)!$.]
#fact[Second kind][$S(n, k)$: partitions of $n$ labelled items into
  $k$ non-empty unlabelled blocks. $S(n + 1, k) = k S(n, k) +
  S(n, k - 1)$; $x^n = sum_k S(n, k) x^(underline(k))$;
  $S(n, k) = 1/(k!) sum_i (-1)^(k - i) binom(k, i) i^n$ (one row by
  convolution).]
#fact[Bell][$B_n = sum_k S(n, k)$: all set partitions.
  $B_(n + 1) = sum_k binom(n, k) B_k$. 1, 1, 2, 5, 15, 52, 203, 877,
  4140, 21147, 115975.]
#fact[Partitions][$p(n)$: 1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42.
  Pentagonal recurrence in
  $O(n sqrt(n))$: $p(n) = sum_(k >= 1) (-1)^(k + 1) (p(n - g_k) +
  p(n - g_(-k)))$, $g_k = k (3k - 1)\/2$, $p(m < 0) = 0$. Exactly
  $k$ parts:
  $p_k (n) = p_(k - 1) (n - 1) + p_k (n - k)$. Distinct parts and
  odd parts are equinumerous.]

=== Twelvefold way
#block(breakable: false)[
$n$ balls into $k$ boxes; L = labelled, U = identical. $p_k (n)$ =
partitions of $n$ into exactly $k$ parts. With $k = 0$ boxes every
count is $[n = 0]$; the U, L binomials assume $k >= 1$.
#table(
  columns: 4,
  [balls, boxes], [any], [at most 1 each], [at least 1 each],
  [L, L], [$k^n$], [$k^(underline(n))$], [$k! S(n, k)$],
  [U, L], [$binom(n + k - 1, n)$], [$binom(k, n)$],
  [$binom(n - 1, k - 1)$],
  [L, U], [$sum_(i <= k) S(n, i)$], [$[n <= k]$], [$S(n, k)$],
  [U, U], [$p_k (n + k)$], [$[n <= k]$], [$p_k (n)$],
)]

=== Permutations and tableaux
#fact[Derangements][$D_n = (n - 1)(D_(n - 1) + D_(n - 2))$
  $= n D_(n - 1) + (-1)^n approx n!\/e$, starting 1, 0, 1, 2, 9, 44,
  265, 1854, 14833. Exactly $k$ fixed points: $binom(n, k) D_(n - k)$.]
#fact[Eulerian][$A(n, k)$: permutations with $k$ descents.
  $A(n, k) = (k + 1) A(n - 1, k) + (n - k) A(n - 1, k - 1)$; rows 1;
  1 1; 1 4 1; 1 11 11 1; 1 26 66 26 1. Worpitzky:
  $x^n = sum_k A(n, k) binom(x + k, n)$.]
#fact[Position bounds][Permutations with $p_i <= b_i$ for all $i$:
  cap the $b_i$ at $n$ and sort them ascending, then
  $product_(i = 1)^n max(0, b_((i)) - i + 1)$ (fill the tightest
  position first). For $p_i >= a_i$ apply it to
  $n + 1 - p_i <= n + 1 - a_i$.]
#fact[Prefix maxima][Insert the values $n, n - 1, dots, 1$: $x$ is a
  prefix maximum iff it goes in front (1 of $n - x + 1$ places). So
  a given set $S in.rev n$ of prefix-maximum values has
  $product_(x in.not S) (n - x)$ permutations.]
#fact[Hook length][Standard Young tableaux of shape $lambda$:
  $f^lambda = n! \/ product "hook"(c)$, the hook of cell $c$ being $c$,
  the cells right of it in its row and below it in its column.
  RSK: permutations biject
  with pairs of tableaux of one shape ($sum (f^lambda)^2 = n!$), with
  LIS $= lambda_1$ and LDS = number of rows. Greene: $k$ disjoint
  increasing subsequences cover at most $lambda_1 + dots + lambda_k$
  elements, and some reach it; for decreasing ones, the $k$ longest
  columns. Involutions biject with
  single tableaux, fixed points = number of odd-length columns.]
#fact[Tree orders][Orders of the vertices of a rooted forest with
  every vertex before its descendants: $n! \/ product_v "size"(v)$.]
#fact[Square roots][$q = p^2$ for some permutation $p$ iff for
  every even $L$ the number of $L$-cycles of $q$ is even (squaring
  keeps an odd cycle and splits an even one in two halves). So the
  reversal $i -> n + 1 - i$ is a square iff $n mod 4 in {0, 1}$.]

=== Inclusion-exclusion and inversion
#fact[Binomial][$g_k = sum_(j >= k) binom(j, k) f_j <=> f_k =
  sum_(j >= k) (-1)^(j - k) binom(j, k) g_j$. "Choose $k$, force
  them" gives $g_k$, not "at least $k$". For all $k$ at once:
  $k! f_k = sum_j j! g_j dot (-1)^(j - k)\/(j - k)!$, one
  convolution after reversing.]
#full[#fact[Subsets][$g(S) = sum_(T subset.eq S) f(T) <=> f(S) =
  sum_(T subset.eq S) (-1)^(|S| - |T|) g(T)$; $O(2^n n)$ with
  `poly/subset-transforms`.]
#fact[Multiples][With $g(d)$ counting objects whose values are all
  divisible by $d$, "gcd exactly $d$" is $f(d) = g(d) -
  sum_(m = 2d, 3d, dots) f(m)$ for $d$ descending, $O(V log V)$; or
  Möbius over multiples (Number theory facts).]]
#fact[Stirling inversion][Count with collisions allowed, then
  remove them: \
  $h(n) = sum_j S(n, j) a(j) <=> a(n) =
  sum_j (-1)^(n - j) c(n, j) h(j)$.]
#fact[Min-max][$max(S) = sum_(emptyset != T subset.eq S)
  (-1)^(|T| + 1) min(T)$, also in expectation (unequal coupons:
  Probability and expectation). The $k$-th largest is
  $sum_(|T| >= k) (-1)^(|T| - k) binom(|T| - 1, k - 1) min(T)$.]

=== Symmetry, trees, determinants
#fact[Burnside][Orbits $= 1/(|G|) sum_(g in G) |"Fix"(g)|$; a
  permutation with $c$ cycles fixes $k^c$ colourings. Rotation by
  $i$ of an $n$-cycle has $gcd(i, n)$ cycles, so necklaces of length
  $n$ with $k$ colours number $1/n sum_(d | n) phi(d) k^(n\/d)$.]
#fact[Bracelets][Up to rotation and reflection: necklaces$\/2$ plus
  $k^((n + 1)\/2)\/2$ for odd $n$, $(k + 1) k^(n\/2)\/4$ for even
  $n$. Cube faces: $(k^6 + 3k^4 + 12k^3 + 8k^2)\/24$.]
#fact[Cayley][$n^(n - 2)$ labelled trees, $n^(n - 1)$ rooted.
  Degrees $d_i$: $(n - 2)! \/ product_i (d_i - 1)!$. Forests with $k$
  trees whose roots are $k$ given vertices: $k n^(n - k - 1)$.
  $K_(n, m)$ has $n^(m - 1) m^(n - 1)$ spanning trees.]
#full[#fact[Prüfer][Remove the smallest leaf $n - 2$ times, writing
  its neighbour: a bijection between labelled trees and
  $[n]^(n - 2)$; vertex $v$ appears $deg(v) - 1$ times. Decode: set
  $deg(v) = 1 +$ (occurrences of $v$); for each entry $x$ in order
  join $x$ to the smallest $v$ with $deg(v) = 1$ and lower both
  degrees; the last two vertices of degree 1 form the last edge
  (a min-heap: $O(n log n)$).]]
#fact[LGV][Vertex-disjoint path systems $a_i -> b_(sigma(i))$ in a
  DAG: $det(M)$ with $M_(i j)$ = number of paths $a_i -> b_j$ counts
  them with sign $"sgn"(sigma)$. If only $sigma = "id"$ can be
  disjoint (planar, ends in order), it is the count.]
#fact[Lucas][$binom(n, k) equiv product_i binom(n_i, k_i)$ mod $p$ for
  base-$p$ digits $n_i, k_i$; odd iff `(k & ~n) == 0`. For fixed
  $k < p$, $binom(m, k) equiv binom(m mod p, k)$, so it, and any
  product of such binomials, is $p$-periodic in $m$.]
#fact[Kummer][The exponent of $p$ in $binom(n, k)$ is the number of
  carries when adding $k$ and $n - k$ in base $p$. Legendre: \
  $v_p (n!) = sum_(i >= 1) floor(n\/p^i) = (n - s_p (n))\/(p - 1)$,
  $s_p$ = digit sum.]

=== Small tables
#table(
  columns: 3,
  [$n$], [$c(n, k)$ for $k = 1..n$], [$S(n, k)$ for $k = 1..n$],
  [1], [1], [1],
  [2], [1 1], [1 1],
  [3], [2 3 1], [1 3 1],
  [4], [6 11 6 1], [1 7 6 1],
  [5], [24 50 35 10 1], [1 15 25 10 1],
  [6], [120 274 225 85 15 1], [1 31 90 65 15 1],
  [7], [720 1764 1624 735 175 21 1], [1 63 301 350 140 21 1],
)
- Involutions $I_n = I_(n - 1) + (n - 1) I_(n - 2)$: 1, 1, 2, 4, 10,
  26, 76, 232, 764. Ordered set partitions (Fubini)
  $sum_k k! S(n, k)$: 1, 1, 3, 13, 75, 541, 4683. Motzkin (paths
  with flat steps): 1, 1, 2, 4, 9, 21, 51, 127, 323. Narayana
  $N(n, k) = 1/n binom(n, k) binom(n, k - 1)$, Dyck paths with $k$
  peaks: 1; 1 1; 1 3 1; 1 6 6 1; 1 10 20 10 1.
- Unlabelled rooted trees: 1, 1, 2, 4, 9, 20, 48, 115, 286; free
  trees: 1, 1, 1, 2, 3, 6, 11, 23, 47, 106 ($n = 1, 2, dots$).
#full[
- $n!$ for $n = 0..12$: 1, 1, 2, 6, 24, 120, 720, 5040, 40320,
  362880, 3628800, 39916800, 479001600.
- Catalan up to $C_(19) = 1767263190$: 58786, 208012, 742900,
  2674440, 9694845, 35357670, 129644790, 477638700 for $n = 11..18$.
- Bell: 678570, 4213597, 27644437, 190899322, 1382958545 for
  $n = 11..15$. Partitions: $p(200) approx 4.0 dot 10^12$. Larger
  values and the 64-bit limits: Numbers, estimates and numerics.

=== More formulas
- Dissections of a convex $n$-gon by $k$ non-crossing diagonals:
  $1/(k + 1) binom(n - 3, k) binom(n + k - 1, k)$; with restricted
  face sizes see the cycle lemma in Generating functions.
- Binary trees and forests: the left child is the first child and
  the right child the next sibling. Binary preorder is forest
  preorder and binary inorder is forest postorder; in the forest's
  bracket string the $k$-th "(" belongs to the $k$-th vertex in
  preorder and the $k$-th ")" to the $k$-th in postorder.
- Tamari lattice: full binary trees with $ell$ leaves under right
  rotations $((A B) C) -> (A (B C))$. For the gap $g$ between
  leaves $g$ and $g + 1$ let $r_g$ be the last leaf below their
  LCA. A rotation raises exactly one $r_g$, and $T'$ is reachable
  from $T$ iff $r_g (T) <= r_g (T')$ for every $g$.
- Monotone paths that must stay on one side of a monotone
  staircase: divide and conquer over columns. A rectangle the
  staircase does not cut passes the counts on its left and bottom
  edges to its right and top edges by convolutions with binomials,
  $O(N log^2 N)$ in all.
- Monotone paths along a chain of rectangles, each sharing a whole
  side with the next, the shared sides of a rectangle meeting at a
  corner: pass first hits along. If $f_s$ paths first meet a side
  at its point $s$ of $0..h$, numbered along the path, the next
  side, of length $w$, gets $"out"_0 = sum_(s <= h) f_s$ and
  $"out"_t = sum_(s < h) f_s binom(h - 1 - s + t, t)$ for $t >= 1$,
  a binomial kernel (_Polynomial techniques_): near-linear in the
  side lengths even if the area is quadratic. Paths from start to
  end of a staircase that touch every run are the paths in the
  chain spanned by consecutive runs: take $f_0 = 1$ on the first
  run and sum the last vector. A coin that picks one of two
  neighbouring positions picks an odd or an even index, so a play
  that empties an array is a monotone path to $(O, E)$, its sums
  over odd and even indices.
- Colourings of the $m$ cells (or edges) of an $n times n$ torus up
  to translation: $1/n^2 sum_(d | n) J_2 (d) k^(m\/d)$. Exactly
  $J_2 (d)$ shifts have order $d$ (Jordan's totient: Totients in
  Number theory facts), and each splits them into cycles of
  length $d$.
- $sum_(m >= r) binom(m, r) x^m = x^r \/ (1 - x)^(r + 1)$, so
  $sum_m binom(m, r) \/ 2^(m + 1) = 1$ for every $r$. With Newton
  series (Generating functions), a polynomial $g$ has
  $sum_(r >= 0) (Delta^r g)(0) = sum_(m >= 0) g(m) \/ 2^(m + 1)$,
  where $(Delta^r g)(0) = sum_m (-1)^(r - m) binom(r, m) g(m)$: if
  $g(m)$ counts objects on $m$ given slots, it counts those that
  use all $r$ slots. E.g. $sum_k k! S(n, k) = sum_m m^n \/ 2^(m + 1)$.
- Kummer corollaries: $binom(2n, n)$ is even for $n >= 1$;
  $binom(n, k)$ is odd for all $k$ iff $n = 2^j - 1$.
- Lucas for prime powers needs Granville's generalization; for a
  composite modulus, compute modulo each prime power and use CRT.

=== Recognition
- "Every value/row/colour appears": PIE over missing ones, or
  surjections $k! S(n, k)$.
- "Exactly $k$ of ...": force $k$, count, binomial inversion.
- "Never below zero", "prefix balance", "first return":
  reflection or Catalan; a step sequence with steps $>= -1$
  summing to $-1$ has exactly one valid rotation (cycle lemma).
  "Maximum subarray sum" of a $plus.minus 1$ sequence: sticky walk.
- "Up to rotation", "distinct colourings": Burnside.
- "Non-intersecting paths", "disjoint routes in a grid": LGV.
- "Value at position $i$ at most $b_i$", "$max_i i a_i <= k$":
  position bounds. LIS or LDS over all permutations or matchings:
  RSK and hook lengths.
- gcd conditions: Möbius over multiples; mod a small prime with
  $n >= p$: Lucas.
- Answer sequence for small $n$ looks familiar: compare with the
  tables above before deriving anything.
]
