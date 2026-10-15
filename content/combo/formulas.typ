#import "/typst/notes.typ": *

Compute values with `combo/binom`; $binom(n, k) = 0$ unless
$0 <= k <= n$. $S(n, k)$ and $c(n, k)$ are Stirling numbers of the
second and (unsigned) first kind, $x^(underline(k)) = x (x - 1) dots
(x - k + 1)$ and $x^(overline(k)) = x (x + 1) dots (x + k - 1)$.

=== Binomials
- $binom(n, k) = binom(n, n - k) = binom(n - 1, k - 1) + binom(n - 1, k)$,
  $k binom(n, k) = n binom(n - 1, k - 1)$, \
  $binom(n, m) binom(m, k) = binom(n, k) binom(n - k, m - k)$.
- Vandermonde: $sum_i binom(m, i) binom(n, k - i) = binom(m + n, k)$;
  squares: $sum_k binom(n, k)^2 = binom(2n, n)$.
- Hockey stick: $sum_(i = r)^n binom(i, r) = binom(n + 1, r + 1)$,
  $sum_(i = 0)^k binom(n + i, i) = binom(n + k + 1, k)$.
- $sum_k binom(n, k) = 2^n$, $sum_k (-1)^k binom(n, k) = [n = 0]$,
  $sum_k k binom(n, k) = n 2^(n - 1)$, \
  $(1 - x)^(-n) = sum_k binom(n + k - 1, k) x^k$.
- Stars and bars for $k >= 1$ ($n = k = 0$ has one solution):
  $x_1 + dots + x_k = n$ has $binom(n + k - 1, k - 1)$ solutions
  with all $x_i >= 0$ and $binom(n - 1, k - 1)$ with all
  $x_i >= 1$. With $x_i <= c$: \
  $sum_j (-1)^j binom(k, j) binom(n - j (c + 1) + k - 1, k - 1)$.
  Grid paths $(0, 0) -> (a, b)$: $binom(a + b, a)$.

=== Paths and Catalan numbers
#fact[Reflection][$plus.minus 1$ walks of $n$ steps from $a$ to $b$:
  $W(b - a)$ with $W(d) = binom(n, (n + d)\/2)$ (0 if parity
  differs). For $a, b >= 0$, those that never go below 0:
  $W(b - a) - W(a + b + 2)$ (reflect the start across $-1$). Grid
  paths to $(n, m)$, $m <= n$, never more ups than rights:
  $binom(n + m, m) - binom(n + m, m - 1)$.]
#fact[Ballot][$A$ gets $a$ votes, $B$ gets $b < a$: $A$ stays strictly
  ahead during a random count with probability $(a - b)\/(a + b)$.]
#fact[Catalan][1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796.
  $C_n = binom(2n, n)\/(n + 1)$ $= binom(2n, n) - binom(2n, n + 1)$,
  $C_(n + 1) = sum_(i = 0)^n C_i C_(n - i)$, and
  $C_(n + 1) = (4n + 2)\/(n + 2) dot C_n$. Counts balanced bracket
  strings with $n$ pairs, Dyck paths, binary trees with $n$ nodes, plane trees with
  $n + 1$ nodes, full binary trees with $n + 1$ leaves,
  triangulations of an $(n + 2)$-gon, non-crossing matchings of
  $2n$ points, non-crossing partitions of $n$, 231-avoiding
  (stack-sortable) permutations.]

=== Stirling, Bell, partitions
#fact[First kind][$c(n, k)$: permutations of $n$ with $k$ cycles.
  $c(n + 1, k) = n c(n, k) + c(n, k - 1)$;
  $x^(overline(n)) = sum_k c(n, k) x^k$; $c(n, 1) = (n - 1)!$.]
#fact[Second kind][$S(n, k)$: partitions of $n$ labelled items into
  $k$ non-empty unlabelled blocks. $S(n + 1, k) = k S(n, k) +
  S(n, k - 1)$; $x^n = sum_k S(n, k) x^(underline(k))$;
  $S(n, k) = 1/(k!) sum_i (-1)^(k - i) binom(k, i) i^n$ (one row by
  convolution). Surjections onto $k$ labelled boxes: $k! S(n, k)$.]
#fact[Bell][$B_n = sum_k S(n, k)$: all set partitions.
  $B_(n + 1) = sum_k binom(n, k) B_k$. 1, 1, 2, 5, 15, 52, 203, 877,
  4140, 21147, 115975.]
#fact[Partitions][$p(n)$: 1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42;
  $p(100) approx 1.9 dot 10^8$. Pentagonal recurrence in
  $O(n sqrt(n))$: $p(n) = sum_(k >= 1) (-1)^(k + 1) (p(n - g_k) +
  p(n - g_(-k)))$, $g_k = k (3k - 1)\/2$, $p(m < 0) = 0$. Exactly
  $k$ parts:
  $p_k (n) = p_(k - 1) (n - 1) + p_k (n - k)$. Distinct parts and
  odd parts are equinumerous.]
#fact[Derangements][$D_n = (n - 1)(D_(n - 1) + D_(n - 2))$
  $= n D_(n - 1) + (-1)^n approx n!\/e$, starting 1, 0, 1, 2, 9, 44,
  265, 1854, 14833. Exactly $k$ fixed points: $binom(n, k) D_(n - k)$.]

=== Twelvefold way
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
)

=== Inclusion-exclusion and inversion
#fact[PIE][Objects in no bad event: $sum_(S) (-1)^(|S|) g(S)$ where
  $g(S)$ counts objects in every event of $S$ (force $S$, leave the
  rest free). If $g(S) = G_(|S|)$: $sum_k (-1)^k binom(n, k) G_k$.]
#fact[Binomial][$g_k = sum_(j >= k) binom(j, k) f_j <=> f_k =
  sum_(j >= k) (-1)^(j - k) binom(j, k) g_j$. "Choose $k$, force
  them" gives $g_k$, not "at least $k$". For all $k$ at once:
  $k! f_k = sum_j j! g_j dot (-1)^(j - k)\/(j - k)!$, one
  convolution after reversing.]
#fact[Subsets][$g(S) = sum_(T subset.eq S) f(T) <=> f(S) =
  sum_(T subset.eq S) (-1)^(|S| - |T|) g(T)$; $O(2^n n)$ with
  `poly/subset-transforms`.]
#fact[Divisors][From $sum_(d | n) mu(d) = [n = 1]$:
  if $g(n) = sum_(d | n) f(d)$ then $f(n) = sum_(d | n) mu(n\/d) g(d)$.
  Over multiples, $g(d) = sum_(d | m) f(m)$ gives
  $f(d) = sum_(d | m) mu(m\/d) g(m)$: "gcd exactly $d$" from "all
  divisible by $d$". Without $mu$: $f(d) = g(d) -
  sum_(m = 2d, 3d, dots) f(m)$ for $d$ descending, $O(V log V)$.]
#fact[Stirling inversion][Count with collisions allowed, then
  remove them: \
  $h(n) = sum_j S(n, j) a(j) <=> a(n) =
  sum_j (-1)^(n - j) c(n, j) h(j)$.]
#fact[Min-max][$max(S) = sum_(emptyset != T subset.eq S)
  (-1)^(|T| + 1) min(T)$, also in expectation: the expected time to
  collect all coupons is $sum_T (-1)^(|T| + 1) \/ P(T)$ with
  $P(T) = sum_(i in T) p_i$.]

=== Symmetry, trees, determinants
#fact[Burnside][Orbits $= 1/(|G|) sum_(g in G) |"Fix"(g)|$; a
  permutation with $c$ cycles fixes $k^c$ colourings. Rotation by
  $i$ of an $n$-cycle has $gcd(i, n)$ cycles, so necklaces of length
  $n$ with $k$ colours number $1/n sum_(d | n) phi(d) k^(n\/d)$.]
#fact[Cayley][$n^(n - 2)$ labelled trees, $n^(n - 1)$ rooted.
  Degrees $d_i$: $(n - 2)! \/ product_i (d_i - 1)!$. Forests with $k$
  trees whose roots are $k$ given vertices: $k n^(n - k - 1)$.
  $K_(n, m)$ has $n^(m - 1) m^(n - 1)$ spanning trees.]
#fact[Prüfer][Remove the smallest leaf $n - 2$ times, writing its
  neighbour: a bijection between labelled trees and $[n]^(n - 2)$;
  vertex $v$ appears $deg(v) - 1$ times.]
#fact[Kirchhoff][Spanning trees = any cofactor of the Laplacian
  $D - A$ (multi-edges allowed); determinant by `linalg/gauss-mod`.]
#fact[LGV][Vertex-disjoint path systems $a_i -> b_(sigma(i))$ in a
  DAG: $det(M)$ with $M_(i j)$ = number of paths $a_i -> b_j$ counts
  them with sign $"sgn"(sigma)$. If only $sigma = "id"$ can be
  disjoint (planar, ends in order), it is the count.]
#fact[Lucas][$binom(n, k) equiv product_i binom(n_i, k_i)$ mod $p$ for
  base-$p$ digits $n_i, k_i$; odd iff `(k & ~n) == 0`.]
#fact[Kummer][The exponent of $p$ in $binom(n, k)$ is the number of
  carries when adding $k$ and $n - k$ in base $p$. Legendre: \
  $v_p (n!) = sum_(i >= 1) floor(n\/p^i) = (n - s_p (n))\/(p - 1)$,
  $s_p$ = digit sum.]

#full[
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
- $n!$ for $n = 0..12$: 1, 1, 2, 6, 24, 120, 720, 5040, 40320,
  362880, 3628800, 39916800, 479001600. $20! approx 2.4 dot 10^18$
  is the last that fits `ll`.
- Catalan up to $C_(19) = 1767263190$: 58786, 208012, 742900,
  2674440, 9694845, 35357670, 129644790, 477638700 for $n = 11..18$.
  $C_(35) approx 3.1 dot 10^18$ is the last in `ll`.
- Bell: 678570, 4213597, 27644437, 190899322, 1382958545 for
  $n = 11..15$; $B_(25) approx 4.6 dot 10^18$ is the last in `ll`.
- Partitions: $p(20) = 627$, $p(50) = 204226$,
  $p(200) approx 4.0 dot 10^12$; $p(405) approx 9.1 dot 10^18$ is the
  last in `ll`.
- Involutions $I_n = I_(n - 1) + (n - 1) I_(n - 2)$: 1, 1, 2, 4, 10,
  26, 76, 232, 764. Ordered set partitions (Fubini)
  $sum_k k! S(n, k)$: 1, 1, 3, 13, 75, 541, 4683. Motzkin (paths
  with flat steps): 1, 1, 2, 4, 9, 21, 51, 127, 323. Narayana rows:
  1; 1 1; 1 3 1; 1 6 6 1; 1 10 20 10 1 (see below).
- Unlabelled rooted trees: 1, 1, 2, 4, 9, 20, 48, 115, 286; free
  trees: 1, 1, 1, 2, 3, 6, 11, 23, 47, 106 ($n = 1, 2, dots$).

=== More formulas
#fact[Eulerian][$A(n, k)$: permutations with $k$ descents.
  $A(n, k) = (k + 1) A(n - 1, k) + (n - k) A(n - 1, k - 1)$; rows 1;
  1 1; 1 4 1; 1 11 11 1; 1 26 66 26 1. Worpitzky:
  $x^n = sum_k A(n, k) binom(x + k, n)$.]
#fact[Narayana][$N(n, k) = 1/n binom(n, k) binom(n, k - 1)$: Dyck
  paths with $k$ peaks; $sum_k N(n, k) = C_n$.]
- Compositions of $n >= 1$ into $k$ positive parts:
  $binom(n - 1, k - 1)$; into any number: $2^(n - 1)$.
- Walks of $n$ steps $plus.minus 1$ from $x$ to $y$ inside $[0, w]$:
  $sum_(j in ZZ) (W(y - x + 2j(w + 2)) - W(-y - 2 - x +
  2j(w + 2)))$, reflecting across $-1$ and $w + 1$ alternately.
- $k$-th largest by min-max:
  $sum_(|T| >= k) (-1)^(|T| - k) binom(|T| - 1, k - 1) min(T)$.
- Bracelets (rotations and reflections): necklaces$\/2$ plus
  $k^((n + 1)\/2)\/2$ for odd $n$, $(k + 1) k^(n\/2)\/4$ for even
  $n$. Cube faces: $(k^6 + 3k^4 + 12k^3 + 8k^2)\/24$.
- $k$-ary trees with $n$ internal nodes:
  $1/((k - 1) n + 1) binom(k n, n)$.
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
- "Up to rotation", "distinct colourings": Burnside.
- "Non-intersecting paths", "disjoint routes in a grid": LGV.
- gcd conditions: Möbius over multiples; mod a small prime with
  $n >= p$: Lucas.
- Answer sequence for small $n$ looks familiar: compare with the
  tables above before deriving anything.
]
