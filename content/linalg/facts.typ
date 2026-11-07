#import "/typst/notes.typ": *

#full[
Over a field ($FF_p$, $FF_2$, $RR$) every question below is
Gaussian elimination: `linalg/gauss-mod`, `linalg/gauss-real`,
`linalg/xor-basis`. Work in the field the problem lives in: a
matrix can have full rank over $QQ$ but not mod $p$, and a 0/1
matrix can have different ranks over $FF_2$ and over $RR$.
]

=== Matrices and subspaces
#full[
#fact[Solutions][For an $n times m$ matrix
  $"rank" A + dim ker A = m$. $A x = b$ is solvable iff
  $"rank" A = "rank" [A | b]$; the solutions are one solution plus
  $ker A$, so there are $p^(m - "rank" A)$ of them over $FF_p$.
  "Count assignments mod $p$ (or 0/1 under parity rules) that
  satisfy linear conditions" is this formula. A basis of $ker A$:
  one vector per free column (`linalg/gauss-mod` Example).]
#fact[Exchange][Replacing $v_i$ in a basis by $w$ gives a basis iff
  the $v_i$-coordinate of $w$ is nonzero. Over $FF_2$, replacing
  $v_i$ by $v_i xor x$ keeps a basis iff the $v_i$-coordinate of
  $x$ is 0, so a rank oracle reveals coordinates one by one.]
]
#fact[Random matrices][A uniform $n times n$ matrix over $FF_q$ is
  invertible with probability $product_(i = 1)^n (1 - q^(-i))$:
  about $0.289$ over $FF_2$, about $1 - 1\/q$ for large $q$.]
#fact[Subspaces][$FF_q^n$ has $product_(i < k) (q^(n - i) - 1) \/
  (q^(k - i) - 1)$ subspaces of dimension $k$, so a plane has
  $q + 1$ lines.#full[ If $q + 1$ pairwise non-proportional linear
  forms span a 2-dimensional space of forms, every vector is a
  zero of at least one of them (with basis $f, g$ they are $f$ and
  $g + t f$).]]
#full[
- Intersection of spans $U$ and $W$ (Zassenhaus): eliminate the
  rows $(u, u)$ for a basis of $U$ and $(w, 0)$ for a basis of $W$;
  the echelon rows of the form $(0, x)$ are a basis of
  $U inter W$.
- Least squares for an overdetermined or noisy $A x approx b$:
  solve $A^T A x = A^T b$ (`linalg/gauss-real`). Measurements
  known only up to sign (unsigned areas that are linear in the
  unknowns): try every sign pattern, keep the smallest residual.
]
#fact[Tridiagonal][$a_i x_(i - 1) + b_i x_i + c_i x_(i + 1) = d_i$
  in $O(n)$ (Thomas): forwards
  $c'_i = c_i \/ (b_i - a_i c'_(i - 1))$ and
  $d'_i = (d_i - a_i d'_(i - 1)) \/ (b_i - a_i c'_(i - 1))$, then
  $x_i = d'_i - c'_i x_(i + 1)$ backwards. Safe without pivoting
  when $|b_i| > |a_i| + |c_i|$; mod $p$ it only needs nonzero
  pivots.#full[ Expectations on a line or level by level are such
  systems.]]

=== Determinants
#full[
- Matrix determinant lemma: $det(A + u v^T) = det A (1 + v^T
  A^(-1) u)$; so $det(I + u v^T) = 1 + v^T u$, and every
  principal minor of $I + bold(1) a^T$ is 1 plus the sum of the
  chosen $a_i$.
- Vandermonde: $det(x_i^j)_(0 <= i, j < n) = product_(i < j)
  (x_j - x_i)$, nonzero iff the $x_i$ are distinct.
]
- The permanent is hard in general, but mod 2 it equals the
  determinant.
- Composite modulus: clear each column with Euclid steps between
  two rows (subtract $floor(a\/b)$ times, swap, repeat; each swap
  flips the sign) in $O(n^3 log M)$.
#full[
#fact[Traces][$op("tr")(A^k) = sum_i lambda_i^k$ are the power sums
  $p_k$ of the eigenvalues. Newton's identities (Generating
  functions) turn $p_1, dots, p_n$ into the elementary symmetric
  sums $e_k$; then $det(t I - A) = sum_k (-1)^k e_k t^(n - k)$ and
  $det A = e_n$ (divisions by $k <= n$: fine modulo a prime
  $> n$). So "sums of cyclic products over all closed walks" give
  determinants.]
- Graham-Pollak: the distance matrix of any tree on $n >= 2$
  vertices has determinant $(-1)^(n - 1) (n - 1) 2^(n - 2)$.
]

=== XOR and $FF_2$
- #full[An XOR basis is Gaussian elimination on the bit vectors of
  the numbers. ]$n$ numbers of rank $r$ give $2^r$ distinct subset
  XORs; each of them, 0 included, is the XOR of exactly
  $2^(n - r)$ subsets. Non-empty subsets with XOR 0:
  $2^(n - r) - 1$. A fixed element lies in none or in exactly
  $2^(n - r - 1)$ of them.
#full[
- Which elements lie in some zero-XOR subset: store with each
  basis vector the mask of independent inserted elements it is
  the XOR of. An element that reduces to 0 forms a zero-XOR subset
  with the mask it reduced by; these circuits span all zero-XOR
  subsets, so OR their masks.
]
- Masking is linear: the largest `(subset XOR) & mask` is the
  largest subset XOR of the values `a_i & mask`.
#full[
- Walks in a graph where the length is the XOR of the edges and
  edges may repeat: take root distances $d$ in a spanning tree;
  every non-tree edge $(u, v, w)$ adds $d_u xor d_v
  xor w$ to a basis, and the reachable values from $s$ to
  $t$ are $d_s xor d_t$ plus the span.
- Cutting an array into segments: a segment XOR is
  $P_j xor P_i$ of prefix XORs, and segment XORs are
  independent iff the prefix XORs at the cuts are.
- Products of numbers that are perfect squares: vectors of prime
  exponents mod 2; a product is a square iff the XOR is 0.
- Toggle puzzles (Lights Out): a button is 0/1, a cell gives one
  equation. Chase the lights row by row so that only the first
  row is unknown, then eliminate with bitsets in $O(n^3 \/ 64)$.
- Conditions on rectangle sums mod 2: take the 2D prefix parities
  $P(x, y)$ as unknowns (a bijection with the grid); a rectangle
  sum is the XOR of four of them. Dually, flipping a rectangle
  toggles four cells of the 2D difference array.
- Two bases span the same space iff their reduced row echelon
  forms are equal: hash that form to compare spans. For a fixed
  right end `r`, the spans of `a[l..r]` take at most $B + 1$
  distinct values over all `l` (`PrefixXorBasis`).
]

#full[
=== Rank tests with random values
Schwartz-Zippel, Freivalds and the Tutte matrix are in
Randomisation and hashing; two disjoint bases (union of linear
matroids) by random lifting in Matroids.
- Perfect matching in a bipartite graph: random $x_(u v)$ in the
  $n times n$ Edmonds matrix; a matching exists iff $det != 0$
  (w.h.p.).
- A random substitution can only lower a rank: repeat with fresh
  values and keep the maximum.
- Graphs are linear: edge $u v$ is the vector $e_u - e_v$ (any
  field), so an edge set has rank $n$ minus its number of
  components.
#fact[Exterior algebra][$v_1, dots, v_d in FF^n$ are independent
  iff $v_1 and dots and v_d != 0$. A wedge of $d$ vectors lives on
  the coordinates $e_S$, $|S| = d$; $e_S and e_u = (-1)^c
  e_(S union {u})$ for $u in.not S$, with $c$ the number of
  elements of $S$ above $u$, and $e_S and e_u = 0$ for $u in S$. So
  the wedge product $product_i (1 + r_i v_i)$ with random $r_i$ mod
  a large prime has a nonzero part of degree $d$ iff the $v_i$
  have rank $>= d$ (w.h.p.). It sums over all subsets at $O(2^n)$
  per edge vector $e_u - e_v$ ($O(n 2^n)$ per dense vector):
  "fewest groups whose edges reach rank $d$" for
  $n <= 16$ is a DP over groups of such sums.]

=== Counting with determinants
- Spanning trees and arborescences are cofactors of a Laplacian
  (matrix-tree theorem, Graph theory facts).
- LGV lemma (Counting formulas), example: pairs of monotone grid
  paths $(0, 0) -> (n, m)$ that meet only at the ends are disjoint
  paths from $(0, 1), (1, 0)$ to $(n - 1, m), (n, m - 1)$. With
  $P(a, b) = binom(a + b, a)$ paths of $a$ right and $b$ up steps,
  the count is $P(n - 1, m - 1)^2 - P(n, m - 2) P(n - 2, m)$.
  Weighted paths work the same way (products of edge weights); the
  $k times k$ matrix of path counts comes from binomials or a DP.
]

=== Matrix powers
#fact[Cayley-Hamilton][$p(A) = 0$ for the characteristic
  polynomial $p(lambda) = det(lambda I - A)$ of degree $n$. So
  $A^k = r(A)$ for $r = x^k mod p$, and $s_k = u^T A^k v$ satisfies
  a linear recurrence of order at most $n$: get $2n + 10$ terms by
  simulation and use `poly/linear-rec`, $O(n^2 log k)$ instead of
  $O(n^3 log k)$.#full[ Walk counts, (max, +) powers and
  $sum_(i < k) A^i$: `linalg/matrix`.]]
