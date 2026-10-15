#import "/typst/notes.typ": *

Over a field ($FF_p$, $FF_2$, $RR$) every question below is
Gaussian elimination: `linalg/gauss-mod`, `linalg/gauss-real`,
`linalg/xor-basis`. Work in the field the problem lives in: a
matrix can have full rank over $QQ$ but not mod $p$, and a 0/1
matrix has different ranks over $FF_2$ and over $RR$.

=== Rank and solutions
#fact[Rank-nullity][For an $n times m$ matrix,
  $"rank" A + dim ker A = m$, and row rank = column rank
  $<= min(n, m)$. The rank is the number of pivots.]
#fact[Solutions][$A x = b$ is solvable iff
  $"rank" A = "rank" [A | b]$; the solutions are one solution plus
  $ker A$, so there are $p^(m - "rank" A)$ of them over $FF_p$.
  "Count assignments mod $p$ (or 0/1 under parity rules) that
  satisfy linear conditions" is this formula. A basis of $ker A$:
  one vector per free column (`linalg/gauss-mod` Example).]
#fact[Square][For $n times n$ matrices these are equivalent:
  $det A != 0$, rank $n$, $A$ invertible, $A x = b$ has exactly one
  solution for every $b$, $A x = 0$ only for $x = 0$.]
#fact[Random matrices][A uniform $n times n$ matrix over $FF_q$ is
  invertible with probability $product_(i = 1)^n (1 - q^(-i))$:
  about $0.289$ over $FF_2$, about $1 - 1\/q$ for large $q$.]

=== Determinants
- $det(A B) = det A det B$, $det A^T = det A$,
  $det(c A) = c^n det A$, $det A^(-1) = 1 \/ det A$.
- A row swap flips the sign, adding a multiple of one row to
  another keeps it; a triangular or block-triangular matrix gives
  the product of the diagonal entries (blocks).
- Matrix determinant lemma: $det(A + u v^T) = det A (1 + v^T
  A^(-1) u)$; so $det(I + u v^T) = 1 + v^T u$, and every
  principal minor of $I + bold(1) a^T$ is 1 plus the sum of the
  chosen $a_i$.
- Vandermonde: $det(x_i^j)_(0 <= i, j < n) = product_(i < j)
  (x_j - x_i)$, nonzero iff the $x_i$ are distinct.
- The permanent is hard in general, but mod 2 it equals the
  determinant.
- Composite modulus: clear each column with Euclid steps between
  two rows (subtract $floor(a\/b)$ times, swap, repeat; each swap
  flips the sign) in $O(n^3 log M)$.

=== XOR and $FF_2$
- An XOR basis is Gaussian elimination on the bit vectors of the
  numbers. $n$ numbers of rank $r$ give $2^r$ distinct subset
  XORs; each of them, 0 included, is the XOR of exactly
  $2^(n - r)$ subsets. Non-empty subsets with XOR 0:
  $2^(n - r) - 1$. A fixed element lies in none or in exactly half
  of the zero-XOR subsets.
- Masking is linear: the largest `(subset XOR) & mask` is the
  largest subset XOR of the values `a_i & mask`.
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
- Two bases span the same space iff their reduced row echelon
  forms are equal: hash that form to compare spans. For a fixed
  right end `r`, the spans of `a[l..r]` take at most $B + 1$
  distinct values over all `l` (`PrefixXorBasis`).

=== Randomized checks
#fact[Schwartz-Zippel][A nonzero polynomial of total degree $d$
  over a field is zero at a uniform random point of $S^k$ with
  probability at most $d \/ |S|$. Work mod a prime near $10^9$
  (or 998244353), not over $FF_2$.]
- Freivalds: $A B = C$ iff $A (B r) = C r$ for random vectors $r$;
  $O(n^2)$ per try, error at most $1\/p$ per try mod $p$.
- Perfect matching in a bipartite graph: random $x_(u v)$ in the
  Edmonds matrix; a matching exists iff $det != 0$ (w.h.p.).
  General graph: Tutte matrix $T_(u v) = x_(u v) = -T_(v u)$ for
  each edge; $"rank" T$ is twice the maximum matching size.
- Rank of hard combinations (two disjoint bases, union of linear
  matroids): replace $v_i$ by $(alpha_i v_i, beta_i v_i)$ with
  random $alpha_i, beta_i$ and take an ordinary rank.
- A random substitution can only lower a rank: repeat with fresh
  values and keep the maximum. Polynomial string hashing is the
  one-variable case.

=== Counting with determinants
#fact[Matrix-tree][Undirected multigraph: the number of spanning
  trees is any cofactor of the Laplacian $L = D - A$ (delete row
  and column $v$, take the determinant). Directed: let $L_(v v)$ be
  the in-degree of $v$ and $L_(u v) = -c(u, v)$ for $u != v$, where
  $c(u, v)$ counts the edges $u -> v$. Deleting row and column $r$
  gives the number of spanning arborescences with every edge
  directed away from the root $r$; out-degrees on the diagonal
  count those directed towards $r$. Self-loops never count. Take
  the determinant mod $p$ with `linalg/gauss-mod`.]
#fact[LGV lemma][In a DAG with sources $a_1, dots, a_k$ and sinks
  $b_1, dots, b_k$, let $M_(i j)$ be the number (or total weight)
  of paths $a_i -> b_j$. Then $det M$ sums $"sgn"(sigma)$ over all
  systems of vertex-disjoint paths $a_i -> b_(sigma(i))$. When only
  $sigma = "id"$ can be disjoint (planar graph, sources and sinks
  in matching order on the boundary, such as grid paths going
  right and up), $det M$ counts the non-intersecting systems.]
- Example: pairs of monotone grid paths $(0, 0) -> (n, m)$ that
  meet only at the ends are disjoint paths from $(0, 1), (1, 0)$
  to $(n - 1, m), (n, m - 1)$. With $P(a, b) = binom(a + b, a)$
  paths of $a$ right and $b$ up steps, the count is
  $P(n - 1, m - 1)^2 - P(n, m - 2) P(n - 2, m)$.
- Weighted counts work the same way (products of edge weights);
  $k$ is usually small, so $M$ comes from binomials or a DP.

=== Matrix powers and recurrences
#fact[Cayley-Hamilton][$p(A) = 0$ for the characteristic
  polynomial $p(lambda) = det(lambda I - A)$ of degree $n$. So
  $s_k = u^T A^k v$ satisfies a linear recurrence of order at most
  $n$: get $2n + 10$ terms by simulation, find the recurrence with
  `poly/linear-rec` (Berlekamp-Massey) and evaluate $s_k$ in
  $O(n^2 log k)$ instead of $O(n^3 log k)$. Also $A^k = r(A)$ for
  $r = x^k mod p$.]
- Walks with exactly $k$ steps: $A^k$ (`linalg/matrix`); shortest
  or longest walks: the same power over (min, +) or (max, +).
- $I + A + dots + A^(k - 1)$: the top-right block of
  $mat(A, I; 0, I)^k$.
