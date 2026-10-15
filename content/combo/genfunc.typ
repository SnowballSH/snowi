#import "/typst/notes.typ": *

A generating function packs a sequence into a power series;
$[x^n] A$ is the coefficient of $x^n$. Model the objects, write an
equation between series, then extract coefficients with `poly/`.
Everything is modulo a prime, so division by $n$ is fine for
$n < p$. $S$ and $c$ are the Stirling numbers of Counting
formulas.

=== OGF or EGF
#fact[OGF][$A(x) = sum_n a_n x^n$ for unlabelled objects glued in
  order: strings, compositions, plane trees, dissections.
  $[x^n] A B = sum_i a_i b_(n - i)$ counts ordered pairs whose
  sizes add to $n$.]
#fact[EGF][$A(x) = sum_n a_n x^n\/n!$ for labelled objects whose $n$
  atoms carry the labels $1..n$: permutations, set partitions,
  labelled graphs. $n! [x^n] A B = sum_i binom(n, i) a_i b_(n - i)$:
  the product also chooses which labels go left.]
#table(
  columns: 3,
  [construction ($A(0) = 0$)], [OGF], [EGF],
  [union, ordered pair], [$A + B$, $A B$], [$A + B$, $A B$],
  [sequence of pieces], [$1\/(1 - A)$], [$1\/(1 - A)$],
  [set of pieces], [$product_k (1 - x^k)^(-a_k)$], [$exp(A)$],
  [set of exactly $k$], [], [$A^k\/k!$],
  [cycle of pieces], [], [$ln 1\/(1 - A)$],
  [substitute $B$ for atoms], [$A(B(x))$], [$A(B(x))$],
)

=== Common series
- $1\/(1 - a x) = sum_n a^n x^n$, $(1 + x)^m = sum_k binom(m, k) x^k$,
  \ $1\/(1 - x)^k = sum_n binom(n + k - 1, k - 1) x^n$.
- $x\/(1 - x)^2 = sum_n n x^n$, $x (1 + x)\/(1 - x)^3 = sum_n n^2 x^n$.
- $ln 1\/(1 - x) = sum_(n >= 1) x^n\/n$; Fibonacci $x\/(1 - x - x^2)$.
- Catalan numbers: $(1 - sqrt(1 - 4x))\/(2x)$; central binomials:
  $1\/sqrt(1 - 4x) = sum_n binom(2n, n) x^n$.
- EGFs: Bell $e^(e^x - 1)$, $S(dot, k)$: $(e^x - 1)^k\/k!$,
  $c(dot, k)$: $(ln 1\/(1 - x))^k\/k!$, derangements $e^(-x)\/(1 - x)$,
  involutions $e^(x + x^2\/2)$, ordered set partitions $1\/(2 - e^x)$,
  Bernoulli numbers $x\/(e^x - 1)$.
- $sum_n S(n, k) x^n = x^k \/ product_(i = 1)^k (1 - i x)$;
  $sum_k c(n, k) x^k = x (x + 1) dots (x + n - 1)$ (divide and
  conquer products for a whole row or column).
- Partitions $product_(k >= 1) 1\/(1 - x^k)$; distinct parts
  $product_k (1 + x^k)$; parts from a set $T$:
  $product_(k in T) 1\/(1 - x^k)$. Pentagonal theorem:
  $product_(k >= 1) (1 - x^k) = sum_(j in ZZ) (-1)^j x^(j (3j - 1)\/2)$, so
  $p(0..n)$ is one series inverse.

=== Operations
- Shifts: $x^m A$ gives $a_(n - m)$; $(A - a_0)\/x$ gives $a_(n + 1)$.
  Prefix sums: $A\/(1 - x)$. Differences: $(1 - x) A$.
- $x A'$ gives $n a_n$ (OGF); for an EGF, $A'$ gives $a_(n + 1)$.
  $A(c x)$ gives $c^n a_n$; $A(x^k)$ spreads the terms.
- Product = convolution (`poly/ntt`); binomial convolution
  $sum_i binom(n, i) a_i b_(n - i)$ = EGF product.
- Binomial transform: $b_n = sum_k binom(n, k) a_k$ iff
  $B = e^x A$ as EGFs; invert with $e^(-x)$.
- Even part $(A(x) + A(-x))\/2$. Sum of $a_n$ over $n equiv r$
  (mod $k$): $1/k sum_(j < k) omega^(-j r) A(omega^j)$ with $omega$ a
  primitive $k$-th root of unity (mod $p$ it exists iff
  $k | p - 1$).
- Taylor shift $A(x + c)$: $n! [x^n] = sum_(m >= n) m! a_m c^(m - n)
  \/ (m - n)!$, one convolution after reversing.
- Composition $A(B(x))$ with $B(0) = 0$ substitutes a $B$-structure
  for every atom. Naive: Horner on truncated series,
  $O(n^2 log n)$ with NTT.

=== Solving recurrences
- $a_n = sum_(i = 1)^d c_i a_(n - i)$ for $n >= d$ iff $A = P\/Q$ with
  $Q = 1 - c_1 x - dots - c_d x^d$ and $P = A Q mod x^d$. The $n$-th
  term: `poly/linear-rec` (Berlekamp-Massey finds $c$ from $2d$
  terms).
- Closed form: if $z^d - c_1 z^(d - 1) - dots - c_d$ has roots $r_i$
  with multiplicities $m_i$, then $a_n = sum_i q_i (n) r_i^n$ with
  $deg q_i < m_i$ (partial fractions). Fibonacci:
  $(phi^n - psi^n)\/sqrt(5)$.
- A polynomial term of degree $k$ ($a_n = c a_(n - 1) + f(n)$)
  multiplies $Q$ by $(1 - x)^(k + 1)$: still linear, order
  $d + k + 1$.
- Coefficients of $Q = P^N$ for a short $P$ with $p_0 != 0$ (else
  factor out $x^v$ first): $P Q' = N P' Q$ gives $q_0 = p_0^N$ and
  $n p_0 q_n = sum_(i = 1)^d ((N + 1) i - n) p_i q_(n - i)$ with
  $d = deg P$. That is $O(d)$ per coefficient with $1\/n$
  precomputed (indices $n < p$); no NTT, good for $N = 10^7$.

=== Exponential formula
#fact[Components][If every labelled structure is a set of
  connected ones: $"All" = exp("Conn")$, $"Conn" = ln "All"$ (EGFs,
  $"All"(0) = 1$). Exactly $k$ components: $"Conn"^k\/k!$. Total
  number of components over all structures: $"Conn" dot "All"$;
  the expected number is $[x^n] "Conn" "All" \/ [x^n] "All"$.]
#fact[Without NTT][Fix the component of label 1:
  $a_n = sum_(k = 1)^n binom(n - 1, k - 1) c_k a_(n - k)$, so
  $c_n = a_n - sum_(k = 1)^(n - 1) binom(n - 1, k - 1) c_k a_(n - k)$
  in $O(n^2)$.]
- Connected labelled graphs: $"All" = sum_n 2^(binom(n, 2)) x^n\/n!$;
  $c = 1, 1, 4, 38, 728$.
- Permutations are sets of cycles: $1\/(1 - x) = exp(ln 1\/(1 - x))$;
  set partitions: $exp(e^x - 1)$; functional graphs: sets of
  cycles of rooted trees.
- Unlabelled multisets (Euler transform): \
  $product_k (1 - x^k)^(-a_k) = exp(sum_(j >= 1) A(x^j)\/j)$; e.g.
  unlabelled rooted trees satisfy $T = x exp(sum_j T(x^j)\/j)$.

=== Lagrange inversion
#fact[Theorem][If $F = x R(F)$ with $R(0) != 0$, then for
  $n, k >= 1$: $[x^n] F^k = k/n [t^(n - k)] R(t)^n$, and
  $[x^n] H(F) = 1/n [t^(n - 1)] H'(t) R(t)^n$. If instead the
  compositional inverse $G$ of $F$ is known, $R(t) = t\/G(t)$.]
- Binary trees $B = x (1 + B)^2$: $[x^n] B = 1/n binom(2n, n - 1) =
  C_n$. Trees with child counts from $D$: $R(t) = sum_(d in D) t^d$;
  $k$-ary trees: $1/((k - 1) n + 1) binom(k n, n)$.
- Rooted labelled trees $T = x e^T$: $n! [x^n] T = n^(n - 1)$;
  forests of $k$ rooted trees: $binom(n - 1, k - 1) n^(n - k)$.
- Cycle lemma: integers $>= -1$ summing to $-1$ have exactly one
  rotation whose proper prefix sums are all $>= 0$. Write a plane
  tree in preorder as (children $- 1$) per node: trees with a
  given multiset of child counts number (arrangements)$\/n$.
- "For every $k$" ($[x^n] F^k$ for all $k$): one power $R^n$ gives
  all of them.

=== Toolbox
- Products: `poly/ntt` ($O(n log n)$); many small factors: divide
  and conquer, $O(n log^2 n)$.
- $1\/A$, $ln A$, $exp A$, $A^k$ mod $x^n$: `poly/poly-ops`, Newton,
  $O(n log n)$. $ln$ needs $a_0 = 1$, $exp$ needs $a_0 = 0$.
- $n$-th term of a rational series: `poly/linear-rec`. Answer is a
  polynomial in a parameter: `poly/lagrange`. OR, AND, XOR
  convolutions: `poly/subset-transforms`.
- $O(n^2)$ is fine up to $n approx 5000$: inverse
  $b_n = -sum_(i = 1)^n a_i b_(n - i)$ ($a_0 = 1$); $E = exp A$ from
  $n e_n = sum_(k = 1)^n k a_k e_(n - k)$; $L = ln A$ from
  $n a_n = sum_(k = 1)^n k l_k a_(n - k)$.

=== Pitfalls
- Labelled vs unlabelled gives plausible wrong numbers: check
  $n <= 4$ by brute force.
- Lagrange uses the power $n$ and coefficient $n - k$; $F = x^2 dots$
  ($R(0) = 0$) does not apply directly.
- Constants: $a^(binom(n, 2))$ reduces the exponent mod $p - 1$, but
  $k dot ln A$ multiplies by $k$ mod $p$.
- Truncate every product to the needed length; resize after each
  multiplication.
- Expected values are (sum over objects)$\/$(number of objects):
  extract both coefficients, divide once at the end.
- Cues: "count objects of size $n$ mod 998244353" with a recursive
  shape; "for each $k$"; "connected"; "expected number of
  components"; "prefix stays non-negative" (cycle lemma); $n$ near
  $10^7$ means an $O(n)$ recurrence, not NTT.
