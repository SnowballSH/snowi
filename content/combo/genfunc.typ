#import "/typst/notes.typ": *

#full[
A generating function packs a sequence into a power series;
$[x^n] A$ is the coefficient of $x^n$. Model the objects, write an
equation between series, then extract coefficients (the algorithms
are in Polynomial techniques). Everything is modulo a prime, so
division by $n$ is fine for $n < p$. $S$ and $c$ are the Stirling
numbers of Counting formulas.
]

=== OGF or EGF
#full[
#fact[OGF][$A(x) = sum_n a_n x^n$ for unlabelled objects glued in
  order: strings, compositions, plane trees, dissections.
  $[x^n] A B = sum_i a_i b_(n - i)$ counts ordered pairs whose
  sizes add to $n$.]
#fact[EGF][$A(x) = sum_n a_n x^n\/n!$ for labelled objects whose $n$
  atoms carry the labels $1..n$: permutations, set partitions,
  labelled graphs. $n! [x^n] A B = sum_i binom(n, i) a_i b_(n - i)$:
  the product also chooses which labels go left.]
]
#block(breakable: false, table(
  columns: 3,
  [construction ($A(0) = 0$)], [OGF], [EGF],
  [union, ordered pair], [$A + B$, $A B$], [$A + B$, $A B$],
  [sequence of pieces], [$1\/(1 - A)$], [$1\/(1 - A)$],
  [set of pieces], [$product_k (1 - x^k)^(-a_k)$], [$exp(A)$],
  [set of exactly $k$], [], [$A^k\/k!$],
  [cycle of pieces], [], [$ln 1\/(1 - A)$],
  [substitute $B$ for atoms], [$A(B(x))$], [$A(B(x))$],
))

=== Common series
#full[
- $1\/(1 - x)^k = sum_n binom(n + k - 1, k - 1) x^n$,
  $x (1 + x)\/(1 - x)^3 = sum_n n^2 x^n$,
  $ln 1\/(1 - x) = sum_(n >= 1) x^n\/n$.
- Catalan numbers: $(1 - sqrt(1 - 4x))\/(2x)$; central binomials:
  $1\/sqrt(1 - 4x) = sum_n binom(2n, n) x^n$.
]
- EGFs: Bell $e^(e^x - 1)$, $S(dot, k)$: $(e^x - 1)^k\/k!$,
  $c(dot, k)$: $(ln 1\/(1 - x))^k\/k!$, derangements $e^(-x)\/(1 - x)$,
  involutions $e^(x + x^2\/2)$, ordered set partitions $1\/(2 - e^x)$,
  Bernoulli numbers $x\/(e^x - 1)$.
- $sum_n S(n, k) x^n = x^k \/ product_(i = 1)^k (1 - i x)$;
  $sum_k c(n, k) x^k = x (x + 1) dots (x + n - 1)$ (divide and
  conquer products for a whole row or column).
#full[
- Partitions $product_(k >= 1) 1\/(1 - x^k)$; distinct parts
  $product_k (1 + x^k)$; parts from a set $T$:
  $product_(k in T) 1\/(1 - x^k)$. Pentagonal theorem:
  $product_(k >= 1) (1 - x^k) = sum_(j in ZZ) (-1)^j x^(j (3j - 1)\/2)$, so
  $p(0..n)$ is one series inverse.
]
- Symmetric sums of $a_1, dots, a_m$ with power sums
  $p_k = sum_i a_i^k$: $sum_k e_k t^k = product_i (1 + a_i t) =
  exp(sum_(k >= 1) (-1)^(k - 1) p_k t^k\/k)$; Newton:
  $k e_k = sum_(i = 1)^k (-1)^(i - 1) e_(k - i) p_i$.

=== Operations
- Binomial transform: $b_n = sum_k binom(n, k) a_k$ iff
  $B = e^x A$ as EGFs; invert with $e^(-x)$. So $a_k = (Delta^k
  b)(0)$ and $b(m) = sum_k binom(m, k) (Delta^k b)(0)$ (Newton
  series).
- Sum of $a_n$ over $n equiv r$ (mod $k$)#full[ (even part:
  $(A(x) + A(-x))\/2$)]: $1/k sum_(j < k) omega^(-j r) A(omega^j)$
  with $omega$ a primitive $k$-th root of unity (mod $p$ it exists
  iff $k | p - 1$).#full[ Several such constraints (every row sum
  and column sum $equiv 0$): give each its own root; if a term
  depends only on how many constraints took each root, sum over
  those counts with multinomial weights.]
- Taylor shift $A(x + c)$: $n! [x^n] = sum_(m >= n) m! a_m c^(m - n)
  \/ (m - n)!$, one convolution after reversing.
#full[
- $sum_m b_m binom(m, k - m) = [x^k] B(x + x^2)$. A quadratic
  inside a polynomial $A$ costs two Taylor shifts:
  $A(a x^2 + b x + c) = h(a (x +
  b\/(2a))^2)$ with $h(z) = A(z + c - b^2\/(4a))$: shift, multiply
  $h_i$ by $a^i$, spread to even powers, shift by $b\/(2a)$.
- Composition $A(B(x))$ with $B(0) = 0$ substitutes a $B$-structure
  for every atom. Naive: Horner on truncated series,
  $O(n^2 log n)$ with NTT.
]

=== Solving recurrences
#full[
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
]
- Coefficients of $Q = P^N$ for a short $P$ with $p_0 != 0$ (else
  factor out $x^v$ first): $P Q' = N P' Q$ gives $q_0 = p_0^N$ and
  $n p_0 q_n = sum_(i = 1)^d ((N + 1) i - n) p_i q_(n - i)$ with
  $d = deg P$. That is $O(d)$ per coefficient with $1\/n$
  precomputed (indices $n < p$); no NTT, good for $N = 10^7$.
- The same works whenever $A'\/A$ is rational. For
  $A = e^(a x) product_i P_i^(N_i)$ with $P_i (0) = 1$: $Q A' = R A$
  with $Q = product_i P_i$, $R = a Q + sum_i N_i P'_i Q\/P_i$, so
  $(n + 1) a_(n + 1) = sum_j r_j a_(n - j) -
  sum_(j >= 1) q_j (n + 1 - j) a_(n + 1 - j)$.
#full[
- A recursion in a parameter, $F_(k + 1) = (a F_k + b)\/(c F_k + d)$
  with series $a, b, c, d$: keep $F_k = S_k\/D_k$, then
  $S_(k + 1) = a S_k + b D_k$, $D_(k + 1) = c S_k + d D_k$, a
  $2 times 2$ matrix of series to the $k$-th power (products mod
  $x^(n + 1)$), $O(n log n log k)$; divide once at the end.
]

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
#full[
- Connected labelled graphs: $"All" = sum_n 2^(binom(n, 2)) x^n\/n!$;
  $c = 1, 1, 4, 38, 728$.
- Permutations are sets of cycles: $1\/(1 - x) = exp(ln 1\/(1 - x))$;
  set partitions: $exp(e^x - 1)$; functional graphs: sets of
  cycles of rooted trees.
]
- Unlabelled multisets (Euler transform): \
  $product_k (1 - x^k)^(-a_k) = exp(sum_(j >= 1) A(x^j)\/j)$; e.g.
  unlabelled rooted trees satisfy $T = x exp(sum_j T(x^j)\/j)$.

=== Lagrange inversion
#fact[Theorem][If $F = x R(F)$ with $R(0) != 0$, then for
  $n, k >= 1$: $[x^n] F^k = k/n [t^(n - k)] R(t)^n$, and
  $[x^n] H(F) = 1/n [t^(n - 1)] H'(t) R(t)^n$. If instead the
  compositional inverse $G$ of $F$ is known, $R(t) = t\/G(t)$.]
#fact[Other forms][For $n >= 0$: $[x^n] H(F) = [t^n] H(t) R(t)^(n - 1)
  (R(t) - t R'(t))$ (no derivative of $H$; for $R = e^t$ it is
  $[t^n] H(t) e^(n t) (1 - t)$). All $n$ at once:
  $sum_n x^n [t^n] H(t) R(t)^n = H(F) \/ (1 - x R'(F))$.]
- Trees with child counts from $D$: $R(t) = sum_(d in D) t^d$;
  $k$-ary trees: $1/((k - 1) n + 1) binom(k n, n)$.
- Cycle lemma: integers $>= -1$ summing to $-1$ have exactly one
  rotation whose proper prefix sums are all $>= 0$. Write a plane
  tree in preorder as (children $- 1$) per node: trees with a
  given multiset of child counts number (arrangements)$\/n$.
#full[
- Rooted labelled trees $T = x e^T$: $n! [x^n] T = n^(n - 1)$;
  forests of $k$ rooted trees: $binom(n - 1, k - 1) n^(n - k)$.
- Sum 0: the rotations of a $plus.minus 1$ sequence with sum 0 whose
  prefix sums all stay $>= 0$ are exactly those starting right
  after a minimum of its prefix sums. So a bracket string with as
  many `(` as `)` has a regular rotation.
- Dissections of an $N$-gon are plane trees: a $c$-gon is a node
  with $c - 1$ children, the $N - 1$ sides other than the root side
  are leaves. With $k$ diagonals and face sizes in $S$:
  $1/(N + k) binom(N + k, k + 1) [x^(N - 2)] (sum_(s in S)
  x^(s - 2))^(k + 1)$.
- "For every $k$" ($[x^n] F^k$ for all $k$): one power $R^n$ gives
  all of them. Without a usable $R$ ($F'(0) = 0$) this is power
  projection: $O(n log^2 n)$ by Kinoshita-Li, the transpose of
  composition (no snippet here).
- $G$ without a closed form: if $F$ satisfies a differential
  equation in $u$, put $u = G(x)$ with $F(G) = x$ and
  $F'(G) = 1\/G'$; the result determines $G$ coefficient by
  coefficient (online convolution when it is quadratic in $G$).
]

=== Without NTT
$O(n^2)$ is fine up to $n approx 5000$: inverse
$b_n = -sum_(i = 1)^n a_i b_(n - i)$ ($a_0 = 1$); $E = exp A$ from
$n e_n = sum_(k = 1)^n k a_k e_(n - k)$ ($a_0 = 0$); $L = ln A$ from
$n a_n = sum_(k = 1)^n k l_k a_(n - k)$ ($a_0 = 1$).

=== Pitfalls
#full[
- Labelled vs unlabelled gives plausible wrong numbers: check
  $n <= 4$ by brute force.
]
- Constants: $a^(binom(n, 2))$ reduces the exponent mod $p - 1$, but
  $k dot ln A$ multiplies by $k$ mod $p$.
#full[
- Expected values are (sum over objects)$\/$(number of objects):
  extract both coefficients, divide once at the end.
- Cues: "count objects of size $n$ mod 998244353" with a recursive
  shape; "for each $k$"; "connected"; "expected number of
  components"; "prefix stays non-negative" (cycle lemma); a
  recursion in a value bound $k$ (2 by 2 matrix power); $n$ near
  $10^7$ means an $O(n)$ recurrence, not NTT.
]
