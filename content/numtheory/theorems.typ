#import "/typst/notes.typ": *

Below `p` is a prime, $v_p (n)$ the exponent of `p` in `n`, and
$d(n)$, $sigma(n)$ the number and the sum of the divisors of `n`.

=== Congruences
#fact[Fermat, Euler][$a^(p - 1) equiv 1 thick (mod p)$ if
  $p divides.not a$; $a^(phi(m)) equiv 1 thick (mod m)$ if
  $gcd(a, m) = 1$. For any `a` and $e >= log_2 m$:
  $a^e equiv a^(e mod phi(m) + phi(m)) thick (mod m)$; apply it
  recursively to power towers.]
#fact[Wilson][$(p - 1)! equiv -1 thick (mod p)$, while
  $(n - 1)! equiv 0 thick (mod n)$ for composite $n > 4$.]
#fact[Bézout, CRT][$a x + b y = c$ is solvable iff
  $gcd(a, b) | c$, and $a x equiv c thick (mod m)$ has
  $gcd(a, m)$ solutions modulo `m` if $gcd(a, m) | c$, none
  otherwise. The system $x equiv r_i thick (mod m_i)$ is
  solvable iff $r_i equiv r_j thick (mod gcd(m_i, m_j))$ for all
  pairs, and then unique modulo the lcm (`numtheory/euclid`). A
  problem modulo $n = product p_i^(e_i)$ splits into independent
  problems modulo each $p_i^(e_i)$.]

=== Orders and primitive roots
#fact[Order][For $gcd(a, m) = 1$, $op("ord")_m (a)$, the least
  $e >= 1$ with $a^e equiv 1 thick (mod m)$, divides every such `e`,
  hence $phi(m)$; $op("ord")(a^k) = op("ord")(a) \/
  gcd(k, op("ord")(a))$.]
#fact[Primitive roots][They exist only modulo 1, 2, 4, $p^k$ and
  $2 p^k$ (odd `p`), $phi(phi(m))$ of them; modulo $2^k$, $k >= 3$,
  the units are $plus.minus 5^j$. Modulo `p` with
  $g = gcd(k, p - 1)$: $x^k equiv a$, $p divides.not a$, has `g`
  solutions if $a^((p - 1)\/g) equiv 1$ and none otherwise; a
  primitive `k`-th root of unity exists iff $k | p - 1$.]

=== Quadratic residues
#fact[Legendre symbol][For odd `p`,
  $(a\/p) = a^((p - 1)\/2) mod p in {0, 1, -1}$ is 1 iff `a` is a
  nonzero square#full[ (root: `numtheory/mod-sqrt`)]. It is
  multiplicative in `a`; $(-1\/p) = 1$ iff $p equiv 1 thick (mod 4)$ and
  $(2\/p) = 1$ iff $p equiv plus.minus 1 thick (mod 8)$.]
#fact[Reciprocity][For odd primes $p != q$, $(p\/q) (q\/p) = -1$ iff
  `p` and `q` are both $3 mod 4$.]
#fact[Sums of squares][$n = x^2 + y^2$ iff every prime
  $equiv 3 thick (mod 4)$ divides `n` an even number of times; there
  are $4 (d_1 (n) - d_3 (n))$ ordered signed pairs, $d_i$ counting
  the divisors $equiv i thick (mod 4)$. Three squares suffice iff
  $n != 4^a (8 b + 7)$; four always do.]
#full[
#fact[Jacobi symbol][For odd $n = product p_i^(e_i)$,
  $(a\/n) = product (a\/p_i)^(e_i)$. It obeys the same rules and is
  computed like a gcd, without factoring `n`. $(a\/n) = -1$ proves
  that `a` is not a square modulo `n`; $+1$ proves nothing when `n`
  is composite.]
#lang-code(
  cpp: "int jacobi(ll a, ll n) { // n odd positive, a >= 0
  int s = 1;
  for (a %= n; a; a %= n) {
    for (; a % 2 == 0; a /= 2)
      if (n % 8 == 3 || n % 8 == 5) s = -s;
    swap(a, n);
    if (a % 4 == 3 && n % 4 == 3) s = -s;
  }
  return n == 1 ? s : 0;
}",
  py: "def jacobi(a, n):  # n odd positive
  a, s = a % n, 1
  while a:
    while a % 2 == 0:
      a //= 2
      if n % 8 in (3, 5):
        s = -s
    a, n = n, a
    if a % 4 == 3 and n % 4 == 3:
      s = -s
    a %= n
  return s if n == 1 else 0",
)
]

=== Factorials and binomials
Legendre's $v_p (n!)$, Kummer and Lucas: see Counting formulas.
#fact[Lucas period][For $0 <= a < p$, $binom(m, a) equiv
  binom(m mod p, a) thick (mod p)$, so $binom(m, a) mod p$, and
  any product of such binomials, has period `p` in `m`.]
#fact[Modulo $p^k$][Only multiples of `p` lack inverses: write
  each factor as $p^v u$ with $p divides.not u$, add and subtract
  the exponents `v` as integers and multiply the units
  $u mod p^k$ (to divide, multiply by $u^(-1)$). The value is
  $p^v u mod p^k$, zero once $v >= k$.]
#fact[Lifting the exponent][Odd `p`, $p | a - b$, $p divides.not a b$:
  $v_p (a^n - b^n) = v_p (a - b) + v_p (n)$; for odd `n` and
  $p | a + b$ also $v_p (a^n + b^n) = v_p (a + b) + v_p (n)$. For
  odd `a`, `b` and even `n`: $v_2 (a^n - b^n) = v_2 (a - b) +
  v_2 (a + b) + v_2 (n) - 1$; for odd `n` it is $v_2 (a - b)$.]

=== Multiplicative functions
For $n = product p_i^(e_i)$: $phi(n) = n product_i (1 - 1\/p_i)$,
$d(n) = product_i (e_i + 1)$,
$sigma(n) = product_i (p_i^(e_i + 1) - 1)\/(p_i - 1)$, and
$mu(n) = (-1)^k$ if `n` is a product of `k` distinct primes, else 0.
`numtheory/sieve` gives $phi$ and $mu$ up to $10^7$; its Tweak
adds $d$ and $sigma$.#py[ In PyPy keep only the tables you need:
all of them take about 340 MB at $10^7$.]
#fact[Totients][$phi(m n) = phi(m) phi(n) dot g\/phi(g)$ with
  $g = gcd(m, n)$; e.g. $phi(3k) = phi(4k)$ when $gcd(k, 6) = 1$.
  For $g | n$, $phi(n\/g)$ of the $x in [0, n)$ have
  $gcd(x, n) = g$.#full[ Jordan's $J_k (n) = n^k product_(p | n)
  (1 - p^(-k))$ counts the $k$-tuples modulo `n` whose gcd with
  `n` is 1; $sum_(d | n) J_k (d) = n^k$, and for $d | n$ exactly
  $J_k (d)$ elements of $(ZZ_n)^k$ have additive order `d`
  (Burnside over the translations of a $k$-dimensional torus).]]
#fact[Dirichlet][$(f * g)(n) = sum_(d | n) f(d) g(n\/d)$ keeps
  multiplicativity. With $1(n) = 1$, $op("id")(n) = n$ and
  $epsilon(n) = [n = 1]$: $mu * 1 = epsilon$, $phi * 1 = op("id")$,
  $phi = mu * op("id")$, $1 * 1 = d$, $op("id") * 1 = sigma$.
  Möbius inversion: $g = f * 1 <=> f = g * mu$; over multiples,
  $g(d) = sum_(d | m) f(m) <=> f(d) = sum_(d | m) mu(m\/d) g(m)$.]
#fact[gcd sums][Write $[gcd(i, j) = 1] = sum_(d | i, d | j) mu(d)$
  or $gcd(i, j) = sum_(d | i, d | j) phi(d)$ and sum over `d` first:
  the pairs $i <= n$, $j <= m$ with $gcd(i, j) = 1$ number
  $sum_d mu(d) floor(n\/d) floor(m\/d)$; for $gcd = g$ divide `n`
  and `m` by `g`. Also $sum_(i <= n) gcd(i, n) = sum_(d | n) d
  phi(n\/d)$, and the $i <= n$ coprime to $n >= 2$ sum to
  $n phi(n)\/2$. As $gcd(x, s - x) = gcd(x, s)$, the coprime
  pairs $x + y = s$ with $x, y >= 1$ number $phi(s)$, $s >= 2$.]
#fact[Floor blocks][$floor(n\/i)$ takes at most $2 sqrt(n)$ values
  and is constant for $i in [l, floor(n\/floor(n\/l))]$, so
  $sum_(i <= n) f(i) g(floor(n\/i))$ costs $O(sqrt(n))$ given prefix
  sums of `f`. With two quotients end each block at the smaller of
  the two ends. Since $floor(floor(n\/a)\/b) = floor(n\/(a b))$,
  recursions on $floor(n\/i)$ stay inside these values.]
#lang-code(
  cpp: "for (ll l = 1, r; l <= n; l = r + 1) {
  r = n / (n / l); // n / i is the same for i in [l, r]
}",
  py: "l = 1
while l <= n:
  r = n // (n // l)  # n // i is the same for i in [l, r]
  l = r + 1",
)
#fact[Divisors][$sum_(i <= n) d(i) = sum_(k <= n) floor(n\/k)$ and
  $sum_(i <= n) sigma(i) = sum_(k <= n) k floor(n\/k)$. The largest
  $d(n)$: the table in Numbers, estimates and numerics.]

=== Primes
#fact[Gaps][There is a prime in $(n, 2n]$ (Bertrand). Largest gaps
  and $pi(n)$: Numbers, estimates and numerics.]
#fact[Near powers of ten][The largest prime below $10^k$ is
  $10^k - delta$, for $k = 1, dots, 18$: $delta =$ 3, 3, 3, 27, 9, 17,
  9, 11, 63, 33, 23, 11, 29, 27, 11, 63, 3, 11. The smallest above is
  $10^k + delta$: 1, 1, 9, 7, 3, 3, 19, 7, 7, 19, 3, 39, 37, 31, 37,
  61, 3, 3.]
#fact[Goldbach][Every even $4 <= n <= 4 dot 10^18$ is a sum of two
  primes (verified; the smaller one is tiny, so try $p = 2, 3, 5,
  dots$), and every odd $n >= 7$ is a sum of three primes (proved).]

=== Continued fractions
#fact[Convergents][$x = [a_0; a_1, a_2, dots]$: $a_0 = floor(x)$,
  then continue with $1\/(x - a_0)$; for $x = p\/q$ the $a_i$ are
  Euclid's quotients. With $h_(-2) = 0$, $h_(-1) = 1$,
  $k_(-2) = 1$, $k_(-1) = 0$, $h_i = a_i h_(i - 1) + h_(i - 2)$ and
  $k_i = a_i k_(i - 1) + k_(i - 2)$, $h_i\/k_i = [a_0; dots, a_i]$
  in lowest terms, $h_i k_(i - 1) - h_(i - 1) k_i = (-1)^(i - 1)$
  and $|x - h_i\/k_i| <= 1\/(k_i k_(i + 1))$. For $i >= 1$ a
  fraction closer to `x` than $h_i\/k_i$ has a denominator
  above $k_i$: truncating gives the best approximations.]
#fact[Square roots][For a non-square $D$, from $m = 0$, $d = 1$,
  $a = a_0 = floor(sqrt(D))$ repeat $m <- d a - m$,
  $d <- (D - m^2)\/d$, $a <- floor((a_0 + m)\/d)$: exact integers,
  periodic after $a_0$.]

=== Miscellany
#fact[Frobenius][For coprime $a, b >= 1$ the largest integer that
  is not $a x + b y$ with $x, y >= 0$ is $a b - a - b$, and
  $(a - 1)(b - 1)\/2$ integers are not. With more coins, Dijkstra on
  residues modulo the smallest coin `c` gives the least
  representable $t_r$ in each class $r$; `n` is representable iff
  $n >= t_(n mod c)$.]
#fact[Pythagorean triples][The primitive ones are
  $(m^2 - n^2, 2 m n, m^2 + n^2)$ with $m > n >= 1$,
  $gcd(m, n) = 1$ and $m - n$ odd, each exactly once; the others are
  their multiples. About $N\/(2 pi)$ primitive triples have
  hypotenuse $<= N$.]
#fact[Zeckendorf][Every $n >= 1$ is uniquely a sum of
  non-consecutive Fibonacci numbers $1, 2, 3, 5, 8, dots$; take the
  largest one $<= n$ greedily.]
#fact[Remainders][$x mod y < x\/2$ for $1 <= y <= x$, so a value
  changes at most $log_2 x + 1$ times under $x <- x mod a_i$, and
  a remainder $c < x$ needs $x > 2c$.]
#fact[Differences][Adding $|x - y|$ for members `x`, `y` until
  nothing changes turns a set of positive integers into exactly
  the multiples of its gcd up to its maximum.]
#full[
#fact[Digit sums][A positive multiple of $b^s - 1$ has base-`b`
  digit sum at least $s (b - 1)$; in binary, at least `s` ones.
  Proof: while the top digit sits at $b^e$, $e >= s$, subtract
  $b^e$ and add $b^(e - s)$: the residue stays, the number
  shrinks, the digit sum never grows, and $b^s - 1$ remains.]
#fact[Small solutions][Among the optimal solutions of: minimise
  $sum b_i x_i$ ($b_i > 0$) over integers $x_i >= 0$ with
  $sum a_i x_i equiv m thick (mod n)$, some has
  $product_i (x_i + 1) <= n$. Otherwise two vectors
  $0 <= V != W <= x$ have equal $sum a_i V_i mod n$; $x - V + W$
  and $x - W + V$ are solutions whose costs average to that of
  `x`, and one of them is lexicographically smaller. With four
  variables, enumerate the small ones and meet in the middle.]
#fact[Balanced words][A cyclic 0/1 word of length `n` with `q`
  ones in which any two windows of equal length differ by at most
  one in their number of ones has
  $w_i = floor((q (i + 1) + theta)\/n) - floor((q i + theta)\/n)$
  for some $0 <= theta < n$: the first `i` letters hold
  $floor((q i + theta)\/n)$ ones.]
#fact[Pell][$x^2 - D y^2 = 1$ for a non-square $D$: the smallest
  solution is the first convergent $h_i\/k_i$ of $sqrt(D)$ with
  $h_i^2 - D k_i^2 = 1$; all solutions are $x_j + y_j sqrt(D) =
  (x_1 + y_1 sqrt(D))^j$. Values explode: $D = 61$ needs
  $x_1 = 1766319049$.]
#fact[Continuants][The numerator $K(a_1, dots, a_n)$ of
  $[a_1; a_2, dots, a_n]$ satisfies $K(a_1, dots, a_n) =
  K(a_n, dots, a_1)$. A rational has exactly two expansions,
  $[dots, a_n]$ with $a_n >= 2$ (or $n = 0$) and
  $[dots, a_n - 1, 1]$. So sequences $(a_1, dots, a_n)$ of
  positive integers, $n >= 0$, match one-to-one the rationals
  $y\/x = [a_n; dots, a_1, 1] >= 1$ (empty: $1\/1$), i.e. the
  coprime pairs $x <= y$.]
#fact[Carmichael][$lambda(m)$, the least `e` with
  $a^e equiv 1$ for all `a` coprime to `m`: $lambda(p^k) = phi(p^k)$
  for odd `p`, $lambda(2) = 1$, $lambda(4) = 2$,
  $lambda(2^k) = 2^(k - 2)$ for $k >= 3$, and the lcm over the prime
  powers of `m`.]
#fact[Floor sums of degree 2][For $i = 0..n$ and
  $q_i = floor((a i + b)\/c)$ ($a, b >= 0$, $c >= 1$) let
  $f = sum q_i$, $g = sum i q_i$, $h = sum q_i^2$, and
  $S_1 = n (n + 1)\/2$, $S_2 = n (n + 1)(2n + 1)\/6$. If $a >= c$
  or $b >= c$, with $a = alpha c + a'$, $b = beta c + b'$ and
  $f', g', h'$ of $(a', b', c, n)$: $f = f' + alpha S_1 +
  beta (n + 1)$, $g = g' + alpha S_2 + beta S_1$, $h = h' +
  2 beta f' + 2 alpha g' + alpha^2 S_2 + beta^2 (n + 1) +
  2 alpha beta S_1$. Otherwise let $M = floor((a n + b)\/c)$; if
  $M = 0$ all are 0, else with $f', g', h'$ of
  $(c, c - b - 1, a, M - 1)$: $f = n M - f'$,
  $g = (M n (n + 1) - h' - f')\/2$, $h = n M (M + 1) - 2 g' -
  2 f' - f$. Depth $O(log)$ like Euclid; use `__int128`.]

=== Recognition
- Answer modulo a composite: solve modulo each prime power, then
  CRT. A huge exponent or a power tower: reduce it modulo $phi(m)$,
  keeping the $+ phi(m)$ when the true exponent is large. Counting
  modulo $2^k$ or 4: (valuation, unit) pairs, as in Modulo $p^k$.
- "gcd", "coprime", "lowest terms" inside a sum: Möbius or
  $phi$ substitution, then floor blocks.
- $a^n - 1$ divisible by a power of `p`, or the exact power of `p`
  in a product: LTE, Legendre, Kummer.
- "Write `n` as a sum of primes or squares": Goldbach, two, three
  and four squares. "Amounts payable with coins": Frobenius.
- $a^x equiv b$: `numtheory/discrete-log`; $x^2 equiv a$:
  `numtheory/mod-sqrt`; $x^k equiv a$ modulo `p`: take logs to a
  primitive root (`numtheory/primitive-root`). Counting or summing
  primes up to $10^11$: `numtheory/prime-count`.
- A process that keeps subtracting the smaller of two quantities,
  or shifts a group of values by a fixed amount per round, is
  subtractive Euclid: replace each run of equal steps by one
  division, $O(log)$ phases.
- $f_(i + 1) = a_i f_i + f_(i - 1)$, a best fraction, a network
  of unit resistors in series and parallel (its resistance is
  $[a_0; a_1, dots]$ with $sum a_i$ resistors): continued
  fractions.
- Sums of floors or of `mod` over an arithmetic progression:
  `numtheory/floor-sum`; with squares or an extra factor `i`, the
  degree-2 floor sums above.
]
