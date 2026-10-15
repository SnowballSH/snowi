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
  $gcd(a, b) | c$. The system $x equiv r_i thick (mod m_i)$ is
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
#icpc[Legendre's $v_p (n!)$, Kummer and Lucas: see Counting
formulas.]
#full[
#fact[Legendre][$v_p (n!) = sum_(i >= 1) floor(n\/p^i) =
  (n - s_p (n))\/(p - 1)$, $s_p$ the base-`p` digit sum. $n!$ ends
  in $v_5 (n!)$ zeros.]
#fact[Kummer, Lucas][$v_p binom(m + n, m)$ is the number of carries
  when adding `m` and `n` in base `p`;
  $binom(n, k) equiv product_i binom(n_i, k_i) thick (mod p)$ over
  the base-`p` digits.]
]
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
  $n phi(n)\/2$.]
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
  $sum_(i <= n) sigma(i) = sum_(k <= n) k floor(n\/k)$; visiting the
  multiples of every $k <= n$ costs $n ln n$. #icpc[Maximum $d(n)$:
  the Numbers table in Formulas.]#full[The most divisors: 240 for
  $n <= 10^6$, 1344 for $10^9$, 6720 for $10^12$, 103680 for
  $10^18$ (at 720720, 735134400, 963761198400,
  897612484786617600).]]

=== Primes
#fact[Gaps][There is a prime in $(n, 2n]$ (Bertrand). #icpc[Largest
  gaps and $pi(n)$: the Numbers table in Formulas.]#full[The largest
  gap between consecutive primes is 114 below $10^6$, 282 below
  $10^9$ and 1442 below $10^18$; the average gap near `n` is
  $ln n$.]]
#fact[Near powers of ten][The largest prime below $10^k$ is
  $10^k - delta$, for $k = 1, dots, 18$: $delta =$ 3, 3, 3, 27, 9, 17,
  9, 11, 63, 33, 23, 11, 29, 27, 11, 63, 3, 11. The smallest above is
  $10^k + delta$: 1, 1, 9, 7, 3, 3, 19, 7, 7, 19, 3, 39, 37, 31, 37,
  61, 3, 3.]
#fact[Goldbach][Every even $4 <= n <= 4 dot 10^18$ is a sum of two
  primes (verified; the smaller one is tiny, so try $p = 2, 3, 5,
  dots$), and every odd $n >= 7$ is a sum of three primes (proved).]

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
#full[
#fact[Pell][$x^2 - D y^2 = 1$ for a non-square $D$: expand
  $sqrt(D)$ as a continued fraction with $a_0 = floor(sqrt(D))$,
  $m <- d a - m$, $d <- (D - m^2)\/d$, $a <- floor((a_0 + m)\/d)$
  from $m = 0$, $d = 1$; the convergents $h_i = a h_(i - 1) +
  h_(i - 2)$, $k_i = a k_(i - 1) + k_(i - 2)$ (from $h_0 = a_0$,
  $k_0 = 1$, $h_(-1) = 1$, $k_(-1) = 0$) reach the smallest
  solution as the first $h^2 - D k^2 = 1$. All solutions:
  $x_j + y_j sqrt(D) = (x_1 + y_1 sqrt(D))^j$. Values explode:
  $D = 61$ needs $x_1 = 1766319049$.]
#fact[Carmichael][$lambda(m)$, the least `e` with
  $a^e equiv 1$ for all `a` coprime to `m`: $lambda(p^k) = phi(p^k)$
  for odd `p`, $lambda(2) = 1$, $lambda(4) = 2$,
  $lambda(2^k) = 2^(k - 2)$ for $k >= 3$, and the lcm over the prime
  powers of `m`.]

=== Recognition
- Answer modulo a composite: solve modulo each prime power, then
  CRT. A huge exponent or a power tower: reduce it modulo $phi(m)$,
  keeping the $+ phi(m)$ when the true exponent is large.
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
- Periodic sequences modulo `m`: Fibonacci modulo `p` repeats with
  a period dividing $p - 1$ when $p equiv plus.minus 1 thick
  (mod 5)$ and $2(p + 1)$ otherwise ($p != 5$).
]
