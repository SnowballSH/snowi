#import "/typst/notes.typ": *

=== Sums and estimates
#full[
- Positive integers with sum $n$ take fewer than $sqrt(2n)$
  distinct values ($k$ distinct ones sum to at least
  $k (k + 1)\/2$): group equal values, then loop over pairs.
]
- A sorted multiset contains a non-degenerate triangle iff three
  consecutive entries do ($a_i + a_(i + 1) > a_(i + 2)$). So
  triangle-free values in $[1, V]$ grow like Fibonacci: at most
  $L$ of them with $F_L <= V$ ($F_1 = F_2 = 1$), i.e. 27 for
  $V = 2 dot 10^5$, 44 for $10^9$, 87 for $10^18$.
- Stirling: $ln n! = n ln n - n + 1/2 ln(2 pi n) + O(1\/n)$;
  $binom(2n, n) approx 4^n\/sqrt(pi n)$.

=== Sizes of numbers
#table(
  columns: 4,
  [$n <=$], [max divisors], [at], [$pi(n)$ (primes)],
  [$10^3$], [32], [840], [168],
  [$10^4$], [64], [7560], [1229],
  [$10^5$], [128], [83160], [9592],
  [$10^6$], [240], [720720], [78498],
  [$10^7$], [448], [8648640], [664579],
  [$10^8$], [768], [73513440], [5761455],
  [$10^9$], [1344], [735134400], [50847534],
  [$10^12$], [6720], [963761198400], [$3.8 dot 10^10$],
  [$10^18$], [103680], [897612484786617600], [$2.5 dot 10^16$],
)
- Distinct prime factors: $omega(n) <= 7$ for $n <= 10^6$
  ($2 dot 3 dots 17 = 510510$), $<= 9$ for $n <= 10^9$
  ($2 dot 3 dots 23 = 223092870$) and $<= 15$ for $n <= 10^18$, so
  inclusion-exclusion over the primes of $n$ has at most $2^15$
  terms.
- Largest prime gap: 114 below $10^6$, 282 below $10^9$, 1442
  below $10^18$ (average $ln n$).
- Primes: $10^9 + 7$, $10^9 + 9$, $998244353 = 119 dot 2^23 + 1$,
  $10^18 + 3$, $10^18 + 9$, $2^61 - 1$.
- Largest that fits in #cpp[`ll`]#py[a machine int] ($2^63 approx
  9.2 dot 10^18$): $20!$, $binom(66, 33)$, Catalan $C_35$, Bell
  $B_25$, partitions $p(405)$. #cpp[`double`]#py[`float`] holds
  integers exactly up to $2^53 approx 9.0 dot 10^15$.
#table(
  columns: 5,
  [$n$], [$n!$], [$C_n$], [$B_n$], [$p(n)$],
  [10], [3628800], [16796], [115975], [42],
  [15], [$1.3 dot 10^12$], [9694845], [$1.4 dot 10^9$], [176],
  [20], [$2.4 dot 10^18$], [$6.6 dot 10^9$], [$5.2 dot 10^13$], [627],
  [25], [$1.6 dot 10^25$], [$4.9 dot 10^12$], [$4.6 dot 10^18$], [1958],
  [30], [$2.7 dot 10^32$], [$3.8 dot 10^15$], [$8.5 dot 10^23$], [5604],
  [50], [$3.0 dot 10^64$], [$2.0 dot 10^27$], [$1.9 dot 10^47$], [204226],
  [100], [$9.3 dot 10^157$], [$9.0 dot 10^56$], [$4.8 dot 10^115$],
  [$1.9 dot 10^8$],
)
#full[
- Numbers $2^a 3^b$: 306 up to $10^9$, 1178 up to $10^18$; few
  enough to enumerate.
- Modulo a Mersenne prime $p = 2^k - 1$ (e.g. $2^31 - 1$,
  $2^61 - 1$) and $x <= (p - 1)^2$: $x mod p$ is
  `(x >> k) + (x & p)`, minus `p` once if that is $>= p$.
]

=== Numerics
#fact[Stable quadratic][$q = -(b + "sgn"(b) sqrt(b^2 - 4 a c))\/2$
  with $"sgn"(0) = 1$, roots $q\/a$ and $c\/q$. The textbook formula
  loses digits of the small root when $b^2 >> |4 a c|$.]
#fact[Simpson][$integral_a^b f approx (b - a)\/6 dot (f(a) +
  4f(m) + f(b))$, $m = (a + b)\/2$, exact for cubics; error
  $O(h^4)$ with $n$ pieces of width $h$. Adaptive: split $[a, b]$
  while $|S_"left" + S_"right" - S| > 15 epsilon$ (halve
  $epsilon$ per level), then return $S_"left" + S_"right" +
  (S_"left" + S_"right" - S)\/15$.#full[ Volume of revolution:
  $pi integral f^2$.]]
#full[
#fact[Quadrature][Over a full period of a smooth periodic function
  the trapezoid rule $h sum_i f(i h)$ converges exponentially
  (for $e^(cos t)$ on $[0, 2 pi]$: error $10^(-6)$ with 8 points,
  where Simpson still has $10^(-2)$); double the points until two
  results agree. Split at kinks, singular points and sharp peaks
  before applying any rule.]
]
#fact[Precision][Clamp `acos`/`asin` arguments into $[-1, 1]$.
  `expm1(x)` and `log1p(x)` keep the digits of $e^x - 1$ and
  $ln(1 + x)$ for tiny `x`; a difference of two large close values
  (say $E[X^2] - E[X]^2$) loses them, so compute it exactly.]
#fact[KKT][Minimise $sum_i f_i (x_i)$, $f_i$ convex, with
  $sum_i x_i = C$, $l_i <= x_i <= u_i$: at the optimum some
  $lambda$ has $f'_i (x_i) = lambda$ where $l_i < x_i < u_i$,
  $f'_i >= lambda$ at $l_i$ and $<= lambda$ at $u_i$. Bisect
  $lambda$, with $x_i$ the root of $f'_i = lambda$ clamped to
  $[l_i, u_i]$, until $sum_i x_i = C$.#full[ In general, to minimise a
  convex differentiable `f` under convex differentiable
  constraints $g_j (x) <= 0$ (one point satisfies all strictly),
  a feasible `x` is optimal iff
  $nabla f(x) + sum_j mu_j nabla g_j (x) = 0$ for some
  $mu_j >= 0$ with $mu_j g_j (x) = 0$; usually one multiplier is
  found by bisection. Projection of `x` onto the box
  $|y_i| <= l_i$ cut by the ball $|y|^2 <= a$: if clamping `x`
  to the box lands in the ball, that is it; otherwise
  $y_i = op("clamp")(s x_i, -l_i, l_i)$ with $s in (0, 1]$ chosen
  so that $|y|^2 = a$.]]
#full[
#fact[Saturating sums][The sign of $sum_(i < n) 2^i a_i$ with all
  $|a_i| <= M$, without big integers: Horner from the top,
  $v <- 2 v + a_i$, clamping $v$ to $[-(M + 1), M + 1]$. Once
  $|v| = M + 1$ the sign never changes. A block of terms maps the
  clamped $v$ to the clamped $v$: store that map as its table of
  $2M + 3$ values and compose tables, so a segment tree can keep
  them.]
]

=== Fibonacci
#full[
- $F_(n + m) = F_n F_(m + 1) + F_(n - 1) F_m$; doubling:
  $F_(2n) = F_n (2 F_(n + 1) - F_n)$,
  $F_(2n + 1) = F_n^2 + F_(n + 1)^2$, so $F_n mod m$ in $O(log n)$.
]
- $gcd(F_n, F_m) = F_(gcd(n, m))$; $F_n | F_m$ iff $n | m$ (for
  $n > 2$).
- Pisano period $pi(m)$ (period of $F_n mod m$) is at most $6m$;
  $pi(op("lcm")(a, b)) = op("lcm")(pi(a), pi(b))$ and
  $pi(p^k) | p^(k - 1) pi(p)$. For a prime
  $p equiv plus.minus 1 (mod 5)$, $pi(p) | p - 1$; for
  $p equiv plus.minus 2 (mod 5)$, $pi(p) | 2(p + 1)$; $pi(5) = 20$.
#full[
- Zeckendorf: Number theory facts.
]
