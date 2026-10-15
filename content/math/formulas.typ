#import "/typst/notes.typ": *

=== Algebra
#fact[Quadratic][$a x^2 + b x + c = 0$, $D = b^2 - 4 a c$:
  $x = (-b plus.minus sqrt(D))\/(2a)$. Without cancellation:
  $q = -(b + "sgn"(b) sqrt(D))\/2$ with $"sgn"(0) = 1$
  (`copysign(sqrt(D), b)`), $x_1 = q\/a$, $x_2 = c\/q$; $q = 0$
  only if $b = c = 0$, a double root 0.]
#fact[Vieta][For $a_n product_i (x - x_i)$: $sum_i x_i = -a_(n - 1)\/a_n$,
  $sum_(i < j) x_i x_j = a_(n - 2)\/a_n$, \
  $product_i x_i = (-1)^n a_0\/a_n$.]
#fact[Cubic][$a x^3 + b x^2 + c x + d$; $x = t - b\/(3a)$ gives
  $t^3 + p t + q = 0$ with \
  $p = (3 a c - b^2)\/(3a^2)$,
  $q = (2b^3 - 9 a b c + 27 a^2 d)\/(27 a^3)$. \
  If $p < 0$ and $4p^3 + 27q^2 <= 0$ all roots are real (repeated
  at equality; clamp the arccos argument to $[-1, 1]$),
  $k = 0, 1, 2$: \
  $t_k = 2 sqrt(-p\/3) cos(1/3 arccos((3q)/(2p) sqrt(-3\/p)) -
  (2 pi k)/3)$. \
  Otherwise one real root: $t = root(3, -q\/2 + sqrt(D)) +
  root(3, -q\/2 - sqrt(D))$, $D = q^2\/4 + p^3\/27$. \
  In practice: bisect for one root, divide it out, solve the
  quadratic.]

=== Sums
- $sum_(i = 1)^n i = n(n + 1)\/2$, $sum i^2 = n(n + 1)(2n + 1)\/6$, \
  $sum i^3 = (n(n + 1)\/2)^2$,
  $sum i^4 = n(n + 1)(2n + 1)(3n^2 + 3n - 1)\/30$; in general a
  polynomial of degree $k + 1$ (interpolate, `poly/lagrange`).
- $sum_(i = a)^b c^i = (c^(b + 1) - c^a)\/(c - 1)$,
  $sum_(i >= 0) i r^i = r\/(1 - r)^2$ for $|r| < 1$.
- $H_n = sum_(i = 1)^n 1/i approx ln n + 0.5772 + 1\/(2n)$,
  $sum 1\/i^2 = pi^2\/6$, \
  $sum_(i <= n) floor(n\/i) approx n ln n$,
  $sum_(p <= n) 1/p approx ln ln n + 0.2615$. \
  $floor(n\/i)$ takes at most $2 sqrt(n)$ distinct values.

=== Series
- $e^x = sum_n x^n\/n!$,
  $(1 + x)^a = sum_k binom(a, k) x^k$ for $|x| < 1$.
- $ln(1 + x) = x - x^2\/2 + x^3\/3 - dots$ for $-1 < x <= 1$.
- $sqrt(1 + x) = 1 + x\/2 - x^2\/8 + x^3\/16 - 5x^4\/128 + dots$
- $sin x = x - x^3\/3! + x^5\/5! - dots$,
  $cos x = 1 - x^2\/2! + x^4\/4! - dots$
- $arctan x = x - x^3\/3 + x^5\/5 - dots$,
  Stirling: $n! approx sqrt(2 pi n) (n\/e)^n$.

=== Trigonometry
- $sin(v plus.minus w) = sin v cos w plus.minus cos v sin w$ \
  $cos(v plus.minus w) = cos v cos w minus.plus sin v sin w$ \
  $tan(v + w) = (tan v + tan w)\/(1 - tan v tan w)$
- $sin 2v = 2 sin v cos v$, $cos 2v = 2cos^2 v - 1 = 1 - 2sin^2 v$
- $sin v + sin w = 2 sin((v + w)\/2) cos((v - w)\/2)$ \
  $cos v + cos w = 2 cos((v + w)\/2) cos((v - w)\/2)$
- $a cos x + b sin x = r cos(x - phi)$,
  $a sin x + b cos x = r sin(x + phi)$ \
  with $r = sqrt(a^2 + b^2)$, $phi = "atan2"(b, a)$;
  $e^(i x) = cos x + i sin x$.

=== Triangles
Sides $a, b, c$ opposite angles $A, B, C$; $s = (a + b + c)\/2$.
- Area: $sqrt(s(s - a)(s - b)(s - c))$, also $1/2 a b sin C$, $r s$
  and $a b c\/(4R)$. Stable Heron with $a >= b >= c$: \
  $1/4 sqrt((a + (b + c))(c - (a - b))(c + (a - b))(a + (b - c)))$.
- Sines: $a\/sin A = b\/sin B = c\/sin C = 2R$. \
  Cosines: $c^2 = a^2 + b^2 - 2 a b cos C$. \
  Tangents: $(a + b)\/(a - b) = tan((A + B)\/2)\/tan((A - B)\/2)$.
- Inradius $r = "area"\/s$, circumradius $R = a b c\/(4 dot "area")$.
- Median to $a$: $1/2 sqrt(2b^2 + 2c^2 - a^2)$. Bisector of $A$:
  $sqrt(b c (1 - (a\/(b + c))^2))$; it meets $a$ at $D$ with
  $B D : D C = c : b$.

=== Quadrilaterals
Sides $a, b, c, d$, diagonals $e, f$ at angle $theta$,
$F = b^2 + d^2 - a^2 - c^2$: \
$4 dot "area" = 2 e f sin theta = F tan theta =
sqrt(4 e^2 f^2 - F^2)$. \
Cyclic (opposite angles sum to $pi$): Ptolemy $e f = a c + b d$,
Brahmagupta area $sqrt((s - a)(s - b)(s - c)(s - d))$.

=== Spherical coordinates
$x = r sin theta cos phi$, $y = r sin theta sin phi$,
$z = r cos theta$; \
$r = sqrt(x^2 + y^2 + z^2)$, $theta = arccos(z\/r)$,
$phi = "atan2"(y, x)$. \
Sphere: volume $4/3 pi r^3$, area $4 pi r^2$; cap of height $h$:
volume $pi h^2 (3r - h)\/3$, area $2 pi r h$. Cone: $pi r^2 h\/3$.

=== Calculus
- $(tan x)' = 1 + tan^2 x$, $(arctan x)' = 1\/(1 + x^2)$, \
  $(arcsin x)' = 1\/sqrt(1 - x^2)$, $(a^x)' = a^x ln a$.
- $integral tan a x dif x = -ln|cos a x|\/a$,
  $integral x e^(a x) dif x = e^(a x) (a x - 1)\/a^2$, \
  $integral x sin a x dif x = (sin a x - a x cos a x)\/a^2$,
  $integral e^(-x^2) dif x = sqrt(pi)\/2 "erf"(x)$. \
  By parts: $integral_a^b f g' = [f g]_a^b - integral_a^b f' g$.
- Simpson: $integral_a^b f approx (b - a)\/6 dot (f(a) +
  4f((a + b)\/2) + f(b))$, exact for cubics; on $n$ pieces the
  error is $O(h^4)$.
- Revolving $y = f(x)$: volume $pi integral f^2 dif x$; arc length
  $integral sqrt(1 + f'^2) dif x$.

=== Numbers
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
- Largest prime gap: 114 below $10^6$, 282 below $10^9$, under
  1500 below $10^18$. Primes: $10^9 + 7$, $10^9 + 9$, $998244353$,
  $10^18 + 3$, $10^18 + 9$, $2^61 - 1$.
- Last value that fits `ll` ($2^63 approx 9.2 dot 10^18$): $20!$,
  $binom(66, 33)$, Catalan $C_35$, Bell $B_25$, partitions
  $p(405)$. `double` is exact for integers up to $2^53 approx
  9.0 dot 10^15$.
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
=== More identities
- $a^n - b^n = (a - b) sum_(i < n) a^(n - 1 - i) b^i$;
  $a^3 + b^3 = (a + b)(a^2 - a b + b^2)$.
- AM-GM: $(x_1 + dots + x_n)\/n >= root(n, x_1 dots x_n)$ for
  $x_i >= 0$. Cauchy-Schwarz: $(sum a_i b_i)^2 <= sum a_i^2 sum b_i^2$.
- $sin v sin w = (cos(v - w) - cos(v + w))\/2$,
  $cos v cos w = (cos(v - w) + cos(v + w))\/2$,
  $sin v cos w = (sin(v + w) + sin(v - w))\/2$,
  $sin 3v = 3 sin v - 4 sin^3 v$, $cos 3v = 4 cos^3 v - 3 cos v$,
  $tan(v\/2) = sin v\/(1 + cos v)$.
- Any quadrilateral (Bretschneider): area$""^2 = (s - a)(s - b)
  (s - c)(s - d) - a b c d cos^2((alpha + gamma)\/2)$ with $alpha,
  gamma$ opposite angles.
- Great-circle angle between latitudes $phi_i$, longitudes
  $lambda_i$ (haversine): $2 arcsin sqrt(sin^2(Delta phi\/2) +
  cos phi_1 cos phi_2 sin^2(Delta lambda\/2))$.
- Newton: $x <- x - f(x)\/f'(x)$ doubles the correct digits near a
  simple root; bisection is slower but never fails.
- Squares: $sum_(i = 1)^n (2i - 1) = n^2$; triangular numbers
  $n(n + 1)\/2$; Fibonacci: $F_(n + m) = F_n F_(m + 1) +
  F_(n - 1) F_m$, $gcd(F_n, F_m) = F_(gcd(n, m))$.
- $n$-th prime $approx n ln n$; $pi(n) approx n\/ln n$; the average
  gap near $n$ is $ln n$, so a random search finds a prime fast.
- Pitfalls: `acos`/`asin` of values slightly outside $[-1, 1]$
  give NaN, so clamp; `pow` with doubles is inexact for integer
  powers; compare doubles with an epsilon relative to the input
  scale.
]
