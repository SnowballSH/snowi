#import "/typst/notes.typ": *

=== Expectation
#fact[Indicators][Linearity needs no independence: write a count as
  a sum of indicators, $E[\#"events"] = sum_t P(t "happens")$. But
  $E[1\/X] != 1\/E[X]$ and $E[max] != max E$.]
#fact[Contribution][Sum of $f$ over all $N$ objects $= N dot E[f]$ for
  a uniform object, and back. Count per element ("$P$(cell $c$ is
  empty)"), per pair ("$P(u, v$ end up inverted$)$"), or per unit of
  length (each gap is a token counted when crossed).]
#fact[Tail sum][For integer $X >= 0$: $E[X] = sum_(k >= 1) P(X >= k)$
  and $E[X^2] = sum_(k >= 1) (2k - 1) P(X >= k)$; for real $X >= 0$:
  $E[X] = integral_0^infinity P(X > t) dif t$. Maxima and minima of
  independent variables: \
  $P(max <= k) = product_i P(Y_i <= k)$,
  $P(min >= k) = product_i P(Y_i >= k)$.]
#fact[Second moment][For indicators $X = sum_t I_t$:
  $E[X^2] = E[X] + sum_(s != t) P(s "and" t)$; variance adds over
  pairwise uncorrelated terms.]

=== Distributions
- Geometric (tries until the first success, probability $p$): mean
  $1\/p$, variance $(1 - p)\/p^2$; memoryless.
- Hypergeometric ($n$ drawn without replacement from $N$, $K$
  marked): mean $n K\/N$, variance
  $n K/N (N - K)/N (N - n)/(N - 1)$.
- Poisson $lambda$: $P(k) = e^(-lambda) lambda^k\/k!$, mean = variance
  $= lambda$; binomial with small $p$ is about Poisson $n p$.
- Exponential races: the minimum of independent exponentials has
  rate $sum lambda_i$, and the $i$-th finishes first with
  probability $lambda_i\/sum lambda_i$.
- Uniform integer on $m$ values: variance $(m^2 - 1)\/12$.

=== Markov chains
- Expected cost to absorption, backwards: $E_s = c_s + sum_t
  P_(s t) E_t$, $E = 0$ on absorbing states. Absorption
  probabilities: the same with $c = 0$ and boundary values 1 or 0.
- Self-loop $P_(s s) < 1$: $E_s = (c_s + sum_(t != s) P_(s t) E_t)\/
  (1 - P_(s s))$; success probability $p$ per try: $1\/p$ tries.
- Structure first: DAG order; tree walk, write $E_v = A_v E_"par" +
  B_v$ bottom-up; everything jumps back to the start, write
  $E_i = a_i E_0 + b_i$ and solve one equation; otherwise Gaussian
  elimination in $O(S^3)$ (`linalg/gauss-mod`).
- Distribution after $t$ steps: $pi_0 P^t$ by matrix power.
  Stationary: $pi P = pi$, $sum pi = 1$ (replace one equation);
  unique for an irreducible chain; expected return time to $s$ is
  $1\/pi_s$. Random walk on an undirected graph: $pi_v = deg v\/(2m)$.
- Tree walk from $v$ to its parent: $2 |"subtree"(v)| - 1$ steps on
  average.

=== Classics
- Coupon collector, $n$ equal coupons: $n H_n$; to see $k$ distinct:
  $sum_(i < k) n\/(n - i)$; distinct after $m$ draws:
  $n (1 - (1 - 1\/n)^m)$. Unequal $p_i$: min-max inclusion-exclusion
  $sum_(T != emptyset) (-1)^(|T| + 1) \/ sum_(i in T) p_i$.
- Birthday: $k$ draws from $n$ all differ with probability
  $product_(i < k) (1 - i\/n) approx e^(-k^2\/(2n))$; 50% at
  $k approx 1.18 sqrt(n)$; first repeat after about
  $sqrt(pi n\/2)$ draws. $m$ hashes modulo $M$ collide with
  probability about $m^2\/(2M)$.
- Gambler's ruin from $i$, absorbed at 0 or $N$, up with probability
  $p$, $r = (1 - p)\/p$: reaches $N$ with probability
  $(1 - r^i)\/(1 - r^N)$. If $p = 1\/2$: probability $i\/N$ and
  $i (N - i)$ steps on average.
- Symmetric $plus.minus 1$ walk: $P(S_(2n) = 0) = binom(2n, n)\/4^n
  approx 1\/sqrt(pi n)$, $E|S_n| approx sqrt(2n\/pi)$; it returns to 0
  surely, but the expected return time is infinite.
- $n$ uniform points in $[0, 1]$: the $k$-th smallest has mean
  $k\/(n + 1)$ ($E[max] = n\/(n + 1)$), all $n + 1$ gaps have mean
  $1\/(n + 1)$, and $P(max <= x) = x^n$.
- Random permutation: fixed points 1, cycles and records $H_n$,
  inversions $n(n - 1)\/4$; the cycle of a given element has a
  uniform length in $1..n$; two given elements share a cycle with
  probability $1\/2$.
- Waiting for a pattern $s$ over an alphabet of size $A$:
  $sum A^k$ over the $k$ with `s[0, k) == s[n - k, n)`, $k = n$
  included (HH: 6 flips, HT: 4).

=== Stopping times
#fact[Wald][$E[sum_(i = 1)^N X_i] = E[N] E[X]$ for i.i.d. $X_i$ and a
  stopping time $N$ with finite mean.]
#fact[Potential][Find $Phi$ of the state with $E[Phi_"next" | "now"] =
  Phi_"now" - 1$ while running; then the expected number of steps
  is $Phi_"start" - E[Phi_"end"]$ (optional stopping, with bounded
  increments and $E[T] < infinity$). Pattern waiting times and
  random acquisitions until one company is left (CF 1025G) are
  solved this way.]
#fact[Optimal stopping][$V(s) = max("stop"(s), -c + sum_t P_(s t)
  V(t))$ in floating point: the max needs real comparisons. Look
  for a threshold rule, or a one-step lookahead when stopping stays
  better once it is; then compute the value mod $p$ if asked.]

=== Modulo a prime
- A rational answer $p\/q$ is printed as $p q^(-1)$ mod $M$; $+, -,
  times, \/$ all commute with reduction. Divide once by the total
  (count favourable outcomes, then multiply by the inverse of the
  number of outcomes).
- Never compare residues or take max or min of them; never drop
  "tiny" probabilities. Decide the strategy first in doubles or by
  proof.
- Inverses of $1..n$ in $O(n)$: $"inv"(1) = 1$,
  $"inv"(i) = (M - floor(M\/i)) dot "inv"(M mod i) mod M$. To read
  a residue $r$ as a fraction, find a small $q$ with $r q mod M$
  small.

#full[
=== Inequalities
Markov: $P(X >= a) <= E[X]\/a$ for $X >= 0$. Chebyshev:
$P(|X - mu| >= k sigma) <= 1\/k^2$. Chernoff for a sum of
independent $[0, 1]$ variables with mean $mu$:
$P(X >= (1 + delta) mu) <= e^(-delta^2 mu\/3)$ ($0 < delta <= 1$),
$P(X <= (1 - delta) mu) <= e^(-delta^2 mu\/2)$. Hoeffding:
$P(X - E X >= t) <= e^(-2t^2\/n)$. Union bound:
$P(union.big A_i) <= sum P(A_i)$. Use them to size random trials:
a test that fails with probability $1\/2$ repeated 60 times fails
with probability $2^(-60)$.

=== Worked examples
- Maximum of two dice: $P(M <= v) = (v\/6)^2$, so
  $E[M] = sum_(v = 0)^5 (1 - v^2\/36) = 161\/36$.
- Walk on the path 0 - 1 - 2 absorbed at 0: $E_2 = 1 + E_1$,
  $E_1 = 1 + (E_0 + E_2)\/2$, so $E_1 = 3 = 2 dot 2 - 1$ (subtree
  of 1 has 2 vertices) and $E_2 = 4$.
- Expected maximum of $n$ dice with $k$ faces:
  $sum_(v = 0)^(k - 1) (1 - (v\/k)^n)$.
- Biased ruin duration ($p != 1\/2$, $q = 1 - p$):
  $i\/(q - p) - N\/(q - p) dot (1 - (q\/p)^i)\/(1 - (q\/p)^N)$.
- Absorbing chain with transient part $Q$ and $R$ to absorbing
  states: $(I - Q)^(-1)$ holds expected visits, its row sums the
  expected steps, $(I - Q)^(-1) R$ the absorption probabilities.

=== Recognition
- A uniformly random process often hides a uniform permutation:
  removing random elements one by one leaves a uniform random
  subset; in a shuffled deck every order of any $k$ given cards is
  equally likely.
- "Expected number of X": indicators, one probability per
  element or pair.
- "Expected maximum / time until all": tail sums with products,
  or min-max inclusion-exclusion.
- "Until the array is empty / sorted / absorbed": backward
  expectation DP; find the structure of the equations.
- "Output $p q^(-1)$ mod 998244353": no comparisons, count and
  divide once.
- "You may stop at any time", "play optimally": Bellman in
  doubles, then look for a threshold.
- Hidden information but a strategy that meets an offline lower
  bound on every outcome: the answer is $E["offline optimum"]$.
- Always check a formula on tiny cases with exact fractions
  (Python `fractions.Fraction`).

=== Pitfalls
- Doubles: compute tiny probabilities like $binom(n, k)\/binom(N, k)$
  with `lgamma` or running ratios. A #cpp[`double`]#py[`float`]
  keeps about 16 significant digits, so a relative error of
  $10^(-6)$ is easy, but an answer near $10^11$ with absolute error
  $10^(-6)$ needs #cpp[`long double`]#py[`Fraction` or a
  high-precision `decimal`].
- $sum_(k >= 0) P(X >= k) = E[X] + 1$: start the tail sum at 1.
- "Uniform random move" differs from "uniform random piece, then
  uniform direction"; read who is uniform over what.
- Draws without replacement are dependent: use hypergeometric
  counts or exchangeability.
- A state with self-loop probability 1, or an unreachable
  absorbing set, has infinite expectation: check reachability.
#cpp[
  - Subtract with `+ M` before `%`: C++ remainders can be
    negative.
]
]
