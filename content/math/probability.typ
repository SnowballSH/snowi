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
  pairwise uncorrelated terms. Total variance:
  $op("Var") X = E[op("Var")(X | Y)] + op("Var")(E[X | Y])$, so a
  mixture varies at least as much as its least varying part.
  Mixing averages the parts' points $(E[X], E[X^2])$ and the
  variance $E[X^2] - E[X]^2$ is concave in them: the largest one
  mixes two parts adjacent on the upper hull of these points (a
  parabola in the weight).]

=== Random subsets and orders
- Removing uniformly random elements one at a time: the $k$
  survivors are a uniform $k$-subset. If the stopping property is
  inherited by subsets (sorted, independent, pairwise far), the
  number $X$ left at the stop has $P(X >= k) = g_k \/ binom(n, k)$
  with $g_k$ the good $k$-subsets, so the expected number of
  removals is $n - sum_(k = 1)^n g_k \/ binom(n, k)$.
- In a uniform permutation of $n$ items the $K$-th of $N$ marked
  items has expected position $K (n + 1)\/(N + 1)$: the other
  $n - N$ items spread evenly over the $N + 1$ gaps.
- $k$ independent draws, outcome $s$ with probability $p_s$: the
  expected number of distinct outcomes is
  $sum_s (1 - (1 - p_s)^k) = sum_(i = 1)^k (-1)^(i + 1) binom(k, i)
  sum_s p_s^i$; the power sums factorise when $p_s$ is a product
  over independent choices.
- Hypergeometric ($n$ drawn without replacement from $N$, $K$
  marked): mean $n K\/N$, variance
  $n K/N (N - K)/N (N - n)/(N - 1)$.
- Exponential races: the minimum of independent exponentials has
  rate $sum lambda_i$, and the $i$-th finishes first with
  probability $lambda_i\/sum lambda_i$.

=== Markov chains
- Expected cost to absorption, backwards: $E_s = c_s + sum_t
  P_(s t) E_t$, $E = 0$ on absorbing states. Absorption
  probabilities: the same with $c = 0$ and boundary values 1 or 0.
- Self-loop $P_(s s) < 1$: $E_s = (c_s + sum_(t != s) P_(s t) E_t)\/
  (1 - P_(s s))$; success probability $p$ per try: $1\/p$ tries.
- Structure first: DAG order; tree walk, write $E_v = A_v E_"par" +
  B_v$ bottom-up; everything jumps back to the start, write
  $E_i = a_i E_0 + b_i$ and solve one equation; states on a line
  (or one level at a time) give tridiagonal systems, solved in
  $O(S)$ (Thomas algorithm, Linear algebra facts); otherwise
  Gaussian elimination in $O(S^3)$ (`linalg/gauss-mod`).
- Stationary: $pi P = pi$, $sum pi = 1$ (replace one equation);
  unique for an irreducible chain; expected return time to $s$ is
  $1\/pi_s$. Random walk on an undirected graph: $pi_v = deg v\/(2m)$.
  If the columns of $P$ also sum to 1 (doubly stochastic: every
  state entered by $d$ arcs and left along one of its $d$ arcs
  uniformly), a uniform distribution stays uniform.
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
  probability $1\/2$. Its longest increasing subsequence has length
  about $2 sqrt(n) - 1.77 n^(1\/6)$ (greedy monotone piles: about
  $2 sqrt(n)$).
- Waiting for a pattern $s$ over an alphabet of size $A$:
  $sum A^k$ over the $k$ with `s[0, k) == s[n - k, n)`, $k = n$
  included (HH: 6 flips, HT: 4).
- Random translation: $integral "area"(P inter (Q + t)) dif t =
  "area"(P) "area"(Q)$ (Fubini). For $t$ uniform over the
  translations with overlap (the polygon $P - Q$,
  `geometry/minkowski`) the expected overlap is
  $"area"(P) "area"(Q) \/ "area"(P - Q)$.

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
#full[
#fact[Index policy][Independent Markov chains, each with a target;
  a step advances one chosen chain at cost 1, and play stops when
  any chain reaches its target (golf with several balls). It is
  optimal to always advance the chain whose current state $s$ has
  the smallest index
  $g(s) = min_C E[T_C] \/ P("target before leaving" C)$
  over state sets $C in.rev s$, where $T_C$ counts the
  steps until the target or the first step out of $C$: the least
  expected time to the target from $s$ if you may restart at $s$
  at any moment (Gittins; Dumitriu, Tetali, Winkler). For a simple
  random walk on a tree with the target at the root the index never
  decreases away from the root, so the best $C$ are the connected
  sets $g <= G$. Grow them from the root: the next vertex is the
  frontier vertex $s$ (parent $p$ in $C$) with the least
  $(deg s + E(p)) \/ a(p)$, which is $g(s)$; $a(p)$ and $E(p)$ are
  the success probability and the expected steps from $p$ of the
  walk stopped when it leaves $C$ ($O(n)$ per step).]
]

=== Modulo a prime
- Count favourable outcomes, then multiply once by the inverse of
  the number of outcomes.
- Never compare residues or take max or min of them; never drop
  "tiny" probabilities. Decide the strategy first in doubles or by
  proof.
- Inverses of $1..n$ in $O(n)$: $"inv"(1) = 1$,
  $"inv"(i) = (M - floor(M\/i)) dot "inv"(M mod i) mod M$. To read
  a residue $r$ as a fraction, find a small $q$ with $r q mod M$
  small.

#full[
=== Concentration
Markov: $P(X >= a) <= E[X]\/a$ for $X >= 0$. Chebyshev:
$P(|X - mu| >= k sigma) <= 1\/k^2$. Chernoff for a sum of
independent $[0, 1]$ variables with mean $mu$:
$P(X >= (1 + delta) mu) <= e^(-delta^2 mu\/3)$ ($0 < delta <= 1$),
$P(X <= (1 - delta) mu) <= e^(-delta^2 mu\/2)$. Hoeffding:
$P(X - E X >= t) <= e^(-2t^2\/n)$, also for $n$ draws without
replacement from a set of numbers in $[0, 1]$. Azuma: a martingale
with $|M_i - M_(i - 1)| <= c_i$ has $P(M_n - M_0 >= t) <=
e^(-t^2 \/ (2 sum c_i^2))$. Union bound:
$P(union.big A_i) <= sum P(A_i)$. A test that fails with
probability $1\/2$ repeated 60 times fails with probability
$2^(-60)$.
- So a random process of $n$ bounded steps rarely strays more
  than a few $sqrt(n)$ from its expected path: a DP over (step,
  state) may skip the states outside that band. Grow the band
  until the answer stops changing (400 balls $+1$ and 250 balls
  $-1$ drawn without replacement, stop at will: error $10^(-7)$
  with half-width 40).
- Total variation $d(P, Q) = max_A |P(A) - Q(A)| = 1/2 sum_x
  |p_x - q_x|$. One sample from $P$ or $Q$ (prior $1\/2$ each) is
  classified correctly with probability at most $(1 + d)\/2$
  (guess the likelier source). $d <= P(X != Y)$ whenever
  $X tilde P$ and $Y tilde Q$ live on one probability space
  (coupling); some coupling attains it.

=== Continuous
- A walk of $n$ steps of length $ell$ in uniformly random
  directions in the plane is close to a Gaussian with independent
  coordinates of variance $n ell^2 \/ 2$ (CLT). A centred 2D
  Gaussian with coordinate variance $sigma^2$ lies within distance
  $R$ of its centre with probability $1 - e^(-R^2 \/ (2 sigma^2))$.
  With $ell^2 = dif t$ (Brownian motion, coordinate variance
  $t\/2$ at time $t$) the expected time within $R$ of the start
  up to $T$ is $integral_0^T (1 - e^(-R^2\/t)) dif t =
  T (1 - e^(-R^2\/T)) + R^2 E_1 (R^2\/T)$, and the time spent per
  unit area at distance $r$ is $E_1 (r^2\/T)\/pi$ (a radial density
  for Green's theorem, Geometry facts), with
  $E_1 (x) = integral_x^infinity e^(-u)\/u dif u$#cpp[
  (`-expint(-x)`)].
- Biased ruin duration ($p != 1\/2$, $q = 1 - p$):
  $i\/(q - p) - N\/(q - p) dot (1 - (q\/p)^i)\/(1 - (q\/p)^N)$.
- Absorbing chain with transient part $Q$ and $R$ to absorbing
  states: $(I - Q)^(-1)$ holds expected visits, its row sums the
  expected steps, $(I - Q)^(-1) R$ the absorption probabilities.
- Distribution after $t$ steps: $pi_0 P^t$ by matrix power.

=== Recognition
- A uniformly random process often hides a uniform permutation:
  in a shuffled deck every order of any $k$ given cards is equally
  likely; random removals leave a uniform random subset.
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
- Several independent processes, one advanced per step, first to
  finish ends the game: index policy.
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
