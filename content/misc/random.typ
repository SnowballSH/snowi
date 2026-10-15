#import "/typst/notes.typ": *
#let xor = sym.plus.o

=== Random numbers
#lang-code(
  cpp: "mt19937_64 rng(chrono::steady_clock::now()
                     .time_since_epoch().count());
uint64_t r = rng();                            // 64 random bits
ll x = uniform_int_distribution<ll>(lo, hi)(rng); // [lo, hi]
double u = uniform_real_distribution<>(0, 1)(rng);
shuffle(all(v), rng);",
  py: "import random  # seeded from the OS at import
r = random.getrandbits(62)   # 62 bits stay a machine int
x = random.randint(lo, hi)   # [lo, hi]
random.shuffle(a)",
)
#cpp[Avoid `rand()` (31 bits here, 15 on some systems) and
  `random_shuffle` (removed in C++17).] Seed with a constant while
debugging so that failures repeat. Where tests may target you
(hacks; tests built against unseeded hashes or first-element
pivots), use an unpredictable seed#cpp[ (the clock)]#py[ (the
default)].

=== Fingerprints
Give every value $v$ an independent random label $R_v$ of $b =$
#cpp[64]#py[62] bits. Two different objects get equal fingerprints
with probability $2^(-b)$; among $n$ objects some pair collides
with probability about $n^2 \/ 2^(b + 1)$, so 32-bit labels
already fail for $10^5$ objects.
#fact[Sets][XOR of the labels. It sees the values that occur an
  odd number of times; adding or removing a value is one XOR.
  Zobrist hashing of a board: XOR of `R[cell][piece]`.]
#fact[Multisets][Sum of the labels#cpp[ in `uint64_t` (it wraps
  mod $2^64$)]#py[; with 62-bit labels keep sums in a machine
  word with `& ((1 << 62) - 1)`]. It keeps multiplicities: `a[l:r]`
  is a permutation of $1, dots, r - l$ iff its sum equals
  $R_1 + dots + R_(r - l)$ (prefix sums for both). If $2^t$ is the
  largest power of 2 dividing all differences of multiplicities,
  a collision has probability $2^(t - b)$.]
#fact[All or nothing]["Every value has all of its occurrences
  inside the range, or none": give each occurrence but the last a
  random label and the last the XOR of the others. A range is
  good iff its XOR is 0, i.e. `pre[l] == pre[r]` (w.h.p.). On a
  tree, path $u$-$v$ has XOR $f(u) xor f(v)$ for root prefixes
  $f$.]
#lang-code(
  cpp: "vector<uint64_t> h(n), acc(V); // a[i] in [0, V)
vi left(V);
for (int v : a) left[v]++;
rep(i, 0, n) {
  int v = a[i];
  h[i] = --left[v] ? rng() : acc[v];
  acc[v] ^= h[i];
}",
  py: "h, acc, left = [0] * n, [0] * V, [0] * V  # a[i] in [0, V)
for v in a:
  left[v] += 1
for i, v in enumerate(a):
  left[v] -= 1
  h[i] = random.getrandbits(62) if left[v] else acc[v]
  acc[v] ^= h[i]",
)
#fact[Counts divisible by $k$][Choose a random half of the values
  and count the range's elements in it. If every count is
  divisible by $k$, so is this number; otherwise it is not with
  probability at least $1 \/ 2$. Repeat 30 times, with one prefix
  count array (or Fenwick tree) per trial.]
Count equal fingerprints by sorting them or in a hash map:
random keys are safe. Keys taken from the input need a seeded
hash#cpp[ (`ds2/hash-map`)]#py[: store `x ^ R` for one random
62-bit `R`].

=== Polynomial identity testing
#fact[Schwartz-Zippel][A nonzero polynomial of total degree $d$
  is 0 at a point with coordinates uniform in $S$ with probability
  at most $d \/ |S|$. Evaluate modulo a prime near $10^9$ or
  $2^61 - 1$, never modulo 2 or $2^64$.]
- Freivalds: to check $A B = C$, compare $A (B x)$ with $C x$ for
  a random vector $x$ mod $p$, in $O(n^2)$. A wrong $C$ passes
  with probability at most $1 \/ p$.
- Matching: the Tutte matrix ($T_(u v) = x_(u v)$,
  $T_(v u) = -x_(u v)$ per edge, random $x$ mod $p$) has rank
  twice the maximum matching size (`linalg/gauss-mod`).
- String hashing is the same lemma: two different strings of
  length $n$ collide for at most $n$ bases mod $p$, so pick the
  base at random (`strings/string-hash`).
- "Does a structure exist": find a determinant or rank that is
  nonzero iff it does, then substitute random values.

=== Randomised algorithms
- Frequent value (more than $n \/ k$ copies) in a range: take $t$
  random positions as candidates and count each exactly (sorted
  position lists per value, binary search). All miss with
  probability below $(1 - 1 \/ k)^t$.
- Witnesses: if a fraction $alpha$ of the candidates prove the
  answer, $t$ random samples all miss with probability
  $(1 - alpha)^t <= e^(-alpha t)$.
- Shuffle first so that bad orders become unlikely: incremental
  minimum enclosing circle in expected $O(n)$, insertion into
  unbalanced search trees, random pivots.
- Random order with skipping: run each candidate $u$ in random
  order unless an earlier run covered it; run $u$ covers
  $s(u) in.rev u$. If $w in s(u)$ implies $s(w) subset.eq s(u)$
  and any two sets are nested or disjoint, an element lying in $t$
  sets is visited by $H_t <= 1 + ln t$ runs in expectation (Kahn
  closures from every start vertex in $O((n + m) log n)$).
#full[
- Colour coding: to find $k$ objects that must all differ (a path
  through $k$ distinct vertices), colour everything with $k$
  random colours and search only for solutions with distinct
  colours, by a DP over subsets of colours. A fixed solution
  survives with probability $k! \/ k^k >= e^(-k)$, so
  $c dot k^k \/ k!$ rounds fail with probability at most $e^(-c)$.
- Las Vegas: retry until a random choice works (random vectors
  that form a basis, see Random matrices in _Linear algebra
  facts_; a random good pivot), checking each result exactly.
- Random answers from the judge: keep the posterior of every
  hypothesis (or per candidate, the number of negative answers)
  and ask queries whose answer is far from certain. If an answer
  carries at most $h$ bits of entropy, telling apart $N$ equally
  likely hypotheses takes at least $log_2 N \/ h$ queries on
  average. Never rule a candidate out on one noisy answer; tune
  the constants with a local simulator.
- Heuristics (random restarts, simulated annealing): keep the
  best answer and loop while the time used is below 90% of the
  limit `TL` in seconds:
  #cpp[`clock() < 0.9 * TL * CLOCKS_PER_SEC`]#py[`time.perf_counter() - start < 0.9 * TL`].

=== Anti-hash and pitfalls
- With hacks, hash tables keyed by input values are attackable.
  #cpp[`unordered_map<ll, ...>` with the default hash puts crafted
  keys into one bucket: use `ds2/hash-map`.]#py[Python dicts and
  sets hash an int to itself, so crafted keys can make them slow:
  store `x ^ R` for a random 62-bit `R`, or sort instead.]
- Polynomial hashing mod $2^64$ is broken for every base: odd
  bases by the Thue-Morse string of length 2048 and its
  complement, even bases by any two strings with the same last 64
  characters. Use a prime modulus with a random base.
- XOR loses multiplicity ($a xor a = 0$): use sums for multisets.
  For "set" semantics with repeated insertions, add a value only
  on its first occurrence.
- Know which side can err: a fingerprint can wrongly say
  "equal", never "different". Do exact checks (lengths, sizes)
  separately, and bound the total error over all comparisons by
  a union bound.
- With many test cases, draw labels only for the values used,
  not a fresh table of $10^6$ labels per test.
]
