#import "/typst/notes.typ": *

Strings are 0-indexed, `s[l, r)` is a substring and $n = |s|$. A
_border_ is a proper prefix that is also a suffix; $p$ is a
_period_ if `s[i] == s[i + p]` for all $i < n - p$.

=== Borders and periods
#fact[Duality][$p$ is a period iff $n - p$ is a border. The
  smallest period is `n - pi[n - 1]`; the periods are $n - b$ for
  the borders `b = pi[n - 1]`, `pi[b - 1]`, ..., 0. With the Z
  array: $p < n$ is a period iff `z[p] == n - p`; $p = n$ always
  is.]
#fact[Primitive root][With smallest period $p$: if $p | n$, `s` is
  `s[0, p)` repeated $n\/p$ times and that root is primitive;
  otherwise `s` is not a power at all ("abaab": $p = 3$). `s`
  equals its rotation by $k$ iff it is a power of
  `s[0, gcd(k, n))`, so a primitive `s` differs from all its other
  rotations.]
#fact[Fine-Wilf][Periods $p, q$ with $p + q - gcd(p, q) <= n$ make
  $gcd(p, q)$ a period. So every period $q <= n - p$ is a multiple
  of the smallest period $p$. Tight: "abaaba" has periods 3 and 5
  but not 1.]
#fact[Gluing][Two windows with period $d$ that share at least $d$
  positions: their union has period $d$ (sharing $d - 1$ is not
  enough).]
#fact[Progressions][Borders of length $>= n\/2$ form one
  arithmetic progression $n - p, n - 2p, ...$, so the border chain
  splits into $O(log n)$ progressions. Two occurrences of `t` at
  distance $d < |t|$ make $d$ a period of `t`; the occurrences in a
  window of length $2|t|$ form one progression.]
#fact[Infinite powers][$u^infinity = v^infinity$ iff $u v = v u$
  iff both are powers of one word. $u^infinity < v^infinity$ iff
  $u v < v u$: sort by `a + b < b + a` to get the smallest
  concatenation.]

=== Prefix function
#fact[Failure tree][Nodes $0..n$; the parent of $L$ is
  `pi[L - 1]`. The ancestors of $L$ are the borders of `s[0, L)`,
  the occurrences of each prefix are subtree sums, and the longest
  common border of prefixes $a$ and $b$ is the LCA of
  `pi[a - 1]` and `pi[b - 1]`. While matching a text, the
  prefixes of `s` ending at the current position are the current
  state and its ancestors. Parents are smaller: loop over $L$
  instead of a DFS.]

=== Palindromes
#fact[Count][At most $n$ distinct palindromic substrings: each end
  position adds at most one new one, the longest palindromic
  suffix. With repeats: the sum of $floor((p + 1)\/2)$ over the
  lengths $p$ in Manacher's `plen`.]
#fact[Borders][The borders of a palindrome are its palindromic
  suffixes; the lengths of the palindromic suffixes of any string
  form $O(log n)$ progressions (eertree series links).]
#fact[Prefixes][Palindromic prefixes of lengths $y < x$ make
  $x - y$ a period of `s[0, x)`. If $x < 2y$, the prefix of length
  $2y - x$ is a palindrome too and may replace $x$: positions
  $y..x - 1$ copy earlier ones. If $x >= 2y$, prefix $y$ lies in
  the first half and $x$ adds $ceil(x\/2) - y$ free positions.
  Repeat on the two longest, batching equal steps like Euclid.]

=== Lyndon words
#fact[Factorization][`s` splits uniquely into Lyndon words
  $w_1 >= w_2 >= dots >= w_k$ in $O(n)$ (Duval). The last factor is
  the smallest suffix; the smallest rotation starts at a factor of
  `s + s`. A Lyndon word has no non-empty border, so it is
  primitive.]

=== Runs
#fact[Runs][A run is a maximal window `s[i, j)` whose smallest
  period $p$ has $j - i >= 2p$; every square or higher power lies
  in one, and a run holds a square of length $2p$ at each start
  in `[i, j - 2p]`. There are fewer than $n$ runs, with exponents
  $(j - i)\/p$ summing to less than $3n$. All runs in
  $O(n log n)$: for each $p$ and each $a = k p$ with $b = a + p < n$,
  take `s[a - lcs, b + lcp)` with `lcp` of the suffixes at $a, b$
  and `lcs` of the prefixes ending there; keep it if its length is
  $>= 2p$ and $p$ is its smallest period, once per distinct window
  (several $a$ find the same run). $O(log n)$ squares with
  a primitive root start at one position.]

=== Suffix array and LCP
#fact[LCP][With `lcp[i]` the LCP of the suffixes ranked $i$ and
  $i + 1$, as in the snippet, the LCP of the suffixes ranked
  $a < b$ is the minimum of `lcp[a, b)` (sparse table). Distinct
  substrings: $n(n+1)\/2 - sum "lcp"[i]$. Longest repeated
  substring: max `lcp`. Longest common substring of `s` and `t`:
  suffix array of `s + '#' + t`, best `lcp` between neighbours
  from different sides.]

=== Hashing
#fact[Safety][A random base modulo $2^61 - 1$ fails with
  probability at most $n\/M$ per comparison. Overflow modulo
  $2^64$ is broken by Thue-Morse strings; a single modulus near
  $10^9$ fails once about $10^5$ hashes meet (birthday bound).
  Compare lengths too.]

=== Wildcards and mismatches by FFT
#fact[Wildcards][Letters are $1..26$, a wildcard is 0. The
  alignment at $i$ matches iff
  $sum_j p_j t_(i+j) (p_j - t_(i+j))^2 = 0$: expand into three
  convolutions with `p` reversed. Doubles lose exactness for large
  $n$; use NTT modulo a prime $q$ with random letter values
  instead (a false match has probability at most $4\/q$).
  Mismatch counts for all alignments: one 0/1 convolution per
  letter.]

#full[
=== Cues
- "Prefix equals suffix", the same piece at both ends: borders,
  failure tree.
- "Made of copies", cyclic strings, "the first $ell$ characters
  repeated": smallest period, Fine-Wilf.
- Several palindromic prefixes: they are borders, hence periods.
- Smallest rotation, suffix or cut: Lyndon factorization.
- Every suffix against a pattern: Z of `t + '#' + s`.
- Many "equal substrings?" queries: hashing; their order: LCP by
  hashing or a suffix array.
- Squares, cubes, patterns like `r r r t t`: runs, primitive
  roots.
- "Delete one letter", "two different index sets": runs of equal
  letters, embeddings that differ in one index.

=== More tools
- Fine-Wilf, weak form ($p < q$, $p + q <= n$): $q - p$ is a period
  (step forward by $q$ and back by $p$, or back by $p$ and forward
  by $q$ near the end). Repeating it is subtractive Euclid.
- Lengths up to $10^18$ built from repeated prefixes shrink like
  Euclid: replace runs of equal subtractions by one division.
- Borders of `s[0, L)` no longer than $L\/2$: keep a second KMP
  pointer, extend it like matching and pull it back with
  `while (2 * k > L) k = pi[k - 1]`. Their number is the depth of
  `k` in the failure tree (NOI 2014 "Zoo").
- The longest Lyndon word starting at $i$ ends just before the
  first later suffix that is smaller than suffix $i$: next smaller
  value over suffix ranks with a monotone stack.
- Lyndon words of length $n$ over $k$ letters:
  $1\/n sum_(d | n) mu(d) k^(n\/d)$ (necklaces: _Counting
  formulas_).
- $k$-th smallest distinct substring, occurrences of every
  substring: suffix array with `lcp`, or the suffix automaton.
- Palindromic suffixes $P, Q$ of one prefix with $|P| < |Q|$: the
  letter just before $P$ is $Q[abs(P)]$. So the letter before
  $"link"(Q)$ depends only on $Q$ and can be stored at $Q$ (quick
  links).
- Subsequences: greedy matching (each letter of `t` at its first
  possible place in `s`) is the leftmost embedding in every index.
  If `t` embeds in two ways, two of its embeddings differ in
  exactly one index.
- Strings read along pointers (tree root paths, next-greater
  links): Upward sequences in _Tree techniques_.
- Deleting one letter: every letter of a run of equal letters
  gives the same string. Deleting from a run followed by a smaller
  letter beats deleting from any later run; from a run followed
  by a larger letter (or the last run) it loses to them. The
  smallest result deletes from the first run followed by a
  smaller letter, else from the last run.
- Strings over $m$ letters with a given `pi`: build the generic
  string (a new letter where `pi[i] = 0`, else `s[pi[i] - 1]`);
  if its `pi` differs, the answer is 0. Otherwise it is
  $m product_(i >= 1, "pi"[i] = 0) max(0, m - z_i)$, where $z_i$
  counts the borders $j$ of `s[0, i)`, 0 included, with
  `pi[j] = 0`: the distinct letters `s[j]` that position $i$ must
  avoid.
- Longest balanced substring of a bracket string, with prefix
  sums $P_0..P_n$ (`(` is $+1$): the largest of $j - 1 - i$ over
  consecutive strict prefix-minimum records $i < j$ from the
  left, the same over records from the right, and the distance
  between the first and the last position of $min P$. For
  substrings, follow next-smaller chains with binary lifting.

=== Pitfalls
- `pi[i]` belongs to the prefix ending at $i$, a border length $L$
  to the prefix of length $L$: the chain step is `k = pi[k - 1]`.
- Some statements count the whole string as a border.
- The smallest period $p$ need not divide $n$; `s` is a power
  $u^k$ with $k >= 2$ iff $p < n$ and $p | n$.
- With text appended, the smallest suffix is not always the best
  start: `a` is the smallest suffix of `aba`, yet `abac < ac`.
- Recursion over the failure tree of `aaa...a` is $n$ deep.
- Characters with value 0 make hashes of different lengths equal;
  key by (hash, length) or map letters to $1..26$.
]
