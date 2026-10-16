#import "/typst/notes.typ": *

#full[
Tools: `poly/ntt` (products), `poly/poly-ops` (inverse, log, exp,
power, division), `poly/lagrange`, `poly/linear-rec`,
`poly/subset-transforms`. Modelling with series is in the
generating functions note; this note is about computing them.
]

=== NTT primes and sizes
#table(
  columns: 4,
  table.header[prime $p$][$p - 1$][longest NTT][root],
  [998244353], [$119 dot 2^23$], [$2^23$], [3],
  [167772161], [$5 dot 2^25$], [$2^25$], [3],
  [469762049], [$7 dot 2^26$], [$2^26$], [3],
  [754974721], [$45 dot 2^24$], [$2^24$], [11],
  [1004535809], [$479 dot 2^21$], [$2^21$], [3],
)
#full[
- A primitive $2^k$-th root of unity is $g^((p - 1) \/ 2^k)$ for
  the root $g$ above. The modulus 998244353 in a statement is the
  setter's hint that NTT is welcome.
]
- Any other modulus ($10^9 + 7$): `poly/fft-mod`, or three primes
  from the table and CRT, exact while $n m^2 < 2^85$ (product of
  the 2nd to 4th primes).
#full[
- The length must be a power of two $>= |a| + |b| - 1$; shorter
  wraps around (cyclic). Reduce negative values first.
]
#full[

=== Divide-and-conquer products
#fact[Many factors][For $P_1 P_2 dots P_m$ of total degree $D$
  multiply in rounds, two at a time: every coefficient takes part
  in $O(log m)$ products, $O(D log D log m)$ in all. Multiplying
  left to right is $O(D^2)$.]
#lang-code(
  cpp: "deque<vl> q(all(P)); // multiply all polynomials in P
while (sz(q) > 1) {
  q.push_back(conv(q[0], q[1]));
  q.pop_front(), q.pop_front();
}",
  py: "q = deque(P)  # multiply all polynomials in P
while len(q) > 1:
  q.append(conv(q.popleft(), q.popleft()))",
)
- Uses: $product_i (1 + y_i x)$ gives all elementary symmetric
  sums $e_k$; one row of Stirling numbers of the first kind is
  $x (x + 1) dots (x + n - 1)$; "each item adds $a_i x^(s_i)$ or
  $b_i x^(t_i)$".
- A DP with $s$ states where each element is skipped or taken
  ($times x$) is a product of $s times s$ matrices of polynomials.
  Multiply them by halves; in a merge, `ntt` all $2 s^2$ entries
  once at a common length, combine pointwise and invert $s^2$
  results.
- Weighted sums of prefix products $sum_t c_t product_(j < t) P_j$:
  each block returns its product $F$ and its weighted sum $H$;
  merge $F = F_L F_R$, $H = H_L + F_L H_R$. $O(D log^2 D)$ for
  total degree $D$.

]

=== Online (relaxed) convolution
When $f_t$ needs $"acc"_t = sum_(0 < i <= t) g_i f_(t - i)$ and then
a non-linear step (a min, a per-index fix-up), $f$ cannot be
multiplied in one go. Split $[l, r)$ at $m$: finish $[l, m)$, add
$f[l, m) times g[0, r - l)$ into $"acc"$ on $[m, r)$, then do
$[m, r)$. Each pair $(i, t)$ meets once, where they are split:
$O(n log^2 n)$.#py[ Depth $log n$, so recursion is fine in PyPy.]
#lang-code(
  cpp: "void cdq(int l, int r) { // f, g, acc: vl of length n
  if (r - l == 1) return void(f[l] = finish(l, acc[l]));
  int m = (l + r) / 2;
  cdq(l, m);
  vl c = conv(vl(f.begin() + l, f.begin() + m),
              vl(g.begin(), g.begin() + r - l));
  rep(t, m, r) acc[t] = (acc[t] + c[t - l]) % MOD;
  cdq(m, r);
}",
  py: "def cdq(l, r):  # f, g, acc: lists of length n
  if r - l == 1:
    f[l] = finish(l, acc[l])
    return
  m = (l + r) // 2
  cdq(l, m)
  c = conv(f[l:m], g[: r - l])
  for t in range(m, r):
    acc[t] = (acc[t] + c[t - l]) % MOD
  cdq(m, r)",
)
- A linear step $f_t = h_t + "acc"_t$ needs no recursion:
  $f = h \/ (1 - g)$ after setting $g_0 = 0$ ($"acc"$ never reads
  $g_0$), one `poly_inv` and one `conv`.
- Both factors unknown ($f_t = sum_(i + j = t - 1) f_i f_j$): at
  $l = 0$ multiply $f[0, m)$ by itself; for $l > 0$ multiply
  $f[l, m)$ by $f[0, r - l)$ and add twice. Pad $n$ to a power of
  two so that $r - l <= l$ there.

#full[
=== Answer polynomial in a parameter
`poly/lagrange` (Notes) covers a DP over values $0, dots, K$ whose
answer is a polynomial in $K$, its degree, and a DP dimension that
only adds under merges replaced by a variable $z$. Further:
- Lattice points in the $t$-times dilated lattice polytope: a
  polynomial in $t$ of degree = dimension (Ehrhart).
- With that $z$ (a component size, a total weight) and an answer
  $sum_s q_s ["count with size" s]$, $s <= n$: precompute weights
  $w_j$ with $sum_j w_j z_j^s = q_s$; the answer is
  $sum_j w_j "DP"(z_j)$, no interpolation. A 2-D tree knapsack in
  $O(n^4)$ becomes $O(n^3)$.
- A rational answer $N(x) \/ D(x)$ with a known denominator (such
  as $product_i (1 - i x)^(c_i)$) and $deg N <= d$: run the scalar
  DP at $d + 1$ points $x_0$ that are not poles, multiply by
  $D(x_0)$, `interpolate` $N$, then expand $N \/ D$ by dividing by
  each factor $1 - i x$: $s_n <- s_n + i s_(n - 1)$ for increasing
  $n$.
- A fixed linear transition repeated $N$ times (fixed matrix or
  automaton) gives a linear recurrence instead:
  `poly/linear-rec`.

=== Coefficient windows
- Middle product: for $a$ with $n$ and $b$ with $m >= n$
  coefficients, a cyclic convolution of any length $L >= m$ gives
  the coefficients $n - 1, dots, m - 1$ of $a b$ exactly (wrapped
  terms land below $n - 1$), at half the NTT length of the full
  product.
- $[y^n] E(y) B(y)^n$ for every $n < N$, $deg B = d$: divide and
  conquer over $n$. For $[l, r)$ only the coefficients of $E B^l$
  with exponents in $[d l - (d - 1) r, r)$ matter. Recurse left on
  the sub-window for $[l, m)$; for the right half multiply the
  window by $B^(m - l)$ (a middle product) and keep
  $[d m - (d - 1) r, r)$. $O(N log^2 N)$ for constant $d$; for
  $B = (1 + y)^2$ the powers are binomials.
- $[x^(c_i)]$ of a prefix times a suffix product for every $i$:
  All but one (Offline techniques), carrying only the window of
  coefficients that the leaves below can still read.

=== Sums as convolutions
- Binomial sums: expand into factorials and split into a factor
  in $i$, one in $k - i$ and one in $k$:
  $sum_i p_i binom(n - i, k - i) = 1/((n - k)!) sum_i (p_i (n - i)!) dot 1/((k - i)!)$,
  one `conv` for all $k$.
- Correlations $sum_i a_i b_(i + d)$: reverse $a$ (`poly/ntt`
  Notes). Cyclic $sum_i a_i b_((i + d) mod n)$: also concatenate
  $b$ with itself; $d$ is then at index $n - 1 + d$.

=== Transforms on bit masks
#table(
  columns: 3,
  table.header[condition][per bit][tool],
  [sum over $T subset.eq S$], [$(u, v) -> (u, u + v)$],
  [`'|'` zeta],
  [sum over $T supset.eq S$], [$(u, v) -> (u + v, v)$],
  [`'&'` zeta],
  [exactly $S$ from within $S$], [$(u, v) -> (u, v - u)$],
  [`'|'` inverse (Möbius)],
  [$sum_T (-1)^(|S inter T|) a_T$], [$(u, v) -> (u + v, u - v)$],
  [`'^'` (Hadamard)],
  [pairs with `i op j == S`], [transform, multiply, invert],
  [`bit_conv`],
  [disjoint $i union j = S$], [zeta per rank, product in rank],
  [`subset_conv`],
)
- Sizes: $n = 20$ gives $n 2^n = 2 dot 10^7$ (easy) and
  $n^2 2^n = 4 dot 10^8$ (subset convolution: 1 s and 200 MB in
  C++, 3.5 s and 450 MB in PyPy). Enumerating all submasks of all
  masks is $3^n$: $4 dot 10^7$ at $n = 16$, $3.5 dot 10^9$ at
  $n = 20$.
- The XOR inverse divides by $2^n$: fine modulo an odd prime,
  impossible modulo $2^64$; with exact integers the division is
  exact if the values fit.
- Wildcard patterns over ${0, 1, ?}^n$: per coordinate
  $(x_0, x_1) -> (x_0, x_1, x_0 + x_1)$ gives all $3^n$ subcube
  sums in $O(3^n)$.
- XOR moves that depend only on the weight ("flip exactly $a$ of
  $n$ bits"): the transform depends only on the weight $j$ of the
  frequency, the Krawtchouk number
  $K_a (j) = [t^a] (1 + t)^(n - j) (1 - t)^j$ with
  $(a + 1) K_(a + 1) = (n - 2j) K_a - (n - a + 1) K_(a - 1)$. Work on
  the $n + 1$ weights in $O(n^2)$: the number of ways to get from
  0 to a fixed $y$ of weight $r$ with moves $a_1, a_2, dots$ is
  $2^(-n) sum_j K_j (r) product_i K_(a_i) (j)$.
]
