#import "/typst/notes.typ": *
#let xor = sym.plus.o

=== Single bits
#lang-code(
  cpp: "x & -x                 // lowest set bit, as a number
x & (x - 1)            // x without its lowest set bit
x >> i & 1             // bit i; set it with x | 1LL << i
__builtin_popcountll(x)  // set bits (without ll: 32 bits)
__lg(x)                // index of the highest bit; x > 0
__builtin_ctzll(x)     // index of the lowest bit; x > 0
bit_width(ull(x))      // bits needed, 0 for 0 (<bit>)
bit_ceil(ull(x))       // smallest power of 2 >= x
countr_zero(ull(x))    // like ctz, but 64 for x = 0",
  py: "x & -x                 # lowest set bit, as a number
x & (x - 1)            # x without its lowest set bit
x >> i & 1             # bit i; set it with x | 1 << i
bin(x).count(\"1\")      # set bits (int.bit_count is 3.10+)
x.bit_length() - 1     # index of the highest bit; x > 0
(x & -x).bit_length() - 1  # index of the lowest bit; x > 0
x.bit_length()         # bits needed, 0 for 0
1 << (x - 1).bit_length()  # smallest power of 2 >= x >= 1",
)
#cpp[The builtins are undefined for 0 (except popcount), and the
  `<bit>` functions take only unsigned types (`ull` is
  `unsigned long long`). `1 << 40` overflows `int`: write
  `1LL << k`.]
#py[PyPy ints are machine words below $2^63$ and big ints above,
  about 10 times slower: keep masks within 62 bits. Negative ints
  act as infinite two's complement: `-1 >> k == -1`,
  `~x == -x - 1`.]

=== Loops over sets
#lang-code(
  cpp: "for (int s = m; s; s = (s - 1) & m) {} // submasks > 0 of m
for (int s = m; s < 1 << n; s = (s + 1) | m) {} // supersets
for (ll t = x; t; t &= t - 1) {
  int i = __builtin_ctzll(t); // each set bit i of x
}
ll next_comb(ll x) { // next larger with the same popcount
  ll c = x & -x, r = x + c;
  return (((r ^ x) >> 2) / c) | r;
} // the k-subsets of n bits in increasing order, k >= 1:
for (ll x = (1LL << k) - 1; x < 1LL << n; x = next_comb(x)) {}",
  py: "s = m                # submasks > 0 of m, decreasing
while s:
  ...                # use s
  s = (s - 1) & m
s = m                # supersets of m within n bits
while s < 1 << n:
  ...
  s = (s + 1) | m
t = x                # each set bit i of x
while t:
  low = t & -t
  i = low.bit_length() - 1
  t ^= low
def next_comb(x):    # next larger with the same popcount
  c = x & -x
  r = x + c
  return ((r ^ x) >> 2) // c | r",
)
Submasks of all masks together: $3^n$ steps. Sums over all
submasks or supersets: `poly/subset-transforms`.#py[ All
$k$-subsets: `itertools.combinations(range(n), k)`.]

=== Identities
- $a + b = (a xor b) + 2 (a "&" b)$: XOR adds without carries, so
  $a + b = a xor b$ iff $a "&" b = 0$. For $a >= b >= 0$:
  $a - b <= a xor b <= a + b$.
- $0 xor 1 xor dots xor n$ is $n, 1, n + 1, 0$ for $n mod 4 = 0,
  1, 2, 3$; the XOR of $[l, r]$ is $f(r) xor f(l - 1)$.
- The smallest XOR of two elements of a set is attained by two
  neighbours in sorted order. Largest XOR pair: put the numbers in
  a binary trie and walk each one down from the top bit, taking
  the opposite bit when it exists#full[ (`ds2/binary-trie`)].
  Largest XOR of a subset: `linalg/xor-basis`.
- Gray code $g(i) = i xor floor(i \/ 2)$ (`i ^ i >> 1`): $g(i)$ and
  $g(i + 1)$ differ exactly in the lowest set bit of $i + 1$, and
  $g(2^n - 1)$ and $g(0)$ in bit $n - 1$, so the codes of
  $0, dots, 2^n - 1$ form a cycle. Inverse: $i$ is the XOR of
  $floor(g \/ 2^j)$ over all $j >= 0$.
- $x xor (x + 1) = 2^(t + 1) - 1$, where $t$ is the number of
  trailing 1 bits of $x$: adding 1 clears them and sets the next
  bit.
#cpp[
  - Precedence: `x & 1 == 0` means `x & (1 == 0)`; shifts bind
    looser than `+`: `1 << n - 1` is `1 << (n - 1)`.
]
#py[
  - Precedence: comparisons bind looser than `&`, so `x & 1 == 0`
    works; shifts bind looser than `+`: `1 << n - 1` is
    `1 << (n - 1)`.
]

#full[
=== More facts and cues
- Gray code, recursively: $0$ followed by $G_(n - 1)$, then $1$
  followed by $G_(n - 1)$ reversed.
- "XOR each vertex with its children", $t$ times, i.e.
  $(I + C)^t$ over GF(2): _Tree techniques_.
- For a fixed right end $r$, the ORs (or ANDs, or GCDs) of
  $a_l, dots, a_r$ over all $l$ take $O(log V)$ distinct values
  (at most $b + 1$ for values below $2^b$): keep the list of
  distinct values with their leftmost $l$ and update it per
  element, $O(n log V)$ in total.
- Items of sizes $2^j$ ($a_j$ of each) fit into bins of capacities
  $H_i$ iff for every $l$:
  $sum_i floor(H_i \/ 2^l) >= sum_(j >= l) a_j 2^(j - l)$.
  Largest first, into any bin with room, is optimal: powers of two
  up to $2^l$ with total at least $2^l$ contain a subset summing to
  exactly $2^l$.
- "Choose numbers so that a bitwise expression is extreme":
  decide bits from the highest down greedily, keeping the
  candidates (or a trie node) that can still reach the prefix.
]
