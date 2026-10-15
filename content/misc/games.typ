#import "/typst/notes.typ": *
#let xor = sym.plus.o

Two players alternate, both know everything, and whoever cannot
move loses (normal play). Impartial: both players have the same
moves from every position.

=== Winning and losing positions
- A position is losing iff every move leads to a winning one; a
  terminal position is losing. DP over positions in an order where
  moves go to smaller states.
- With cycles, run a BFS backwards from the terminal positions:
  a position becomes winning once one successor is losing, and
  losing once all its successors are winning (count them down).
  Positions never decided are draws.
- Scoring games: with $f(s)$ = (my score) $-$ (opponent's score)
  from $s$ on, $f(s) = max_(s -> s') ("gain" - f(s'))$.

=== Sprague-Grundy
#fact[Theorem][For an impartial game,
  $g(x) = "mex" {g(y) : x -> y}$ (the smallest value $>= 0$ not
  among the moves). A sum of games, where a move changes one
  component, has $g = g_1 xor dots xor g_k$, and the player to
  move loses iff $g = 0$.]
- A move that splits a component into independent parts leads to
  the XOR of their values.
- $g(x)$ is at most the number of moves from $x$: size mex arrays
  by that.
- Different moves for the two players (partizan games): no
  Grundy values; solve win/lose on the whole state.

=== Nim and relatives
#fact[Nim][The mover wins iff $s = a_1 xor dots xor a_n != 0$.
  Winning move: a pile with $a_i xor s < a_i$, reduced to
  $a_i xor s$.]
#fact[Misère Nim][Whoever makes the last move loses. If every pile
  has at most 1 stone, the mover wins iff the number of piles of
  size 1 is even; otherwise as in normal Nim ($s != 0$).]
#fact[Subtraction games][Take 1 to $k$: $g(n) = n mod thin (k + 1)$.
  For any finite set of allowed amounts, $g$ is eventually
  periodic: tabulate it and find the period.]
#fact[Staircase Nim][Coins on steps $0, 1, dots$; a move shifts
  coins from step $i > 0$ to step $i - 1$. It is Nim on the odd
  steps: moves from even steps are undone by the opponent.]
#fact[Moore's Nim][Up to $k$ piles may be reduced in one move. The
  mover loses iff for every bit, the number of piles with that bit
  set is divisible by $k + 1$.]
#fact[Wythoff][Two piles; take any amount from one pile or the
  same amount from both. The losing positions $(a, b)$ with
  $a <= b$ are $(floor(j phi), floor(j phi) + j)$ for $j >= 0$,
  $phi = (1 + sqrt(5)) \/ 2$ (sort the piles first); exactly,
  $floor(j phi) = floor((j + floor(sqrt(5 j^2))) \/ 2)$ with an
  integer square root.]
#fact[Green Hackenbush][On a rooted tree a move deletes an edge
  and everything no longer connected to the root.
  $g(v) = xor.big_c (g(c) + 1)$ over the children $c$ of $v$; a
  path of $n$ edges is worth $n$.]

=== Strategy arguments
- Mirror: on a symmetric position the second player copies the
  first through the centre (two equal piles, a symmetric board);
  if the centre is a legal move, the first player takes it and
  then mirrors.
- Pairing: split the moves into pairs and answer each move with
  its partner.
- Strategy stealing: if an extra move never hurts, the first
  player wins (Chomp, Hex). It shows who wins, not how.
- Fixed length: if every play has the same number of moves (the
  final state is forced, whatever the order), only its parity
  matters.
- Token on an undirected graph, visited vertices deleted: the
  mover wins iff the start vertex is in every maximum matching,
  i.e. $nu(G - s) < nu(G)$.

=== Patterns from brute force
Large inputs with a simple rule hide a period or a formula. Write
the brute force, print tables over each parameter, guess, then
check the guess against the brute force on every small case.
Reduce objects to a small invariant first: if only "is the
product a square mod $p$" matters, each number is one Legendre bit
($a^((p - 1) \/ 2) mod p = plus.minus 1$, multiplicative), so the
state is a few counts and a parity; tabulate those.
#lang-code(
  cpp: "// g[n] for n < N; opts(n, g) lists the values of the moves
template<class F> vi grundy(int N, F opts) {
  vi g(N);
  rep(n, 0, N) {
    vi s = opts(n, g);
    sort(all(s));
    for (int v : s) g[n] += v == g[n]; // mex
  }
  return g;
}
auto opts = [](int n, const vi& g) { // take 1, 3 or 4
  vi s;
  for (int d : {1, 3, 4}) if (d <= n) s.push_back(g[n - d]);
  return s; // splitting into a, b: push g[a] ^ g[b]
};
vi g = grundy(100, opts); // 0 1 0 1 2 3 2, period 7",
  py: "def grundy(N, opts):  # opts(n, g): values of the moves
  g = [0] * N
  for n in range(N):
    s = set(opts(n, g))
    while g[n] in s:
      g[n] += 1
  return g
# take 1, 3 or 4; splitting into a, b: g[a] ^ g[b]
opts = lambda n, g: [g[n - d] for d in (1, 3, 4) if d <= n]
g = grundy(100, opts)  # 0 1 0 1 2 3 2, period 7",
)

#full[
=== More
- Lasker's Nim (take any, or split a pile into two nonempty
  piles): $g(4 i + 1) = 4 i + 1$, $g(4 i + 2) = 4 i + 2$,
  $g(4 i + 3) = 4 i + 4$, $g(4 i + 4) = 4 i + 3$.
- Green Hackenbush with cycles: contract each 2-edge-connected
  part into one vertex; its $E$ edges become loops, together worth
  $E mod 2$ (fusion principle). Then apply the tree rule.
- Games with chance: value of a chance node = probability-weighted
  average; a cycle in the states gives an equation (solve it) or
  value iteration until changes are below $10^(-12)$.
- Hidden information that the opponent need not commit to: play
  against the set of positions still consistent with the moves.
- Recognition: "both play optimally" plus independent piles,
  subtrees or segments: Sprague-Grundy. Huge counts with a simple
  rule: tables and periods. A modification step before the game
  ("change the fewest cells so Alice wins"): characterise the
  winner exactly first, then optimise.
]
