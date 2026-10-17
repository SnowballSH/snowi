#import "/typst/notes.typ": *

#full[
When no clean $O(log n)$ structure exists, split the work into a
part that is cheap per item and a part that is cheap per group,
and balance the two. Pick the block size by equating the costs
($a n / B = b B$), then time a few sizes around it: C++ manages a
few times $10^8$ simple steps per second, PyPy a few times fewer.
Offline orders (sorting by a threshold, last occurrences, reversed
time) are in _Offline techniques_.

]
=== Blocks
#full[
#fact[Range update, range query][Cut the array into blocks of
  $B approx sqrt(n)$. An operation touches $O(n / B)$ whole blocks
  through a per-block summary and lazy tag, and $O(B)$ single
  elements in the two partial blocks, after which those blocks
  are rebuilt. Use it when a block summary is easy but a segment
  tree node is not (mode, values `>= x` under adds). A block may
  also tabulate its effect for each of the $O(B^2)$ contiguous
  ranges of its own elements or operations, when a query always
  uses such a range.]
]
#fact[Sorted blocks][Range add and "count values `>= x` in
  `[l, r)`": each block keeps a sorted copy and an add tag; a whole
  block answers by binary search for $x - "tag"$ in $O(log B)$. In
  a partial block, add to the affected elements in place; they
  and the other elements each stay in sorted order, so one merge
  rebuilds the sorted copy in $O(B)$. Then
  $B approx sqrt(n log n)$ gives $O(sqrt(n log n))$ per
  operation; a plain re-sort costs $O(B log B)$ and wants
  $B approx sqrt(n)$, for $O(sqrt(n) log n)$.]
#full[
#fact[Periodic rebuild][Keep a static structure (prefix sums,
  sorted array, a built tree) plus a list of the last $B$
  updates. A query asks the static part and scans the list; after
  $B$ updates rebuild in $O(n)$. Total $O(q B + q n / B)$, so
  $B = sqrt(n)$. Also works over queries in time order: process
  $B$ queries at a time, with only the $O(B)$ "active" elements
  handled by brute force.]
#fact[Heavy and light][Sum of degrees $2 m$ means at most
  $2 m / D$ vertices have degree $> D$. Heavy vertices get
  precomputed answers or are updated eagerly, light ones are
  scanned: neighbour sums under point updates in
  $O(sqrt(m))$ per operation; all triangles in $O(m sqrt(m))$
  (_Graph theory facts_). The same split works for values by
  frequency: at most $n / D$ values occur more than $D$ times.]
]
#fact[Small and large parameters][Queries with a step, modulus or
  length $d$: for $d <= sqrt(n)$ precompute a table per $d$
  ($O(n sqrt(n))$ in total), for larger $d$ walk the
  $n / d <= sqrt(n)$ positions directly. Updates along a
  progression `s, s + d, ..., s + (c - 1)d` with small $d$: in a
  difference array for $d$ add $x$ at `s` and $-x$ at `s + c d`,
  then propagate once, `D[i] += D[i - d]`, $O(n)$ per small $d$.
  If the weights grow by a factor $q$ per step, propagate with
  `D[i] += q * D[i - d]` and put $-x q^c$ at `s + c d`.]
#fact[Split the bits][Keys of $2 b$ bits where an update and a
  query must meet on all keys, and the cost splits into a high
  and a low half (XOR distance $(y xor x)$, subset or superset
  tests): an update fills one row `T[hi(x)][*]` over all
  $2^b$ low halves, a query scans all $2^b$ rows at its own low
  half. Example, min over stored $(x, "val")$ of
  $"val" + (y xor x)$: update
  `T[hi(x)][l] = min(.., val + (lo(x) ^ l))` for all `l`; query
  `min over h of T[h][lo(y)] + ((hi(y) ^ h) << b)`. Both cost
  $O(2^b) = O(sqrt(V))$ for keys below $V = 4^b$.]

=== Monotone stack and queue
#full[
#fact[Next greater element][Scan left to right with a stack of
  indices whose values decrease; an element pops everything
  smaller, and it is the answer for each popped index. Each index
  is pushed and popped once: $O(n)$. Gives, for each `i`, the
  maximal range where `a[i]` is the maximum (count subarrays by
  their maximum: $(i - L)(R - i)$, ties broken by index on one
  side) and the largest rectangle in a histogram.]
#lang-code(
  cpp: "vi nxt(n, n), st; // nxt[i]: first j > i, a[j] > a[i]
rep(i, 0, n) {
  while (sz(st) && a[st.back()] < a[i])
    nxt[st.back()] = i, st.pop_back();
  st.push_back(i);
}",
  py: "nxt, st = [n] * n, []  # nxt[i]: first j > i, a[j] > a[i]
for i, x in enumerate(a):
  while st and a[st[-1]] < x:
    nxt[st.pop()] = i
  st.append(i)",
)
#fact[Next-greater forest][Link each `i` to `nxt[i]`. The prefix
  maxima of `a[s, n)` are the chain `s, nxt[s], nxt[nxt[s]], ...`,
  so binary lifting on `nxt` counts the records of `a[s, r)`,
  finds the first record above a value, or folds anything over
  them (sums, the longest gap) in $O(log n)$. Prefix minima use
  "next smaller"; cyclic arrays use the doubled array; tree paths:
  Records in _Tree techniques_.]
]
#fact[Two smallest][In a permutation, let $a$ and $b$ hold the
  smallest and second smallest values of `[l, r]`. Then $a$ is
  the nearest position left or right of $b$ with a value below
  $p_b$, so at most $2 n$ pairs $(a, b)$ occur. With $p$ and $q$
  the nearest positions outside $[min(a, b), max(a, b)]$ with
  values below $p_b$, exactly $(min(a, b) - p)(q - max(a, b))$
  intervals have this pair.]
#fact[Queue of any monoid][For gcd, matrix products or other
  folds without a deque trick, use two stacks. Push onto the back
  stack, which keeps the fold of all its elements. Pop from the
  front stack, whose entries store the fold from themselves down
  to its bottom; when it is empty, move the back stack over,
  recomputing those folds. The window fold is
  `op(front top fold, back fold)`, amortized $O(1)$ per
  operation. With two pointers this finds, for every `l`, the
  largest `r` whose fold still satisfies a monotone condition.]
#full[
#fact[Sliding window minimum][A deque of indices with increasing
  values: pop the back while it is not smaller than the new value,
  pop the front once it leaves the window. $O(n)$ for all windows;
  the same deque speeds up a DP of the form
  $"dp"[i] = c_i + min_(i - k <= j < i) "dp"[j]$.]
#lang-code(
  cpp: "deque<int> dq; vi mins; // window [i - k + 1, i]
rep(i, 0, n) {
  while (sz(dq) && a[dq.back()] >= a[i]) dq.pop_back();
  dq.push_back(i);
  if (dq[0] <= i - k) dq.pop_front();
  if (i >= k - 1) mins.push_back(a[dq[0]]);
}",
  py: "from collections import deque
dq, mins = deque(), []  # window [i - k + 1, i]
for i, x in enumerate(a):
  while dq and a[dq[-1]] >= x:
    dq.pop()
  dq.append(i)
  if dq[0] <= i - k:
    dq.popleft()
  if i >= k - 1:
    mins.append(a[dq[0]])",
)
#fact[Aggregates on the stack][An entry may carry an aggregate of
  the indices it absorbed, such as the least DP value under a
  shared range maximum: Block maximum in _DP techniques_.]
#fact[Cursor][A sequence edited at a moving cursor: keep the part
  left of the cursor and the part right of it as two stacks;
  moving the cursor pops one and pushes the other. Any structure
  with "push" and "undo the last push" serves each side (LIS
  tails with an undo log, prefix folds); combine the two at the
  cursor.]
#fact[Contribution][Instead of summing an answer over all
  subarrays or pairs, ask for each element in how many of them it
  is the minimum, the first occurrence, or the crossing pair; the
  monotone stack or a sort gives those counts.]

]
=== Small to large
#full[
#fact[Merging sets][Always insert the smaller container into the
  larger. An element only moves into a container at least twice
  its old one, so it moves $O(log n)$ times: $O(n log n)$ moves,
  times the cost of one insert. Swapping two `std` containers is
  $O(1)$; in Python swap the references. Statistics of every
  subtree: Small to large in _Tree techniques_.]
#lang-code(
  cpp: "if (sz(a) < sz(b)) swap(a, b); // map<int, int> counts
for (auto [x, c] : b) a[x] += c; // sets: a.insert(all(b))",
  py: "if len(a) < len(b):  # a, b: dicts of counts
  a, b = b, a
for x, c in b.items():  # sets: a |= b
  a[x] = a.get(x, 0) + c",
)
]
#fact[Shorter side][For each `s`, the subarrays where `a[s]` is
  the maximum are `l` in `(L_s, s]`, `r` in `[s, R_s)`.
  Enumerating only the shorter of the two sides for every `s` costs
  $O(n log n)$ in total (small to large on the Cartesian tree):
  one counting query per enumerated end.]

=== Union-find tricks
#full[
#fact[Skip pointers][Paint (or delete) all unpainted cells of
  arbitrary ranges, each cell once: a DSU where `find(i)` is the
  first unpainted cell `>= i`; painting `i` joins it to `i + 1`.
  $O(n alpha(n))$ in total. The same gives "next alive
  neighbour" after deletions.]
]
#fact[One DSU per mask][Edges carry masks $w$ of $B$ bits and
  questions use only the edges whose mask contains $m$: keep a DSU
  $D(m)$ for every $m$. Insert edge $(u, v, w)$ into $D(w)$, then
  recurse to $w$ without one bit, for each bit, stopping where $u$
  and $v$ are already joined (they are then joined in every
  submask). A visit joins or stops, and a join spawns at most $B$
  visits: $O(q + B 2^B n)$ finds in total.]

=== Bitsets and their costs
#fact[Word parallelism][One shift-or of $10^6$ bits takes about
  #cpp[0.05 ms (`bitset`)]#py[0.2 ms (a big int)].#full[ DAG
  reachability: OR the successors' sets in reverse topological
  order, $O(n m / 64)$. Typical uses: subset sums (_DP
  techniques_), common neighbours and triangles ($O(n^3 / 64)$
  with row ANDs), Gaussian elimination mod 2, shift-and string
  matching.] To rebuild a subset sum, record for each sum the item
  that first made it reachable (the bits of `nb & ~b`); walking
  back from the target uses each item once.]
#fact[Too many bits][If $n$ bitsets of $N$ bits do not fit
  (reachability to $10^5$ targets), take the targets in chunks of
  $64 k$ bits with one pass over the graph per chunk: the same
  time, $O(n k)$ words of memory. Bitsets cannot delete: give a
  changed set a fresh id and AND with a mask of the ids valid at
  that moment.]
#cpp[#fact[C++][`bitset<N>` needs a compile-time `N`: use the
  maximum, or a `template<int N>` that recurses to the next power
  of two. `_Find_first()` and `_Find_next(i)` walk the set bits
  in $O(N / 64 + "count")$.]]
#py[#fact[PyPy][A big int is the bitset: `b | b << w`, `&`, `^` run
  about 4 times slower than C++ bitsets. Masks keep the size
  bounded: `b = (b | b << w) & M`. Popcount (no `bit_count`
  before 3.10) is `bin(b).count("1")`, 100 times slower than an
  AND. Reading one bit with `b >> i & 1` copies the whole number,
  $O(N / 64)$: to read many bits convert once with `bin(b)`, or
  keep a list of ints of 60 bits instead.]]
