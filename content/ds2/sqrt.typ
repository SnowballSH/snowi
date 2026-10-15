#import "/typst/notes.typ": *

When no clean $O(log n)$ structure exists, split the work into a
part that is cheap per item and a part that is cheap per group,
and balance the two. Pick the block size by equating the costs,
then check the total: C++ manages a few times $10^8$ simple steps
per second, PyPy a few times fewer.

=== Blocks
#fact[Range update, range query][Cut the array into blocks of
  $B approx sqrt(n)$. An operation touches $O(n / B)$ whole blocks
  through a per-block summary and lazy tag, and $O(B)$ single
  elements in the two partial blocks, after which those blocks
  are rebuilt. Use it when a block summary is easy but a segment
  tree node is not (mode, values `>= x` under adds).]
#fact[Sorted blocks][Range add and "count values `>= x` in
  `[l, r)`": each block keeps a sorted copy and an add tag; a whole
  block answers by binary search for $x - "tag"$ in $O(log B)$. In
  a partial block, add to the affected elements in place; they
  and the other elements each stay in sorted order, so one merge
  rebuilds the sorted copy in $O(B)$. Then
  $B approx sqrt(n log n)$ gives $O(sqrt(n log n))$ per
  operation; a plain re-sort costs $O(B log B)$ and wants
  $B approx sqrt(n)$, for $O(sqrt(n) log n)$.]
#fact[Balancing][If an operation costs $a dot n / B + b dot B$,
  take $B = sqrt(a n / b)$ for $O(sqrt(a b n))$. Constants are
  rough guesses: time a few block sizes around the formula.]
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
  $O(sqrt(m))$ per operation; all triangles in $O(m sqrt(m))$ by
  orienting each edge from lower to higher (degree, id) and
  intersecting out-lists. The same split works for values by
  frequency: at most $n / D$ values occur more than $D$ times.]
#fact[Small and large parameters][Queries with a step, modulus or
  length $d$: for $d <= sqrt(n)$ precompute a table per $d$
  ($O(n sqrt(n))$ in total), for larger $d$ walk the
  $n / d <= sqrt(n)$ positions directly.]

=== Monotone stack and queue
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
#fact[Queue of any monoid][For gcd, matrix products or other
  folds without a deque trick, use two stacks. Push onto the back
  stack, which keeps the fold of all its elements. Pop from the
  front stack, whose entries store the fold from themselves down
  to its bottom; when it is empty, move the back stack over,
  recomputing those folds. The window fold is
  `op(front top fold, back fold)`, amortized $O(1)$ per
  operation. With two pointers this finds, for every `l`, the
  largest `r` whose fold still satisfies a monotone condition.]

=== Small to large
#fact[Merging sets][Always insert the smaller container into the
  larger. An element only moves into a container at least twice
  its old one, so it moves $O(log n)$ times: $O(n log n)$ moves,
  times the cost of one insert. Swapping two `std` containers is
  $O(1)$; in Python swap the references.]
#lang-code(
  cpp: "if (sz(a) < sz(b)) swap(a, b); // map<int, int> counts
for (auto [x, c] : b) a[x] += c; // sets: a.insert(all(b))",
  py: "if len(a) < len(b):  # a, b: dicts of counts
  a, b = b, a
for x, c in b.items():  # sets: a |= b
  a[x] = a.get(x, 0) + c",
)
#fact[Subtree data][DSU on tree: keep the container of the heavy
  (largest) child and add the light subtrees into it, element by
  element: $O(n log n)$ insertions. Iteratively, process vertices
  in reverse BFS order and keep one container id per vertex.
  Typical: distinct colours, the most frequent colour, or the
  count of each depth, for every subtree.]

=== Offline orders
#fact[Sort by threshold][Queries "elements `>= x` in `[l, r)`"
  or "edges of weight `<= w`": sort queries and elements by the
  threshold, insert elements as the threshold passes them
  (Fenwick tree over positions, DSU over edges) and answer each
  query when its turn comes. Turns a 2-D question into 1-D.]
#fact[Distinct values in ranges][Sort queries by `r`; when `a[r]`
  arrives, add 1 at `r` and remove the 1 at its previous
  occurrence. A query `[l, r)` is the Fenwick sum over `[l, r)`.
  The same "last occurrence" trick counts other per-value events.]
#fact[Reverse time][Deletions are hard but insertions easy (DSU,
  reachability): read all operations, then process them backwards
  so deletions become insertions; answers are output in reverse.]
#fact[Contribution][Instead of summing an answer over all
  subarrays or pairs, ask for each element in how many of them it
  is the minimum, the first occurrence, or the crossing pair; the
  monotone stack or a sort gives those counts.]

=== Bitsets and their costs
#fact[Word parallelism][A bitset operation on $N$ bits costs
  $N / 64$ word steps. One shift-or of $10^6$ bits takes about
  0.05 ms in C++ (`bitset`) and 0.2 ms in PyPy (a big int).
  Typical uses: subset sums (`dp/techniques`), reachability
  in a DAG ($O(n m / 64)$: OR the successors' sets in reverse
  topological order), common neighbours and triangles
  ($O(n^3 / 64)$ with row ANDs), Gaussian elimination mod 2,
  shift-and string matching.]
#cpp[#fact[C++][`bitset<N>` needs a compile-time `N`: use the
  maximum, or a `template<int N>` that recurses to the next power
  of two. `count()` is a popcount; `_Find_first()` and
  `_Find_next(i)` walk the set bits in $O(N / 64 + "count")$.
  `vector<bool>` is packed too but has no word-level operations.]]
#py[#fact[PyPy][A big int is the bitset: `b | b << w`, `&`, `^` run
  about 4 times slower than C++ bitsets. Masks keep the size
  bounded: `b = (b | b << w) & M`. Popcount (no `bit_count`
  before 3.10) is `bin(b).count("1")`, 100 times slower than an
  AND. Reading one
  bit with `b >> i & 1` copies the whole number, $O(N / 64)$: to
  read many bits convert once with `bin(b)`, or keep a list of
  ints of 60 bits instead.]]
