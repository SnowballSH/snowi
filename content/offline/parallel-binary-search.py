"""
Title: Parallel binary search
Description: For each query `i`, finds the fewest updates `t`
  (from `0..T`) after which `pred(i)` holds, or `T + 1` if it never
  does. All the binary searches advance together: a round resets
  the state, applies the updates in order and tests each unfinished
  query at its own midpoint.
Precondition: `pred(i)` is monotone: once true, it stays true as
  more updates are applied. `reset()` restores the state before any
  update, `apply(t)` applies update `t` (`0 <= t < T`), and
  `pred(i)` only reads the state.
Time: $ceil(log_2(T + 2))$ rounds, each with one `reset`, $T$
  calls of `apply`, one `pred` per unfinished query and $O(T + Q)$
  extra work. PyPy, with ds/fenwick: the Usage takes about 1 s for
  $T = Q = 10^5$ and 2.2 s for $2 dot 10^5$. The METEORS Example
  (up to 3 adds per shower and one query per sector each round)
  takes 1.5 s for $n = m = k = 10^5$ but about 7 s at the full
  $3 dot 10^5$, where C++ needs 1.5 s.
Usage:
  # ds/fenwick; update t: a[p[t]] += w[t] with w[t] >= 0
  fw = Fenwick(n)
  def reset():
    fw.t = [0] * n
  res = parallel_bs(T, Q, reset, lambda t: fw.add(p[t], w[t]),
                    lambda i: fw.query(r[i]) >= need[i])
Tweak: Over values (k-th smallest of `a[l:r]`): sort the positions
  by value, update `t` marks the `t`-th of them, and `pred(i)` is
  "at least `k` marks in `[l, r)`"; the answer is `a` at the
  `res[i] - 1`-th sorted position. DSU state (first edge prefix
  joining `u` and `v`): `reset` makes a fresh DSU, `apply(t)` joins
  edge `t`. Faster rounds: stop applying after the largest
  midpoint.
Example:
  # Answers: 0 if true before any update, else the number of
  # updates needed, or T + 1 if never. The state is a total here.
  >>> w, need, s = [2, 0, 3], [0, 1, 5, 6], [0]
  >>> def reset(): s[0] = 0
  >>> def apply(t): s[0] += w[t]
  >>> parallel_bs(3, 4, reset, apply, lambda i: s[0] >= need[i])
  [0, 1, 3, 4]

  # METEORS (POI 2011) with ds/fenwick: sector j of a ring of m
  # belongs to state own[j]; shower t adds x[t] to sectors
  # l[t]..r[t] (wrapping if l[t] > r[t]); when does state i
  # collect need[i]?
  sec = [[] for _ in range(Q)]
  for j in range(m):
    sec[own[j]].append(j)
  fw = Fenwick(m)  # differences: value of j is query(j + 1)
  def reset():
    fw.t = [0] * m
  def apply(t):
    fw.add(l[t], x[t])
    fw.add(r[t] + 1, -x[t])
    if l[t] > r[t]:  # wraps: [l[t], m) and [0, r[t]]
      fw.add(0, x[t])
  def pred(i):
    s = 0
    for j in sec[i]:  # stop early: fewer queries
      s += fw.query(j + 1)
      if s >= need[i]:
        break
    return s >= need[i]
  res = parallel_bs(T, Q, reset, apply, pred)
  # res[i] = T + 1: "NIE"; else shower number res[i] (1-based)
Notes: Think of it for "for each query, the first day (edge,
  update) after which ..." when one check needs the state after `t`
  updates and rebuilding it per query is too slow; also for "the
  smallest value `v` such that ..." with values sorted into
  updates. A round costs one sweep, so the total is
  $O((T + Q) log T)$ structure operations: about 18 rounds for
  $T = 2 dot 10^5$. Allow both "already true" (0, which is why `lo`
  starts at $-1$) and "never" ($T + 1$). Negative updates break
  monotonicity. The recursive form ("overall binary search")
  splits the queries and the updates at the middle value together;
  it needs no `reset`, only an undo or subtracting what the left
  half gave. Online alternatives: a Kruskal reconstruction tree
  (DSU), a persistent segment tree or wavelet matrix (k-th
  smallest).
Source: UCup S4-S5 Training Packet 14; POI 2011 Meteors.
Status: stress-tested against brute force; rounds bounded
"""


def parallel_bs(T, Q, reset, apply, pred):
  lo, hi = [-1] * Q, [T + 1] * Q  # answer in (lo, hi]
  nxt = [-1] * Q  # bucket lists by midpoint: head[mid], nxt[i]
  while True:
    head, more = [-1] * (T + 1), False
    for i in range(Q):
      if hi[i] - lo[i] > 1:
        mid = (lo[i] + hi[i]) >> 1
        nxt[i], head[mid], more = head[mid], i, True
    if not more:
      return hi
    reset()
    for t in range(T + 1):
      i = head[t]
      while i >= 0:
        if pred(i):
          hi[i] = t
        else:
          lo[i] = t
        i = nxt[i]
      if t < T:
        apply(t)
