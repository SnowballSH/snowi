"""
Title: Mo's algorithm
Description: Offline subarray queries `[l, r)`: `add(i)` /
  `remove(i)` move `a[i]` into or out of the window; `answer(q)`
  runs when it equals `qs[q]`.
Precondition: $0 <= l <= r <= n < 2^21$. `add` and `remove` are
  called at both ends. The window grows before it shrinks,
  so its left end never passes its right end.
Time: $O(n sqrt(q))$ calls of `add` and `remove` with block size
  $n / sqrt(q)$, plus $O(q log q)$ to sort. In PyPy, distinct
  values for $n = q = 10^5$ take about 1 s, for $2 dot 10^5$ about
  3 s (Hilbert order 20% less).
Usage:
  cnt, res, cur = [0] * V, [0] * len(qs), [0]  # qs: [(l, r)]
  def add(i):
    cur[0] += cnt[a[i]] == 0; cnt[a[i]] += 1
  def remove(i):
    cnt[a[i]] -= 1; cur[0] -= cnt[a[i]] == 0
  def answer(q):
    res[q] = cur[0]
  mo(n, qs, add, remove, answer)  # hil=True: Hilbert order
Tweak: Inlining `add` and `remove` into the four `while` loops
  saves about 25% in PyPy. Updates: give each query its time
  `t` (updates before it), sort by `(l // B, r // B, t)` with
  $B = n^(2/3)$ and move `t` too by applying or undoing one
  update. Tree paths: tour with `v` at `tin[v]` and `tout[v]` ($2n$
  slots), `add` and `remove` both toggle `v`; for
  `tin[u] <= tin[v]` the window is `[tin[u], tin[v] + 1)` if `u` is
  an ancestor of `v`, else `[tout[u], tin[v] + 1)` plus the lca
  (toggle it in `answer` and back). Rollback Mo when `remove` is
  hard (max, DSU): per block of `l`, sort by `r`, grow `r` from
  the block end, add `[l, block end)` per query, answer, undo;
  brute-force queries inside one block.
Example:
  >>> a, qs = [1, 2, 1, 3, 2, 1], [(0, 3), (1, 6), (2, 2)]
  >>> cnt, cur, res = [0] * 4, [0], [0] * len(qs)
  >>> def add(i):  # cur = sum over values v of cnt[v]^2
  ...   cur[0] += 2 * cnt[a[i]] + 1; cnt[a[i]] += 1
  >>> def remove(i):
  ...   cnt[a[i]] -= 1; cur[0] -= 2 * cnt[a[i]] + 1
  >>> def answer(q):
  ...   res[q] = cur[0]
  >>> mo(len(a), qs, add, remove, answer); res
  [5, 9, 0]
Notes: Think of it when queries are offline, the answer does not
  merge from two halves (distinct values, mode, pairs of equal
  values, "values that occur exactly k times") and moving the
  window by one costs $O(1)$. In PyPy, $n = q = 10^5$ is
  comfortable and $2 dot 10^5$ is tight. With an $O(log n)$ move
  it is $O(n sqrt(q) log n)$: prefer an $O(1)$ move with an
  $O(sqrt(n))$ query (counts per value block). Mo with updates:
  moving `t` swaps the stored value at the update's position, and
  if that position is in the window removes the old value and
  adds the new one; the total is $O(n^(5/3))$ for $n approx q$.
  Side-dependent answers (inversions): track `[L, R)` yourself;
  `add(i)` is on the left iff `L < R` and `i < L` (an add to an
  empty window gives `[i, i + 1)`), `remove(i)` iff `i == L`.
Source: Mo Tao; kactl (CC0); Hilbert order by gepardo.
Status: stress-tested against brute force; moves bounded
"""


def hilbert(x, y):  # position on a 2^21 Hilbert curve
  d, s = 0, 1 << 20
  while s:
    rx, ry = x & s > 0, y & s > 0
    d += s * s * ((3 * rx) ^ ry)
    if not ry:
      if rx:
        x, y = x ^ (s - 1), y ^ (s - 1)
      x, y = y, x
    s >>= 1
  return d


def mo(n, qs, add, remove, answer, hil=False):
  B = max(1, int(n / (len(qs) + 1) ** 0.5))
  key = []
  for l, r in qs:
    b = l // B
    if hil:
      key.append(hilbert(l, r))
    else:
      key.append(b << 32 | (n - r if b & 1 else r))
  L = R = 0
  for i in sorted(range(len(qs)), key=key.__getitem__):
    l, r = qs[i]
    while L > l:
      L -= 1
      add(L)
    while R < r:
      add(R)
      R += 1
    while L < l:
      remove(L)
      L += 1
    while R > r:
      R -= 1
      remove(R)
    answer(i)
