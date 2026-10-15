/**
 * Title: Mo's algorithm
 * Description: Answers offline queries on subarrays `[l, r)` by
 *   sliding one window over the array. `add(i)` and `remove(i)`
 *   put `a[i]` into or take it out of the window, and `answer(q)`
 *   runs when the window equals `qs[q]`. Queries are sorted by
 *   block of `l`, then by `r` (descending in odd blocks), or along
 *   a Hilbert curve when `hil` is set.
 * Precondition: $0 <= l <= r <= n < 2^21$. `add` and `remove` act
 *   on either end of the window. The window grows before it
 *   shrinks, so its left end never passes its right end.
 * Time: $O(n sqrt(q))$ calls of `add` and `remove` with block size
 *   $n / sqrt(q)$, plus $O(q log q)$ to sort.
 * Usage:
 *   vi cnt(V), res(sz(qs)); // qs: vector<pii> of [l, r)
 *   int distinct = 0;
 *   mo(n, qs, [&](int i) { distinct += !cnt[a[i]]++; },
 *      [&](int i) { distinct -= !--cnt[a[i]]; },
 *      [&](int q) { res[q] = distinct; });
 * Tweak: `hil = true` (Hilbert order) is often faster when $q$ is
 *   near $n$. If the state cares which end moved, give `add` and
 *   `remove` the end as a second argument. Updates: give each
 *   query its time `t` (updates before it), sort by
 *   `(l / B, r / B, t)` with $B = n^(2/3)$ and move `t` too by
 *   applying or undoing one update. Tree paths: run on an Euler
 *   tour listing each vertex twice, with toggling `add`/`remove`.
 * Example:
 *   // Sum over values of cnt[v]^2 in each a[l, r):
 *   vector<ll> cnt(V), res(sz(qs));
 *   ll cur = 0;
 *   mo(n, qs, [&](int i) { cur += 2 * cnt[a[i]]++ + 1; },
 *      [&](int i) { cur -= 2 * --cnt[a[i]] + 1; },
 *      [&](int q) { res[q] = cur; });
 * Notes: Think of it when queries are offline, the answer does not
 *   merge from two halves (distinct values, mode, pairs of equal
 *   values, "values that occur exactly k times") and moving the
 *   window by one costs $O(1)$; $n = q = 2 dot 10^5$ is fine in
 *   C++. With an $O(log n)$ move it is $O(n sqrt(q) log n)$:
 *   prefer an $O(1)$ move with an $O(sqrt(n))$ query (counts per
 *   value block). Rollback Mo when `remove` is hard (max, DSU):
 *   take the queries with `l` in one block, sorted by `r`; grow
 *   `r` from the block's end only, and for each query add
 *   `[l, "block end")` on top, answer, then undo those adds;
 *   queries inside one block are brute-forced. Mo on tree paths:
 *   in an Euler tour with `v` at `tin[v]` and `tout[v]`, the path
 *   `u-v` (with `tin[u] <= tin[v]`) is the window
 *   `[tin[u], tin[v] + 1)` if `u` is an ancestor of `v`, else
 *   `[tout[u], tin[v] + 1)` plus the lca on the side; a vertex
 *   seen twice is off the path, so `add` and `remove` both toggle
 *   it. Mo with updates: moving `t` swaps the stored value at the
 *   update's position, and if that position is in the window
 *   removes the old value and adds the new one; the total is
 *   $O(n^(5/3))$ for $n approx q$.
 * Source: Mo Tao; kactl (CC0); Hilbert order by gepardo.
 * Status: stress-tested against brute force; moves bounded
 */
#pragma once

ll hilbert(int x, int y) { // position on a 2^21 Hilbert curve
  ll d = 0;
  for (int s = 1 << 20; s; s /= 2) {
    int rx = (x & s) > 0, ry = (y & s) > 0;
    d += (ll)s * s * ((3 * rx) ^ ry);
    if (!ry) {
      if (rx) x ^= s - 1, y ^= s - 1;
      swap(x, y);
    }
  }
  return d;
}

void mo(int n, const vector<pii>& qs, auto add, auto remove,
        auto answer, bool hil = false) {
  int q = sz(qs), B = max(1, int(n / sqrt(q + 1.0)));
  vector<ll> key(q);
  rep(i, 0, q) {
    auto [l, r] = qs[i];
    int b = l / B;
    if (hil) key[i] = hilbert(l, r);
    else key[i] = (ll)b << 32 | (b & 1 ? n - r : r);
  }
  vi ord(q);
  iota(all(ord), 0);
  sort(all(ord), [&](int i, int j) { return key[i] < key[j]; });
  int L = 0, R = 0;
  for (int i : ord) {
    auto [l, r] = qs[i];
    while (L > l) add(--L);
    while (R < r) add(R++);
    while (L < l) remove(L++);
    while (R > r) remove(--R);
    answer(i);
  }
}
