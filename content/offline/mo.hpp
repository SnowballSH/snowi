/**
 * Title: Mo's algorithm
 * Description: Offline subarray queries `[l, r)`: `add(i)` /
 *   `remove(i)` move `a[i]` into or out of the window; `answer(q)`
 *   runs when it equals `qs[q]`.
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
 *   near $n$. Updates: give each query its time `t` (updates
 *   before it), sort by `(l / B, r / B, t)` with $B = n^(2/3)$
 *   and move `t` too by applying or undoing one update. Mex or
 *   `k`-th value of the window: counts per value and per block of
 *   $sqrt(V)$ values, $O(1)$ per move and $O(sqrt(V))$ per answer.
 *   Tree paths: tour with `v` at `tin[v]` and `tout[v]` ($2n$
 *   slots), `add` and `remove` both toggle `v`; for
 *   `tin[u] <= tin[v]` the window is
 *   `[tin[u], tin[v] + 1)` if `u` is an ancestor of `v`, else
 *   `[tout[u], tin[v] + 1)` plus the lca (toggle it in `answer`
 *   and back). Rollback Mo when `remove` is hard (max, DSU): per
 *   block of `l`, sort by `r`, grow `r` from the block end, add
 *   `[l, block end)` per query, answer, undo; brute-force queries
 *   inside one block.
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
 *   prefer value blocks (Tweak). Mo with updates: moving `t` swaps the stored
 *   value at the update's position, and if that position is in
 *   the window removes the old value and adds the new one; the
 *   total is $O(n^(5/3))$ for $n approx q$.
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
