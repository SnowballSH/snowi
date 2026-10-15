/**
 * Title: Parallel binary search
 * Description: For each query `i`, finds the fewest updates `t`
 *   (from `0..T`) after which `pred(i)` holds, or `T + 1` if it
 *   never does. All the binary searches advance together: a round
 *   resets the state, applies the updates in order and tests each
 *   unfinished query at its own midpoint.
 * Precondition: `pred(i)` is monotone: once true, it stays true as
 *   more updates are applied. `reset()` restores the state before
 *   any update, `apply(t)` applies update `t` (`0 <= t < T`), and
 *   `pred(i)` only reads the state.
 * Time: $ceil(log_2(T + 2))$ rounds, each with one `reset`, $T$
 *   calls of `apply`, one `pred` per unfinished query and
 *   $O(T + Q)$ extra work.
 * Usage:
 *   // ds/fenwick; update t: a[p[t]] += w[t] with w[t] >= 0
 *   Fenwick fw(n);
 *   vi res = parallel_bs(T, Q, [&] { fw = Fenwick(n); },
 *       [&](int t) { fw.add(p[t], w[t]); },
 *       [&](int i) { return fw.query(r[i]) >= need[i]; });
 * Tweak: Over values (k-th smallest of `a[l, r)`): sort the
 *   positions by value, update `t` marks the `t`-th of them, and
 *   `pred(i)` is "at least `k` marks in `[l, r)`"; the answer is
 *   `a` at the `res[i] - 1`-th sorted position. DSU state (first
 *   edge prefix joining `u` and `v`): `reset` makes a fresh DSU,
 *   `apply(t)` joins edge `t`. Faster rounds: stop applying after
 *   the largest midpoint.
 * Example:
 *   // METEORS (POI 2011): sector j of a ring of m belongs to
 *   // state own[j]; shower t adds x[t] to sectors l[t]..r[t]
 *   // (maybe wrapping); when does state i collect need[i]?
 *   vector<vi> sec(Q);
 *   rep(j, 0, m) sec[own[j]].push_back(j);
 *   Fenwick fw(m); // differences: value of j is query(j + 1)
 *   vi res = parallel_bs(T, Q, [&] { fw = Fenwick(m); },
 *       [&](int t) {
 *         fw.add(l[t], x[t]), fw.add(r[t] + 1, -x[t]);
 *         if (l[t] > r[t]) fw.add(0, x[t]); // wraps around
 *       },
 *       [&](int i) {
 *         ll s = 0; // stop early: the full sum can overflow
 *         for (int j : sec[i])
 *           if ((s += fw.query(j + 1)) >= need[i]) break;
 *         return s >= need[i];
 *       });
 *   // res[i] = T + 1: "NIE"; else shower number res[i] (1-based)
 * Notes: Think of it for "for each query, the first day (edge,
 *   update) after which ..." when one check needs the state after
 *   `t` updates and rebuilding it per query is too slow; also for
 *   "the smallest value `v` such that ..." with values sorted into
 *   updates. A round costs one sweep, so the total is
 *   $O((T + Q) log T)$ structure operations: about 18 rounds for
 *   $T = 2 dot 10^5$. Allow both "already true" (0, which is why
 *   `lo` starts at $-1$) and "never" ($T + 1$). Negative updates
 *   break monotonicity. The recursive form ("overall binary
 *   search") splits the queries and the updates at the middle
 *   value together; it needs no `reset`, only an undo or
 *   subtracting what the left half gave. Online alternatives: a
 *   Kruskal reconstruction tree (DSU), a persistent segment tree
 *   or wavelet matrix (k-th smallest).
 * Source: UCup S4-S5 Training Packet 14; POI 2011 Meteors.
 * Status: stress-tested against brute force; rounds bounded
 */
#pragma once

vi parallel_bs(int T, int Q, auto reset, auto apply,
               auto pred) {
  vi lo(Q, -1), hi(Q, T + 1); // answer in (lo, hi]
  while (1) {
    vector<vi> at(T + 1);
    bool more = false;
    rep(i, 0, Q)
      if (hi[i] - lo[i] > 1) {
        at[(lo[i] + hi[i]) / 2].push_back(i);
        more = true;
      }
    if (!more) return hi;
    reset();
    rep(t, 0, T + 1) {
      for (int i : at[t]) (pred(i) ? hi : lo)[i] = t;
      if (t < T) apply(t);
    }
  }
}
