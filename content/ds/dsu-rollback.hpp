/**
 * Title: Union-find with rollback
 * Description: Union by size without path compression, so every
 *   `join` can be undone. `time()` returns a checkpoint;
 *   `rollback(t)` undoes all joins made after checkpoint `t`.
 * Precondition: Elements are `0..n-1`. Rollbacks are nested: `t`
 *   comes from an earlier `time()` not yet rolled past.
 * Time: $O(log n)$ per `find`/`join`, $O(1)$ amortized per undone
 *   join.
 * Usage:
 *   RollbackDSU d(n);
 *   int t = d.time();
 *   d.join(a, b); d.find(a); d.size(a);
 *   d.rollback(t); // a and b are separate again
 * Tweak: Components: `n - time() / 2` (a join pushes two entries).
 *   Other state (a parity bit, per-set sums) needs its own stack:
 *   push once per `join`, pop back to size `t / 2` in `rollback`.
 * Example:
 *   // Offline dynamic connectivity: an edge alive during [l, r)
 *   // is stored in O(log q) nodes `at[v]` of a segment tree over
 *   // query times; each leaf answers its query.
 *   void dfs(int v, int lo, int hi) {
 *     int t = d.time();
 *     for (auto [a, b] : at[v]) d.join(a, b);
 *     if (hi - lo == 1) answer(lo);
 *     else {
 *       int mid = (lo + hi) / 2;
 *       dfs(2 * v, lo, mid), dfs(2 * v + 1, mid, hi);
 *     }
 *     d.rollback(t);
 *   }
 * Notes: Use it whenever a search or divide and conquer must undo
 *   unions: offline dynamic connectivity (segment tree over time),
 *   Mo's algorithm with rollback, or trying an edge set and then
 *   restoring. No path compression: it would make undo expensive,
 *   and union by size alone keeps depth below $log_2 n$.
 * Source: kactl RollbackUF (CC0).
 * Status: stress-tested against brute force
 */
#pragma once

struct RollbackDSU {
  vi e;
  vector<pii> st; // (index, old value)
  RollbackDSU(int n) : e(n, -1) {}
  int find(int x) { return e[x] < 0 ? x : find(e[x]); }
  int size(int x) { return -e[find(x)]; }
  int time() { return sz(st); }
  void rollback(int t) {
    for (int i = time(); i-- > t;)
      e[st[i].first] = st[i].second;
    st.resize(t);
  }
  bool join(int a, int b) {
    a = find(a), b = find(b);
    if (a == b) return false;
    if (e[a] > e[b]) swap(a, b);
    st.push_back({a, e[a]}), st.push_back({b, e[b]});
    e[a] += e[b], e[b] = a;
    return true;
  }
};
