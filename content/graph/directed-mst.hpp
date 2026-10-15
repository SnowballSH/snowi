/**
 * Title: Directed MST
 * Description: Minimum arborescence rooted at `root`
 *   (Chu-Liu/Edmonds): every vertex except `root` takes one incoming
 *   edge `par[v] -> v` and all vertices are reached from `root`, at
 *   the smallest total weight. Returns the cost and `par`
 *   (`par[root] = -1`), or `{-1, {}}` if some vertex cannot be
 *   reached from `root`.
 * Precondition: Vertices are 0-indexed, edges `a -> b`. Weights
 *   need $n dot max|w| <= 4 dot 10^18$ (costs and reduced weights
 *   must fit in `ll`). With negative weights a cost of -1 is
 *   possible, so test `par.empty()`. Self-loops and parallel edges
 *   are fine.
 * Time: $O(m log m)$.
 * Usage:
 *   DirectedMST g(n);
 *   g.add_edge(a, b, w); // a -> b
 *   auto [cost, par] = g.solve(root);
 * Tweak: Maximum arborescence: negate the weights. Free root: add
 *   a vertex `n` with edges `n -> v` of weight $B = 1 + sum |w|$
 *   and solve `DirectedMST(n + 1)` from `n`. One root suffices iff
 *   exactly one `v` has `par[v] == n`; that `v` is the best root
 *   and the cost is `cost - B`. Each extra root adds `B` to the
 *   cost, so the Precondition becomes $(n + 1) B <= 4 dot 10^18$.
 * Example:
 *   // Each house builds a well for c[v] or gets a one-way pipe
 *   // u -> v for w. Vertex n is a super root.
 *   DirectedMST g(n + 1);
 *   rep(v, 0, n) g.add_edge(n, v, c[v]);
 *   for (auto [u, v, w] : pipes) g.add_edge(u, v, w);
 *   auto [cost, par] = g.solve(n); // par[v] == n: a well at v
 * Notes: Think of it when every vertex needs exactly one supplier
 *   along directed edges and all must be reachable from a source:
 *   "broadcast from 1", "each node copies from one other node",
 *   "a one-way pipe from a neighbour". An undirected MST is wrong
 *   for one-way edges. "Build your own well for `c[v]` or take a
 *   pipe": a super root with edges of weight `c[v]` (several roots
 *   allowed, no `B` needed). Idea: each vertex picks its cheapest
 *   incoming edge; a cycle among the picks is contracted, and an
 *   edge entering it is discounted by the cycle edge it replaces
 *   (the lazy add on the heap). The number of arborescences comes
 *   from the directed matrix-tree theorem instead.
 * Source: kactl DirectedMST (CC0), with index-based skew heaps and
 *   an iterative merge.
 * Status: stress-tested against brute force over all parent choices
 */
#pragma once
#include "../ds/dsu-rollback.hpp"

struct DirectedMST {
  int n;
  vi from, to, l, r;
  vector<ll> w, key, lz;
  DirectedMST(int n_) : n(n_) {}
  void add_edge(int a, int b, ll c) {
    from.push_back(a), to.push_back(b), w.push_back(c);
  }
  void prop(int x) {
    key[x] += lz[x];
    if (l[x] >= 0) lz[l[x]] += lz[x];
    if (r[x] >= 0) lz[r[x]] += lz[x];
    lz[x] = 0;
  }
  int merge(int a, int b) { // skew heaps; nodes are edge ids
    int res = -1, *hole = &res; // top-down: no deep recursion
    while (a >= 0 && b >= 0) {
      prop(a), prop(b);
      if (key[a] > key[b]) swap(a, b);
      *hole = a, swap(l[a], r[a]), hole = &l[a], a = *hole;
    }
    *hole = max(a, b);
    return res;
  }
  pair<ll, vi> solve(int root) {
    int m = sz(w);
    key = w, lz.assign(m, 0), l.assign(m, -1), r = l;
    RollbackDSU uf(n);
    vi heap(n, -1), seen(n, -1), path(n), q(n), in(n, -1);
    rep(e, 0, m) heap[to[e]] = merge(heap[to[e]], e);
    vector<tuple<int, int, vi>> cycles;
    ll cost = 0;
    seen[root] = root;
    rep(s, 0, n) {
      int u = s, k = 0, x;
      while (seen[u] < 0) {
        int e = heap[u];
        if (e < 0) return {-1, {}};
        prop(e), cost += key[e], lz[e] = -key[e], prop(e);
        heap[u] = merge(l[e], r[e]);
        q[k] = e, path[k++] = u, seen[u] = s;
        u = uf.find(from[e]);
        if (seen[u] == s) { // a cycle: contract it
          int cyc = -1, end = k, t = uf.time();
          do cyc = merge(cyc, heap[x = path[--k]]);
          while (uf.join(u, x));
          u = uf.find(u), heap[u] = cyc, seen[u] = -1;
          cycles.push_back(
              {u, t, vi(q.begin() + k, q.begin() + end)});
        }
      }
      rep(i, 0, k) in[uf.find(to[q[i]])] = q[i];
    }
    for (auto& [u, t, cyc] : cycles | views::reverse) {
      uf.rollback(t);
      int e = in[u];
      for (int c : cyc) in[uf.find(to[c])] = c;
      in[uf.find(to[e])] = e;
    }
    vi par(n, -1);
    rep(v, 0, n)
      if (in[v] >= 0) par[v] = from[in[v]];
    return {cost, par};
  }
};
