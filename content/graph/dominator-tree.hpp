/**
 * Title: Dominator tree
 * Description: Lengauer-Tarjan. Returns `idom`: the immediate
 *   dominator of each vertex (the closest strict dominator), with
 *   `idom[s] = s` and `idom[v] = -1` if `v` is unreachable from
 *   `s`. The edges `idom[v] -> v` form a tree rooted at `s`; the
 *   dominators of `v` are its ancestors.
 * Precondition: `adj[u]` lists the heads of the edges leaving `u`;
 *   vertices 0-indexed. Self-loops and parallel edges are fine.
 * Time: $O((n + m) log n)$, no recursion.
 * Usage:
 *   vi idom = dominator_tree(adj, s);
 * Tweak: Paths into a sink `t` instead (post-dominators): run it
 *   on the reversed graph from `t`. A DAG needs no code: in
 *   topological order, `idom[v]` is the LCA of the predecessors of
 *   `v` that are reachable from `s` (binary lifting). Edges on
 *   every path from `s` to `w`: subdivide each edge `u -> v` by a
 *   new vertex `x`; the edge is one iff `x` dominates `w`.
 * Example:
 *   // Critical vertices: those on every s -> t path.
 *   vi idom = dominator_tree(adj, s), crit;
 *   if (idom[t] >= 0)
 *     for (int v = t;; v = idom[v]) {
 *       crit.push_back(v);
 *       if (v == s) break;
 *     }
 *   // lost[v]: vertices cut off from s when v fails, i.e. the
 *   // subtree sizes; children come after parents in BFS order.
 *   vector<vi> ch(n);
 *   rep(v, 0, n)
 *     if (idom[v] >= 0 && v != s) ch[idom[v]].push_back(v);
 *   vi order = {s}, lost(n);
 *   rep(i, 0, sz(order))
 *     for (int c : ch[order[i]]) order.push_back(c);
 *   for (int i = sz(order); i--;) {
 *     int v = order[i];
 *     lost[v]++;
 *     if (v != s) lost[idom[v]] += lost[v];
 *   }
 * Notes: Think of it for "which vertices (or edges) must every
 *   route from `s` use", "removing which node disconnects `v` from
 *   the capital", or "for each vertex, how many become unreachable
 *   if it fails": the subtree sizes of the dominator tree. On an
 *   undirected graph use articulation points (`graph/bcc`) instead.
 * Source: Lengauer and Tarjan (1979), the simple $O(m log n)$
 *   version with path compression.
 * Status: stress-tested against removal-based brute force
 */
#pragma once

vi dominator_tree(const vector<vi>& adj, int s) {
  int n = sz(adj);
  vi tin(n, -1), ord = {s}, par = {-1}, it(n), st = {s};
  tin[s] = 0; // DFS numbers; below, everything is in them
  while (sz(st)) {
    int u = st.back();
    if (it[u] == sz(adj[u])) {
      st.pop_back();
      continue;
    }
    int v = adj[u][it[u]++];
    if (tin[v] < 0)
      tin[v] = sz(ord), ord.push_back(v), par.push_back(tin[u]),
      st.push_back(v);
  }
  int k = sz(ord);
  vector<vi> pred(k), bucket(k);
  for (int u : ord)
    for (int v : adj[u])
      if (tin[v] >= 0) pred[tin[v]].push_back(tin[u]);
  vi semi(k), idom(k), anc(k), best(k), path;
  iota(all(semi), 0), anc = best = semi;
  auto eval = [&](int v) { // min semi on v's forest path
    int x = v;
    for (; anc[x] != x; x = anc[x]) path.push_back(x);
    for (; sz(path); path.pop_back()) {
      int y = path.back(), a = anc[y];
      if (semi[best[a]] < semi[best[y]]) best[y] = best[a];
      anc[y] = x;
    }
    return best[v];
  };
  for (int w = k - 1; w > 0; w--) {
    for (int v : pred[w]) semi[w] = min(semi[w], semi[eval(v)]);
    bucket[semi[w]].push_back(w);
    int p = anc[w] = par[w];
    for (int v : bucket[p]) {
      int u = eval(v);
      idom[v] = semi[u] < semi[v] ? u : p;
    }
    bucket[p].clear();
  }
  vi res(n, -1);
  res[s] = s;
  rep(i, 1, k) {
    if (idom[i] != semi[i]) idom[i] = idom[idom[i]];
    res[ord[i]] = ord[idom[i]];
  }
  return res;
}
