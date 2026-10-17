/**
 * Title: Euler path
 * Description: Hierholzer's algorithm: a walk from `src` that uses
 *   every edge exactly once. Returns the vertices `vs` ($m + 1$ of
 *   them, ending where the walk ends) and the edge ids `es`, where
 *   edge `es[i]` goes from `vs[i]` to `vs[i + 1]`. Returns empty
 *   vectors if no such walk exists.
 * Precondition: Edges `(u, v)` with ids `0..m-1`, directed `u -> v`
 *   or undirected; vertices 0-indexed. Self-loops and parallel
 *   edges are fine; isolated vertices are ignored.
 * Time: $O(n + m)$, iterative.
 * Usage:
 *   auto [vs, es] = euler_path(n, edges, directed, src);
 *   if (vs.empty()) ... // impossible from src
 * Tweak: Circuit only: also require `vs.back() == src`. Orienting
 *   so that every $|"in" - "out"| <= 1$: pair the odd vertices by
 *   dummy edges, then orient along a circuit of each component.
 * Example:
 *   // De Bruijn sequence: every length-k word over {0..a-1} once.
 *   // Vertices: (k-1)-words; word w*a + c is the edge from w to
 *   // (w*a + c) % a^(k-1). Directed Euler circuit from 0.
 *   int nv = 1; rep(i, 0, k - 1) nv *= a;
 *   vector<pii> edges;
 *   rep(w, 0, nv) rep(c, 0, a)
 *     edges.push_back({w, (w * a + c) % nv});
 *   auto [vs, es] = euler_path(nv, edges, true, 0);
 *   string s(k - 1, '0'); // then append each edge's last letter
 *   for (int e : es) s += char('0' + e % a);
 * Notes: Undirected: a path exists iff the edges are connected and
 *   0 or 2 vertices have odd degree (a circuit iff 0); it must start
 *   at an odd vertex if there are two. Directed: the edges are
 *   connected (ignoring directions) and every vertex has
 *   $"in" = "out"$ (circuit), or the start has $"out" = "in" + 1$,
 *   the end $"in" = "out" + 1$ and all others balance. Think of it
 *   for "use every edge/domino/word once", "reconstruct a string
 *   from its overlapping pieces" (pieces are edges), De Bruijn
 *   sequences, and "draw without lifting the pen". Covering all
 *   edges of an undirected component with $2k$ odd vertices needs
 *   $max(1, k)$ paths: pair odd vertices with dummy edges, run a
 *   circuit, cut at the dummies. For the orientation in Tweak
 *   make one call per component: a walk never leaves its
 *   component.
 * Source: kactl EulerWalk (CC0), with explicit degree checks.
 * Status: stress-tested against exhaustive search
 */
#pragma once

pair<vi, vi> euler_path(int n, const vector<pii>& edges,
                        bool directed, int src) {
  int m = sz(edges), bad = 0;
  vector<vector<pii>> adj(n);
  vi deg(n), it(n), used(m), vs, es;
  rep(e, 0, m) {
    auto [u, v] = edges[e];
    adj[u].push_back({v, e}), deg[u]++;
    if (directed) deg[v]--;
    else adj[v].push_back({u, e}), deg[v]++;
  }
  // Pretend the walk enters src once more: then at most one
  // vertex (the end) may be unbalanced.
  deg[src] += directed ? -1 : 1;
  for (int d : deg) bad += directed ? d != 0 : d & 1;
  if (bad > 1) return {};
  vector<pii> st = {{src, -1}};
  while (sz(st)) {
    auto [u, e] = st.back();
    if (it[u] == sz(adj[u])) {
      vs.push_back(u), es.push_back(e), st.pop_back();
      continue;
    }
    auto [v, f] = adj[u][it[u]++];
    if (!used[f]) used[f] = 1, st.push_back({v, f});
  }
  if (sz(vs) != m + 1) return {}; // edges not connected to src
  es.pop_back();
  reverse(all(vs)), reverse(all(es));
  return {vs, es};
}
