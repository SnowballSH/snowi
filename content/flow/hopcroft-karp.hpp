/**
 * Title: Hopcroft-Karp matching
 * Description: Maximum matching in a bipartite graph with `n` left
 *   and `m` right vertices. `ml[u]` is the right partner of left
 *   `u` and `mr[v]` the left partner of right `v` (-1 if free).
 *   `cover()` then gives a minimum vertex cover (Kőnig).
 * Precondition: Edges go from left `u` in `[0, n)` to right `v` in
 *   `[0, m)`; both sides are numbered from 0 separately.
 * Time: $O(m' sqrt(n + m))$ for $m'$ edges.
 * Usage:
 *   HopcroftKarp hk(n, m);
 *   hk.add_edge(u, v);
 *   int k = hk.max_matching(); // then hk.ml, hk.mr
 *   auto [cl, cr] = hk.cover(); // after max_matching
 * Example:
 *   // Minimum path cover of a DAG (vertex-disjoint paths):
 *   // match out-copy u to in-copy v for each edge u -> v.
 *   HopcroftKarp hk(n, n);
 *   for (auto [u, v] : edges) hk.add_edge(u, v);
 *   int paths = n - hk.max_matching();
 *   // Paths start at the v with hk.mr[v] < 0 and continue
 *   // v -> hk.ml[v] while that is not -1.
 * Notes: Matching facts (bipartite only): minimum vertex cover =
 *   maximum matching (Kőnig); maximum independent set and minimum
 *   edge cover (no isolated vertices) are both $n + m - k$. Hall:
 *   all left vertices can be matched iff every left set $X$ has
 *   $|N(X)| >= |X|$, and in general
 *   $k = n - max_X (|X| - |N(X)|)$. Cues: "pair up", "assign each
 *   to one", "maximum set with no two conflicting" when conflicts
 *   only join two classes (grid cells by color parity, even and
 *   odd numbers), "fewest rows and columns covering all marks"
 *   (Kőnig). The cover is read off the last search: left vertices
 *   not reached by alternating paths from free left vertices,
 *   plus right vertices that were reached. The search recurses as
 *   deep as the longest augmenting path.
 * Source: kactl hopcroftKarp (CC0), restructured.
 * Status: stress-tested against bitmask DP and Kőnig certificates
 */
#pragma once

struct HopcroftKarp {
  vector<vi> g;
  vi ml, mr, lvl;
  int lim = 0; // layer of the shortest augmenting paths
  HopcroftKarp(int n, int m) : g(n), ml(n, -1), mr(m, -1) {}
  void add_edge(int u, int v) { g[u].push_back(v); }
  bool bfs() {
    lvl.assign(sz(g), -1), lim = INT_MAX;
    vi q;
    rep(u, 0, sz(g))
      if (ml[u] < 0) lvl[u] = 0, q.push_back(u);
    for (int i = 0; i < sz(q) && lvl[q[i]] < lim; i++)
      for (int v : g[q[i]]) {
        int w = mr[v];
        if (w < 0) lim = lvl[q[i]];
        else if (lvl[w] < 0)
          lvl[w] = lvl[q[i]] + 1, q.push_back(w);
      }
    return lim < INT_MAX;
  }
  bool dfs(int u) {
    int d = lvl[u];
    lvl[u] = -1; // visit once per phase
    for (int v : g[u]) {
      int w = mr[v];
      if (w < 0 ? d == lim : lvl[w] == d + 1 && dfs(w))
        return ml[u] = v, mr[v] = u, true;
    }
    return false;
  }
  int max_matching() {
    while (bfs())
      rep(u, 0, sz(g))
        if (ml[u] < 0) dfs(u);
    return sz(g) - (int)count(all(ml), -1);
  }
  pair<vi, vi> cover() {
    vi cl, cr;
    rep(u, 0, sz(g))
      if (lvl[u] < 0) cl.push_back(u);
    rep(v, 0, sz(mr))
      if (mr[v] >= 0 && lvl[mr[v]] >= 0) cr.push_back(v);
    return {cl, cr};
  }
};
