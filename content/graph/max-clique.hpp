/**
 * Title: Maximum clique
 * Description: Returns the vertices of one maximum clique (sorted).
 *   Branch and bound: a greedy colouring of the candidates bounds
 *   the clique size (a clique uses each colour at most once), and
 *   branches start at the highest colours.
 * Precondition: `adj[u][v] = adj[v][u] = 1` for each edge;
 *   `adj[v][v] = 0`. At most 128 vertices (the bitset size).
 * Time: Exponential. Random graphs with edge probability 0.9:
 *   $n = 128$ in about 0.05 s, $n = 150$ (larger `B`) about 1 s.
 *   Sparser graphs are much faster.
 * Usage:
 *   vector<MaxClique::B> adj(n);
 *   adj[u][v] = adj[v][u] = 1;
 *   vi clique = MaxClique(adj).solve();
 * Tweak: More vertices: enlarge `B`. Maximum independent set: run
 *   on the complement (`adj[u][v] = !g[u][v]` for `u != v`).
 * Example:
 *   // Most mutually compatible items; cmp(i, j) is symmetric.
 *   vector<MaxClique::B> adj(n);
 *   rep(i, 0, n) rep(j, 0, n) adj[i][j] = i != j && cmp(i, j);
 *   vi best = MaxClique(adj).solve();
 *   // Independent set: adj[i][j] = i != j && !conflict(i, j).
 * Notes: Clique, independent set and vertex cover are one problem:
 *   $S$ is independent iff it is a clique of the complement iff $V
 *   without S$ is a vertex cover, so the minimum vertex cover is
 *   $n$ minus the maximum independent set. Check the structure
 *   before going exponential: in a bipartite graph the maximum
 *   independent set is $n$ minus the maximum matching (Kőnig,
 *   `flow/hopcroft-karp`); trees take a DP, interval graphs a
 *   greedy, partial orders Dilworth. Cues: "pairwise compatible" or
 *   "pairwise conflicting", "largest group where everyone knows
 *   everyone", a general graph with $n <= 100$. For $n <= 40$, meet
 *   in the middle ($2^(n\/2)$ subsets per half) also counts the
 *   independent sets. Dense graphs are the slow case, and
 *   complements of sparse inputs are dense.
 * Source: Tomita and Seki (MCQ) with bitsets as in San Segundo's
 *   BBMC; compare kactl MaxClique (CC0).
 * Status: stress-tested against subset enumeration
 */
#pragma once

struct MaxClique {
  typedef bitset<128> B;
  vector<B> adj;
  vi ord, cur, best;
  MaxClique(const vector<B>& g) : ord(sz(g)) {
    iota(all(ord), 0); // relabel by decreasing degree
    sort(all(ord), [&](int a, int b) {
      return g[a].count() > g[b].count();
    });
    for (int u : ord) {
      adj.emplace_back();
      rep(j, 0, sz(g)) adj.back()[j] = g[u][ord[j]];
    }
  }
  void expand(B cand) {
    vi vs, col;
    int k = 0;
    for (B left = cand; left.any(); k++) // colour classes
      for (B q = left; q.any();) {
        int v = (int)q._Find_first();
        q &= ~adj[v], q[v] = 0, left[v] = 0;
        vs.push_back(v), col.push_back(k + 1);
      }
    for (int i = sz(vs); i--;) {
      if (sz(cur) + col[i] <= sz(best)) return;
      int v = vs[i];
      cur.push_back(v);
      B next = cand & adj[v];
      if (next.any()) expand(next);
      else if (sz(cur) > sz(best)) best = cur;
      cur.pop_back(), cand[v] = 0;
    }
  }
  vi solve() {
    B cand;
    rep(i, 0, sz(adj)) cand[i] = 1;
    expand(cand);
    vi res;
    for (int v : best) res.push_back(ord[v]);
    sort(all(res));
    return res;
  }
};
