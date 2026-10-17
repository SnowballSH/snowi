/**
 * Title: Edge colouring
 * Description: `edge_colouring`: a bipartite multigraph with
 *   $Delta$ colours, the maximum degree (Kőnig). `vizing`: a simple
 *   graph with at most $Delta + 1$ colours (Misra-Gries).
 * Precondition: `edge_colouring`: edges `(u, v)` with `u` in
 *   `[0, L)` on the left and `v` in `[0, R)` on the right; parallel
 *   edges are fine. `vizing`: vertices `0..n-1`, no loops or
 *   parallel edges.
 * Time: $O(m (n + Delta))$ time and $4 n Delta$ bytes, $n$ the
 *   number of vertices.
 * Usage:
 *   vi col = edge_colouring(L, R, edges); // in [0, Delta)
 *   vi col2 = vizing(n, edges); // in [0, Delta]
 * Example:
 *   // Games between a_i of league A and b_i of league B, each
 *   // team in at most one game per round: Delta rounds, game i
 *   // in round col[i].
 *   // An n x n matrix of nonnegative integers with all row and
 *   // column sums k, as k permutation matrices: one edge (i, j)
 *   // per unit of a[i][j]; colour c is a permutation.
 * Notes: For each edge, take a colour `a` free at `u` and `b` free
 *   at `v`. If they differ, flip `a` and `b` along the alternating
 *   path from `v` that starts with colour `a`; in a bipartite graph
 *   it cannot reach `u`, and afterwards `a` is free at both ends.
 *   `vizing` recolours a fan of neighbours of `u` and one
 *   alternating path per edge. Deciding whether a general graph
 *   needs $Delta$ or $Delta + 1$ colours is NP-hard.
 * Source: Kőnig's proof by alternating paths; kactl EdgeColoring
 *   (CC0).
 * Status: stress-tested against brute force
 */
#pragma once

vi edge_colouring(int L, int R, const vector<pii>& edges) {
  int D = 0;
  vi deg(L + R), col(sz(edges)), ends; // XOR of the two ends
  for (auto [u, v] : edges) {
    D = max({D, ++deg[u], ++deg[L + v]});
    ends.push_back(u ^ (L + v));
  }
  vector<vi> at(L + R, vi(D, -1)); // at[x][c]: edge coloured c
  rep(e, 0, sz(edges)) {
    int u = edges[e].first, v = L + edges[e].second;
    int a = 0, b = 0;
    while (at[u][a] >= 0) a++;
    while (at[v][b] >= 0) b++;
    // Flip a and b on the path from v, so that a is free at v.
    for (int x = v, c = a; a != b && x >= 0; c ^= a ^ b) {
      int f = at[x][c];
      swap(at[x][a], at[x][b]);
      x = f < 0 ? -1 : ends[f] ^ x;
    }
    at[u][a] = at[v][a] = e;
  }
  rep(x, 0, L)
    rep(c, 0, D)
      if (at[x][c] >= 0) col[at[x][c]] = c;
  return col;
}
vi vizing(int n, const vector<pii>& edges) {
  vi cc(n + 1), res(sz(edges)), fan(n), fr(n), loc;
  for (auto [a, b] : edges) ++cc[a], ++cc[b];
  int u, v, ncol = *max_element(all(cc)) + 1;
  vector<vi> adj(n, vi(ncol, -1)); // adj[x][c]: neighbour via c
  for (auto e : edges) {
    tie(u, v) = e;
    fan[0] = v;
    loc.assign(ncol, 0);
    int at = u, end = u, d, c = fr[u], ind = 0, i = 0;
    while (d = fr[v], !loc[d] && (v = adj[u][d]) != -1)
      loc[d] = ++ind, cc[ind] = d, fan[ind] = v;
    cc[loc[d]] = c;
    for (int cd = d; at != -1; cd ^= c ^ d, at = adj[at][cd])
      swap(adj[at][cd], adj[end = at][cd ^ c ^ d]);
    while (adj[fan[i]][d] != -1) {
      int l = fan[i], r = fan[++i], f = cc[i];
      adj[u][f] = l, adj[l][f] = u, adj[r][f] = -1, fr[r] = f;
    }
    adj[u][d] = fan[i], adj[fan[i]][d] = u;
    for (int y : {fan[0], u, end})
      for (int& z = fr[y] = 0; adj[y][z] != -1;) z++;
  }
  rep(i, 0, sz(edges))
    for (tie(u, v) = edges[i]; adj[u][res[i]] != v;) ++res[i];
  return res;
}
