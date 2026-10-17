/**
 * Title: Tree isomorphism
 * Description: Canonical ids (AHU), equal iff isomorphic.
 *   `rooted(adj, r)`: the id of every subtree, rooted at `r`.
 *   `unrooted(adj)`: the smallest id over the one or two
 *   `centers(adj)`, the middle vertices of every longest path.
 *   `par` holds the parents from the last `rooted` call. Ids are
 *   1, 2, ... in order of first appearance.
 * Precondition: `adj` holds both directions of the $n - 1$ edges of
 *   a tree on vertices `0..n-1`, $n >= 1$. Compare only ids from the
 *   same `TreeHash` object.
 * Time: $O(n log n)$ per call; memory $O(n)$ per tree hashed.
 * Usage:
 *   TreeHash th; // one object for all trees compared
 *   vi id = th.rooted(adj, r); // id[v]: subtree of v
 *   bool iso = th.unrooted(adj1) == th.unrooted(adj2);
 * Tweak: Coloured vertices: put `col[q[i]]` at the front of `key`
 *   after sorting. Edge labels: push `{label, child id}` pairs
 *   (`map<vector<pii>, int>`). Forests: hash each tree, then sort
 *   the ids.
 * Example:
 *   // Distinct subtrees of a rooted tree, up to isomorphism:
 *   vi id = th.rooted(adj, r);
 *   int distinct = sz(set<int>(all(id)));
 *   // Automorphisms fixing r, mod a prime: siblings with equal
 *   // ids can be permuted, k equal ones in k! ways.
 *   map<pii, int> cnt; // (parent, id) -> such siblings so far
 *   ll aut = 1;
 *   rep(v, 0, n) if (v != r)
 *     aut = aut * ++cnt[{th.par[v], id[v]}] % MOD;
 *   // Unrooted: r = a center; with two centers c, d, double it
 *   // if th.rooted(adj, c)[c] == th.rooted(adj, d)[d].
 * Notes: Cues: "are these trees the same up to relabelling", "number
 *   of different subtrees", "which roots give the same rooted tree"
 *   (group roots by `rooted(adj, r)[r]` in $O(n^2 log n)$, or by
 *   rerooted hashes below in $O(n)$), symmetric trees, counting
 *   automorphisms (Example) or distinct labellings ($n!$ over the
 *   automorphisms). An isomorphism maps centers to centers, hence
 *   the rooting at a center; two centers are the middle edge of
 *   every diameter. Random hashes instead of exact ids: $h(v) =
 *   product_c (x_k + h(c))$ mod a large prime $p$, with $k$ the
 *   height of `v` and a random $x_k$ per height; two different trees
 *   collide with probability about $n \/ p$ (Schwartz-Zippel). The
 *   sum form $h(v) = 1 + sum_c f(h(c))$, with $f$ a random 64-bit
 *   mixer (splitmix64 of $h$ plus a random constant), is a `merge`
 *   (+), `lift` ($f$) and `fin` (+1), so `tree/rerooting` gives the
 *   hash for every root in $O(n)$. Never use a plain sum or XOR of
 *   child hashes: different trees collide easily.
 * Source: Aho, Hopcroft, Ullman; centers by peeling leaves.
 * Status: stress-tested against isomorphism by permutations
 */
#pragma once

struct TreeHash {
  map<vi, int> ids;
  vi par; // parents in the last rooted() call
  vi rooted(const vector<vi>& adj, int r) {
    int n = sz(adj);
    vi q{r}, lo(n + 1, n), id(n);
    par.assign(n, -1);
    rep(i, 0, n) { // BFS; kids of q[i]: q[lo[i]..lo[i+1])
      lo[i] = sz(q);
      for (int v : adj[q[i]])
        if (v != par[q[i]]) par[v] = q[i], q.push_back(v);
    }
    for (int i = n; i--;) {
      vi key;
      rep(j, lo[i], lo[i + 1]) key.push_back(id[q[j]]);
      sort(all(key));
      int& x = ids[key];
      if (!x) x = sz(ids);
      id[q[i]] = x;
    }
    return id;
  }
  vi centers(const vector<vi>& adj) {
    int n = sz(adj), left = n;
    vi deg(n), q;
    rep(v, 0, n)
      if ((deg[v] = sz(adj[v])) <= 1) q.push_back(v);
    while (left > 2) { // peel all leaves
      left -= sz(q);
      vi nq;
      for (int v : q)
        for (int u : adj[v])
          if (--deg[u] == 1) nq.push_back(u);
      q = nq;
    }
    return q;
  }
  int unrooted(const vector<vi>& adj) {
    int best = INT_MAX;
    for (int c : centers(adj))
      best = min(best, rooted(adj, c)[c]);
    return best;
  }
};
