/**
 * Title: Centroid decomposition
 * Description: Repeatedly removes a centroid (a vertex whose
 *   removal leaves parts of at most half the piece) and recurses
 *   into the parts. `par[v]` is the parent of `v` in the centroid
 *   tree (-1 for the first centroid), `level[v]` its depth there,
 *   and `dist[k][v]` the distance in the tree from `v` to its
 *   centroid ancestor at level `k` (for `k <= level[v]`). The
 *   path between `u` and `v` passes through their lowest common
 *   centroid ancestor `c`, so `dist(u, v) = dist(u, c) + dist(c, v)`.
 * Precondition: `adj` holds both directions of the $n - 1$ edges of
 *   a tree on vertices `0..n-1`, $n >= 1$.
 * Time: $O(n log n)$ time and memory, iterative; at most
 *   $floor(log_2 n) + 1$ levels.
 * Usage:
 *   Centroid cd(adj); // nearest marked vertex, O(log n) per call:
 *   vi best(n, INT_MAX / 2); // best[c]: nearest in c's piece
 *   auto mark = [&](int v) {
 *     for (int c = v; c >= 0; c = cd.par[c])
 *       best[c] = min(best[c], cd.dist[cd.level[c]][v]);
 *   };
 *   auto nearest = [&](int v) { // INT_MAX / 2 if none
 *     int r = INT_MAX / 2;
 *     for (int c = v; c >= 0; c = cd.par[c])
 *       r = min(r, best[c] + cd.dist[cd.level[c]][v]);
 *     return r;
 *   };
 * Tweak: Weighted edges: store `ll` distances and add the weight in
 *   the BFS. Sum of distances from `v` to a changing marked set:
 *   per centroid `c` keep `cnt[c]` and `sum[c]`, the number and
 *   distance sum (to `c`) of marked vertices in its piece, and
 *   `sub[c]`, their distance sum to `par[c]`; walking up from `v`
 *   with `ch` the previous centroid, add `sum[c] - sub[ch] +
 *   (cnt[c] - cnt[ch]) * d(v, c)` (without the `ch` terms at `c =
 *   v`).
 * Example:
 *   // Count unordered pairs {u, v} with dist(u, v) == K >= 1.
 *   // cnt[c][d]: vertices added so far in c's piece at distance d
 *   // from c; sub[c][d]: the same for c's piece, measured from
 *   // par[c] (pairs inside one part must not count at par[c]).
 *   vi size(n); // piece sizes
 *   rep(v, 0, n) for (int c = v; c >= 0; c = cd.par[c]) size[c]++;
 *   vector<vi> cnt(n), sub(n);
 *   rep(c, 0, n) cnt[c].resize(size[c]);
 *   rep(c, 0, n) sub[c].resize(size[c] + 1);
 *   ll pairs = 0;
 *   rep(v, 0, n) {
 *     int ch = -1; // the previous centroid on v's way up
 *     for (int c = v; c >= 0; ch = c, c = cd.par[c]) {
 *       int d = cd.dist[cd.level[c]][v], e = K - d;
 *       if (0 <= e && e < size[c]) pairs += cnt[c][e];
 *       if (ch >= 0 && 0 <= e && e <= size[ch]) pairs -= sub[ch][e];
 *       cnt[c][d]++;
 *       if (ch >= 0) sub[ch][d]++;
 *     }
 *   }
 * Notes: Every path is counted at exactly one centroid: the first
 *   one removed on it, the lowest common centroid ancestor of its
 *   ends. Each vertex lies in at most $floor(log_2 n) + 1$ pieces,
 *   so $O(s)$ work per piece of size $s$ costs $O(n log n)$ in
 *   total, and $O(s log s)$ work (sorting depths, an FFT of the
 *   depth counts) $O(n log^2 n)$. For nearest or minimum queries
 *   (Usage), meeting a pair again at a higher common centroid only
 *   overestimates their distance, so no correction is needed; for
 *   counts and sums it is (the Example's `sub`). Cues: "count paths
 *   of length `K` or at most `K`", "number of pairs whose path has
 *   a property", "distance to the nearest red vertex with
 *   recolouring", "sum of distances to a changing set", and DP
 *   transitions $min_v ("val"_v + d(u, v))$ over active `v` (keep
 *   $"val"_v + d(v, c)$ in a multiset per centroid `c`). Size
 *   arrays per centroid by the piece size, never by `n`. Paths of
 *   length at most `K`: sort each piece's depths and count pairs
 *   with two pointers, minus the same inside each part.
 * Source: Classic; Codeforces 342E pattern.
 * Status: stress-tested against brute force
 */
#pragma once

struct Centroid {
  vi par, level;
  vector<vi> dist; // dist[k][v]: v to its level-k ancestor
  Centroid(const vector<vi>& adj)
      : par(sz(adj), -1), level(sz(adj), -1) {
    int n = sz(adj);
    vi bpar(n), size(n), q;
    // BFS over the piece of s: parents, sizes, depths d
    auto bfs = [&](int s, vi& d) {
      q = {s}, bpar[s] = -1, d[s] = 0;
      rep(i, 0, sz(q))
        for (int v : adj[q[i]])
          if (v != bpar[q[i]] && level[v] < 0)
            bpar[v] = q[i], d[v] = d[q[i]] + 1, q.push_back(v);
      for (int v : q) size[v] = 1;
      for (int i = sz(q); --i;) size[bpar[q[i]]] += size[q[i]];
    };
    auto centroid = [&](int s) { // walk towards the big part
      for (int p = -1, tot = size[s]; p != s;) {
        p = s;
        for (int v : adj[p])
          if (v != bpar[p] && level[v] < 0 && 2 * size[v] > tot)
            s = v;
      }
      return s;
    };
    dist.emplace_back(n);
    bfs(0, dist[0]); // the first centroid redoes dist[0]
    vector<pii> st{{centroid(0), -1}};
    while (sz(st)) {
      auto [c, p] = st.back();
      st.pop_back();
      par[c] = p, level[c] = p < 0 ? 0 : level[p] + 1;
      if (sz(dist) == level[c]) dist.emplace_back(n);
      bfs(c, dist[level[c]]);
      for (int v : adj[c])
        if (level[v] < 0) st.push_back({centroid(v), c});
    }
  }
};
