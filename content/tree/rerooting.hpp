/**
 * Title: Rerooting
 * Description: A tree DP for every root at once. Rooted anywhere,
 *   the value of the subtree of `v` is `fin(acc, v)`, where `acc`
 *   merges `lift(x, c, v)` over the children `c` of `v` (`x` is the
 *   value of the subtree of `c`), starting from `e`. `dp[v]` is the
 *   value of the whole tree rooted at `v`. The values on directed
 *   edges come too: with `par` from root 0 and `p = par[v]`, `dn[v]`
 *   is `lift(x, v, p)` for `v`'s side of the edge and `up[v]` is
 *   `lift(y, p, v)` for `p`'s side (`e` at vertex 0).
 * Precondition: `adj` holds both directions of the $n - 1$ edges of
 *   a tree on vertices `0..n-1`, $n >= 1$. `merge` is associative
 *   and commutative with identity `e`; it needs no inverse.
 * Time: $O(n)$ calls of `merge`, `lift` and `fin`; iterative.
 * Usage:
 *   // eccentricity: the largest distance from each vertex
 *   Reroot<int> rr(adj, 0, [](int a, int b) { return max(a, b); },
 *       [](int x, int, int) { return x + 1; },
 *       [](int acc, int) { return acc; });
 *   rr.dp[v]; // also rr.par, rr.dn, rr.up
 * Tweak: `lift(x, u, v)` carries `u`'s side across the edge to `v`
 *   and may depend on the direction (edges that must be reversed).
 *   Edge weights: let `adj` hold `{v, w}` pairs, keep `pw[v] = w`
 *   (the edge to `par[v]`) in the BFS, and pass `pw[q[i]]` and
 *   `pw[c]` to the two `lift` calls.
 * Example:
 *   // Sum of distances from every vertex; T = {size, sum}.
 *   typedef pair<ll, ll> P;
 *   auto merge = [](P a, P b) {
 *     return P{a.first + b.first, a.second + b.second};
 *   };
 *   auto lift = [](P x, int, int) { // all one edge farther
 *     return P{x.first, x.second + x.first};
 *   };
 *   auto fin = [](P x, int) { return P{x.first + 1, x.second}; };
 *   Reroot<P> rr(adj, P{0, 0}, merge, lift, fin);
 *   // rr.dp[v].second; rr.dn[v].first = size of v's subtree
 *   // Connected vertex sets containing v, mod a prime: product
 *   // over neighbours of (1 + count on their side); zeros are
 *   // fine as nothing is divided.
 *   Reroot<ll> cs(adj, 1LL,
 *       [](ll a, ll b) { return a * b % MOD; },
 *       [](ll x, int, int) { return (x + 1) % MOD; },
 *       [](ll acc, int) { return acc; });
 * Notes: Cues: "for every vertex as the root", "sum of distances (or
 *   farthest vertex) from each vertex", "the best vertex to root the
 *   tree at", answers on both sides of every edge (`dn`, `up`: a
 *   value per directed edge). Write the rooted DP as `fin` of merged
 *   `lift`s first and check it against a single root; rerooting then
 *   needs no new ideas. When `merge` has an inverse (sums), the
 *   direct formula is shorter: sum of distances is `ans[c] = ans[p]
 *   + n - 2 size[c]`. Each vertex costs $O(deg)$ merges, so a heavy
 *   `T` (an array of size $s$) costs $O(n s)$; keep it small or
 *   shift indices instead of copying. Values like $n^2$ distances
 *   need `ll`.
 * Source: Folklore; prefix and suffix merges as in maspy's library
 *   (CC0).
 * Status: stress-tested against a DP from every root
 */
#pragma once

template<class T> struct Reroot {
  vi par;
  vector<T> dp, dn, up;
  Reroot(const vector<vi>& adj, T e, auto merge, auto lift,
         auto fin)
      : par(sz(adj), -1), dp(sz(adj), e), dn(dp), up(dp) {
    int n = sz(adj);
    vi q{0}, lo(n + 1, n); // kids of q[i]: q[lo[i]..lo[i+1])
    rep(i, 0, n) {
      lo[i] = sz(q);
      for (int v : adj[q[i]])
        if (v != par[q[i]]) par[v] = q[i], q.push_back(v);
    }
    for (int i = n; --i;) {
      T acc = e;
      rep(j, lo[i], lo[i + 1]) acc = merge(acc, dn[q[j]]);
      dn[q[i]] = lift(fin(acc, q[i]), q[i], par[q[i]]);
    }
    rep(i, 0, n) {
      int v = q[i];
      T acc = up[v]; // prefix merges, then suffix merges
      rep(j, lo[i], lo[i + 1]) {
        up[q[j]] = acc;
        acc = merge(acc, dn[q[j]]);
      }
      dp[v] = fin(acc, v), acc = e;
      for (int j = lo[i + 1]; j-- > lo[i];) {
        int c = q[j];
        up[c] = lift(fin(merge(up[c], acc), v), v, c);
        acc = merge(dn[c], acc);
      }
    }
  }
};
