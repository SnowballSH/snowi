/**
 * Title: Virtual tree
 * Description: Compresses a tree to a set of key vertices plus the
 *   LCAs of all pairs of them (at most $2k - 1$ vertices for `k`
 *   keys), keeping ancestry. Returns `vs` sorted by DFS preorder
 *   `tin`, and `par[i]`, the index in `vs` of the parent of `vs[i]`
 *   in the compressed tree (-1 for `vs[0]`, the root).
 * Precondition: Keys may repeat; no keys gives empty lists.
 * Time: $O(k (log k + log n))$ for `k` keys.
 * Usage:
 *   LCA lc(adj);
 *   auto [vs, par] = virtual_tree(lc, keys);
 *   for (int i = sz(vs) - 1; i > 0; i--) // children first
 *     dp[par[i]] += dp[i]; // dp indexed like vs
 * Example:
 *   // Sum of distances over all pairs of distinct keys.
 *   // is_key: vi of size n, all 0 between queries.
 *   auto [vs, par] = virtual_tree(lc, keys);
 *   int m = sz(vs), k = 0;
 *   vector<ll> cnt(m); // keys in the subtree of vs[i]
 *   for (int v : keys) k += !is_key[v], is_key[v] = 1;
 *   rep(i, 0, m) cnt[i] = is_key[vs[i]];
 *   ll total = 0;
 *   for (int i = m - 1; i > 0; i--) {
 *     ll len = lc.depth[vs[i]] - lc.depth[vs[par[i]]];
 *     total += len * cnt[i] * (k - cnt[i]);
 *     cnt[par[i]] += cnt[i];
 *   }
 *   for (int v : keys) is_key[v] = 0;
 * Notes: Only the LCAs of keys adjacent in preorder need to be
 *   added. An `RmqLCA` (tree/lca-rmq) also works: $O(k log k)$.
 *   Use it when queries each name a few special vertices with
 *   the total size bounded ($sum k_i <= 2 dot 10^5$): run the
 *   $O(n)$ tree DP on the compressed tree instead, in $O(k)$ per
 *   query after sorting. Keep per-query arrays indexed by position
 *   in `vs` (or clear only the touched vertices), never refill an
 *   array of size `n`. If the DP needs the real root, add it to the
 *   keys. Non-key vertices of `vs` are branching points: each has
 *   at least two children. The edges of the smallest subtree that
 *   connects the keys (Steiner tree) number half the closed tour
 *   `dist(v_1, v_2) + ... + dist(v_k, v_1)` in preorder; with
 *   vertex counts `depth + 1`, the union of the root paths of the
 *   keys has $sum_i ("depth"(v_i) + 1) - sum_i ("depth"(l_i) + 1)$
 *   vertices for $l_i = "lca"(v_i, v_(i+1))$, $i < k$. Cues:
 *   "$q$ queries, each with a set of $k_i$ vertices", "cut the
 *   tree to separate the marked cities" (Codeforces 613D).
 * Source: kactl CompressTree (CC0).
 * Status: stress-tested against brute force over all key pairs
 */
#pragma once
#include "lca.hpp"

pair<vi, vi> virtual_tree(auto& tree, vi keys) {
  auto& tin = tree.tin;
  auto cmp = [&](int a, int b) { return tin[a] < tin[b]; };
  sort(all(keys), cmp);
  int k = sz(keys);
  rep(i, 1, k) keys.push_back(tree.lca(keys[i - 1], keys[i]));
  sort(all(keys), cmp);
  keys.erase(unique(all(keys)), keys.end());
  vi par(sz(keys), -1);
  rep(i, 1, sz(keys)) {
    int p = tree.lca(keys[i - 1], keys[i]);
    par[i] = lower_bound(all(keys), p, cmp) - keys.begin();
  }
  return {keys, par};
}
