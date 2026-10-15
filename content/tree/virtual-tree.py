"""
Title: Virtual tree
Description: Compresses a tree to a set of key vertices plus the
  LCAs of all pairs of them (at most $2k - 1$ vertices for `k`
  keys), keeping ancestry. Returns `vs` sorted by DFS preorder
  `tin`, and `par[i]`, the index in `vs` of the parent of `vs[i]`
  in the compressed tree (-1 for `vs[0]`, the root). The edge
  into `vs[i]` stands for a path of
  `depth[vs[i]] - depth[vs[par[i]]]` tree edges. Only the LCAs of
  keys adjacent in preorder need to be added.
Precondition: `tree` is a tree/lca `LCA` (an `RmqLCA` from tree/lca-rmq,
  full packet, also works) of the tree. Keys may repeat; no keys gives empty
  lists.
Time: $O(k (log k + log n))$ for `k` keys; $O(k log k)$ with
  `RmqLCA`. PyPy, $sum k = 2 dot 10^5$ keys in total: about 0.9 s
  with `LCA`, 0.4 s with `RmqLCA`.
Usage:
  lc = LCA(adj)
  vs, par = virtual_tree(lc, keys)
  for i in range(len(vs) - 1, 0, -1):  # children first
    dp[par[i]] += dp[i]  # dp indexed like vs
Uses: tree/lca
Example:
  >>> adj = [[1, 2], [0, 3, 4], [0, 5], [1], [1, 6], [2], [4]]
  >>> lc = LCA(adj)  # edges 0-1 0-2 1-3 1-4 2-5 4-6
  >>> vs, par = virtual_tree(lc, [6, 3, 5, 3])
  >>> lc.order  # the preorder that sorts vs
  [0, 2, 5, 1, 4, 6, 3]
  >>> vs, par  # adds 0 = lca(5, 6) and 1 = lca(6, 3)
  ([0, 5, 1, 6, 3], [-1, 0, 0, 2, 2])
  >>> # sum of distances over all pairs of distinct keys
  >>> keys = {3, 5, 6}
  >>> cnt = [int(v in keys) for v in vs]  # keys below vs[i]
  >>> total = 0
  >>> for i in range(len(vs) - 1, 0, -1):
  ...   d = lc.depth[vs[i]] - lc.depth[vs[par[i]]]
  ...   total += d * cnt[i] * (len(keys) - cnt[i])
  ...   cnt[par[i]] += cnt[i]
  >>> total  # d(3, 6) + d(3, 5) + d(6, 5) = 3 + 4 + 5
  12
Notes: Use it when queries each name a few special vertices with
  the total size bounded ($sum k_i <= 2 dot 10^5$): run the
  $O(n)$ tree DP on the compressed tree instead, in $O(k)$ per
  query after sorting. Keep per-query arrays indexed by position
  in `vs` (or clear only the touched vertices), never refill a
  list of size `n`. If the DP needs the real root, add it to the
  keys. Non-key vertices of `vs` are branching points: each has
  at least two children. The edges of the smallest subtree that
  connects the keys (Steiner tree) number half the closed tour
  `dist(v_1, v_2) + ... + dist(v_k, v_1)` in preorder; with
  vertex counts `depth + 1`, the union of the root paths of the
  keys has $sum_i ("depth"(v_i) + 1) - sum_i ("depth"(l_i) + 1)$
  vertices for $l_i = "lca"(v_i, v_(i+1))$, $i < k$. Cues:
  "$q$ queries, each with a set of $k_i$ vertices", "cut the
  tree to separate the marked cities" (Codeforces 613D).
Source: kactl CompressTree (CC0).
Status: stress-tested against brute force over all key pairs
"""


def virtual_tree(tree, keys):
  from bisect import bisect_left

  tin, order, lca = tree.tin, tree.order, tree.lca
  ts = sorted(tin[v] for v in keys)  # sort by preorder time
  vs = [order[t] for t in ts]
  ts += [tin[lca(u, v)] for u, v in zip(vs, vs[1:])]
  ts = sorted(set(ts))
  vs = [order[t] for t in ts]
  par = [-1] * len(vs)
  for i in range(1, len(vs)):
    par[i] = bisect_left(ts, tin[lca(vs[i - 1], vs[i])])
  return vs, par
