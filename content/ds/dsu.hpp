/**
 * Title: Union-find (DSU)
 * Description: Disjoint sets with union by size and path
 *   compression. `e[x]` is the parent of `x`, or minus the size of
 *   its set when `x` is a root.
 * Precondition: Elements are `0..n-1`.
 * Time: $O(alpha(n))$ amortized per operation.
 * Usage:
 *   DSU d(n);
 *   d.join(a, b); // false if already in the same set
 *   d.same(a, b); d.find(a); d.size(a);
 * Example:
 *   // Kruskal: edges sorted by weight
 *   ll mst = 0;
 *   for (auto [w, u, v] : edges)
 *     if (d.join(u, v)) mst += w;
 * Notes: Think DSU when edges only get added and the questions
 *   are about connectivity, component sizes or counts: Kruskal,
 *   offline "delete edges" read backwards as additions, merging
 *   equal cells of a grid. For deletions in arbitrary order use
 *   the rollback version with offline divide and conquer over
 *   time. For parity or offset constraints ("x and y differ"),
 *   store the distance to the parent and add it up in `find`
 *   (weighted DSU), or use $2n$ nodes for "same/different".
 *   `find` recurses $O(log n)$ deep at most, since sets are joined
 *   by size.
 * Source: kactl UnionFind (CC0).
 * Status: stress-tested against brute force
 */
#pragma once

struct DSU {
  vi e;
  DSU(int n) : e(n, -1) {}
  int find(int x) { return e[x] < 0 ? x : e[x] = find(e[x]); }
  bool same(int a, int b) { return find(a) == find(b); }
  int size(int x) { return -e[find(x)]; }
  bool join(int a, int b) {
    a = find(a), b = find(b);
    if (a == b) return false;
    if (e[a] > e[b]) swap(a, b);
    e[a] += e[b], e[b] = a;
    return true;
  }
};
