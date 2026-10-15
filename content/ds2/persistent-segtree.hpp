/**
 * Title: Persistent segment tree
 * Description: Point add and range sum over `[0, n)`, where every
 *   update leaves the old tree intact and returns the root of a new
 *   version. Versions share all untouched nodes. Root 0 is the
 *   all-zero array. `kth(a, b, k)` reads version `b` minus version
 *   `a` as counts and returns the smallest `i` whose count over
 *   `[0, i]` exceeds `k`: the `k`-th smallest (0-indexed) value of
 *   a subarray when version `i` counts the compressed values of
 *   `a[0, i)`.
 * Precondition: 0-indexed, half-open ranges, `n >= 1`. For `kth`,
 *   every count of `b` minus `a` is $>= 0$ and `0 <= k <` their
 *   total.
 * Time: $O(log n)$ per operation. Memory: each `add` creates at
 *   most $ceil(log_2 n) + 1$ nodes of 16 bytes.
 * Usage:
 *   PersistentSegTree st(n);
 *   int v1 = st.add(0, i, 5), v2 = st.add(v1, j, 3);
 *   ll s = st.query(v1, l, r); // sum of [l, r) in version v1
 *   int idx = st.kth(v1, v2, k); // k-th unit of v2 - v1
 * Tweak: Point set: `add(v, i, x - query(v, i, i + 1))`. Max or
 *   min: `add` adds `delta` along the path, which is right only
 *   for sums (or xor), so set the new leaf and recompute the new
 *   path nodes from their children bottom-up; `query` then uses
 *   `max` with identity $-infinity$. Starting from an array `a`:
 *   `n` adds into root 0. Near the memory limit call
 *   `t.reserve(nodes)` first.
 * Example:
 *   // k-th smallest (0-indexed) of a[l, r), values compressed.
 *   vector<ll> vals = a;
 *   sort(all(vals)), vals.erase(unique(all(vals)), vals.end());
 *   PersistentSegTree st(sz(vals));
 *   vi root{0}; // root[i]: counts of the values in a[0, i)
 *   for (ll x : a) {
 *     int c = int(lower_bound(all(vals), x) - vals.begin());
 *     root.push_back(st.add(root.back(), c, 1));
 *   }
 *   ll kth = vals[st.kth(root[l], root[r], k)];
 *   // number of x in a[l, r) with x < y:
 *   int c = int(lower_bound(all(vals), y) - vals.begin());
 *   ll less = st.query(root[r], 0, c) - st.query(root[l], 0, c);
 * Notes: Think of it for "k-th smallest in a range", "count values
 *   in `[x, y)` in a subarray" online, "answer for the array as it
 *   was after update `t`", or when a sweep must be queried at any
 *   earlier moment (versions over time or over a sorted key). The
 *   trick is prefix versions over the value axis: version `r`
 *   minus version `l` is the histogram of `a[l, r)`, and one walk
 *   down both trees finds the k-th unit. On a tree, a version per
 *   vertex built from its parent gives path histograms as
 *   $v + u - "lca" - "parent"("lca")$; walk four roots at once.
 *   Memory is the usual limit: $2 dot 10^5$ updates on
 *   $n = 2 dot 10^5$ make about $3.8 dot 10^6$ nodes (61 MB). On a
 *   static array a wavelet matrix answers the same k-th smallest
 *   and count queries online with far less memory; offline,
 *   parallel binary search also works.
 * Source: Folklore ("chairman tree"); path copying.
 * Status: stress-tested against brute force over all versions
 */
#pragma once

struct PersistentSegTree {
  struct Node {
    int l = 0, r = 0;
    ll sum = 0;
  };
  int n;
  vector<Node> t{Node()}; // node 0: empty tree, children 0
  PersistentSegTree(int len) : n(len) {}
  int add(int x, int i, ll delta) { // returns the new root
    int root = sz(t);
    for (int lo = 0, hi = n;;) {
      Node c = t[x];
      c.sum += delta, t.push_back(c);
      if (hi - lo == 1) return root;
      int mid = (lo + hi) / 2, y = sz(t) - 1;
      if (i < mid) t[y].l = y + 1, x = c.l, hi = mid;
      else t[y].r = y + 1, x = c.r, lo = mid;
    }
  }
  ll query(int x, int l, int r) { return query(x, l, r, 0, n); }
  ll query(int x, int l, int r, int lo, int hi) {
    if (!x || r <= lo || hi <= l) return 0;
    if (l <= lo && hi <= r) return t[x].sum;
    int mid = (lo + hi) / 2;
    return query(t[x].l, l, r, lo, mid) +
           query(t[x].r, l, r, mid, hi);
  }
  int kth(int a, int b, ll k) { // k-th unit of b minus a
    int lo = 0, hi = n;
    while (hi - lo > 1) {
      int mid = (lo + hi) / 2;
      ll c = t[t[b].l].sum - t[t[a].l].sum;
      if (k < c) a = t[a].l, b = t[b].l, hi = mid;
      else k -= c, a = t[a].r, b = t[b].r, lo = mid;
    }
    return lo;
  }
};
