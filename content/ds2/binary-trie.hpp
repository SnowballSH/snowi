/**
 * Title: Binary trie
 * Description: Multiset of integers in $[0, 2^B)$ stored as a trie of
 *   their bits, highest bit first, with a count per node.
 *   `kth(k, x)` is the `k`-th smallest (0-indexed) of the values
 *   $y xor x$ over stored `y`; so `kth(k)` is the `k`-th smallest
 *   value, and `kth(0, x)` and `kth(size() - 1, x)` are the minimum
 *   and maximum xor with `x`. `count_less(x, k)` counts stored `y`
 *   with $(x xor y) < k$.
 * Precondition: $0 <= x, y < 2^B$; `add(y, -1)` needs `y` stored;
 *   `kth` needs `0 <= k < size()`; `count_less` takes any
 *   $k >= 0$.
 * Time: $O(B)$ per operation. Memory: up to $n B$ nodes of 12
 *   bytes (72 MB for $2 dot 10^5$ values with $B = 30$).
 * Usage:
 *   BinaryTrie tr;
 *   tr.add(y, 1); tr.add(y, -1); // insert, erase one copy
 *   ll lo = tr.min_xor(x), hi = tr.max_xor(x);
 *   ll v = tr.kth(k); int c = tr.count_less(x, k);
 *   int copies = tr.count_less(y, 1); // count of y
 * Tweak: Values up to $2^60$: `B = 60`. Near the memory limit,
 *   `reserve` `n * B + 2` nodes in `ch` and `cnt` first: vector
 *   growth briefly needs twice the memory. For many tries (one
 *   per prefix or subtree), share the pool with one root each,
 *   copying the path on `add` for a persistent version.
 * Example:
 *   // Max xor of a subarray: best p[i] ^ p[j] over prefix xors.
 *   BinaryTrie tr;
 *   ll p = 0, best = 0;
 *   tr.add(0, 1);
 *   for (ll a : arr) {
 *     p ^= a;
 *     best = max(best, tr.max_xor(p));
 *     tr.add(p, 1);
 *   }
 *   // Pairs i < j with arr[i] ^ arr[j] < k.
 *   BinaryTrie seen;
 *   ll pairs = 0;
 *   for (ll a : arr) pairs += seen.count_less(a, k), seen.add(a, 1);
 * Notes: Think of it for "maximum xor of a pair", "maximum xor of
 *   a subarray" (prefix xors), "k-th smallest pairwise xor" (binary
 *   search on the answer with `count_less`), and online "minimum
 *   xor with any stored value". A walk follows the bit of `x` to
 *   keep a 0 in the answer and crosses to the other child for a 1;
 *   `kth` crosses when the 0-side holds at most `k` values.
 *   Erasing only decrements counts, so empty nodes stay; a count
 *   of 0 means missing. When a query may only use the `y` from an
 *   index range `[l, r)`, use a persistent trie (one root per
 *   prefix, subtract the counts of two versions) or sort the
 *   queries offline. "Max xor of any subset" is a different
 *   problem: use a linear basis over $FF_2$.
 * Source: Folklore.
 * Status: stress-tested against brute force
 */
#pragma once

struct BinaryTrie {
  static constexpr int B = 30;
  vector<array<int, 2>> ch; // node 0 is null, node 1 the root
  vi cnt;
  BinaryTrie() : ch(2), cnt(2) {}
  void add(ll y, int c) { // c copies; c < 0 erases
    int v = 1;
    cnt[v] += c;
    for (int i = B; i--;) {
      int b = y >> i & 1;
      if (!ch[v][b])
        ch[v][b] = sz(cnt), ch.push_back({}), cnt.push_back(0);
      cnt[v = ch[v][b]] += c;
    }
  }
  int size() { return cnt[1]; }
  ll kth(int k, ll x = 0) { // k-th smallest y ^ x
    ll res = 0;
    for (int i = B, v = 1; i--;) {
      int b = x >> i & 1, c = cnt[ch[v][b]];
      if (c <= k) k -= c, b ^= 1, res |= 1LL << i;
      v = ch[v][b];
    }
    return res;
  }
  ll min_xor(ll x) { return kth(0, x); }
  ll max_xor(ll x) { return kth(size() - 1, x); }
  int count_less(ll x, ll k) { // #y with (x ^ y) < k
    if (k >> B) return size();
    int res = 0;
    for (int i = B, v = 1; i-- && v;) {
      int b = x >> i & 1;
      if (k >> i & 1) res += cnt[ch[v][b]], b ^= 1;
      v = ch[v][b];
    }
    return res;
  }
};
