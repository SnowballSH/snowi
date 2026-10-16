/**
 * Title: XOR basis
 * Description: `max_xor(x)` / `min_xor(x)`: the largest / smallest
 *   `x ^ v` over the span of the inserted numbers; `kth(k)`: the
 *   `k`-th smallest span element. `PrefixXorBasis`: max subset XOR
 *   of `a[l..r]` as `a` is inserted left to right.
 * Precondition: $0 <= x < 2^60$. `kth`: $0 <= k < 2^"rank"$ (`k = 0`
 *   gives 0). `PrefixXorBasis`: insert `a[r]` with `i = r` for
 *   `r = 0, 1, ...`; `max_xor(l)` then covers `a[l..r]`, inclusive.
 * Time: $O(B)$ per operation, $B = 60$.
 * Usage:
 *   XorBasis xb;
 *   bool indep = xb.insert(x);
 *   ll best = xb.max_xor(), low = xb.min_xor(y);
 *   bool in = xb.contains(y);
 *   ll v = xb.kth(k); // xb.rank: basis size
 *   PrefixXorBasis pb;
 *   rep(r, 0, n) pb.insert(a[r], r); // then pb.max_xor(l)
 * Tweak: Set `B` to the bit width: at most 63 with `ll`, and 30 for
 *   `int` input is faster.
 * Example:
 *   // Range max-XOR queries [l, r], offline by r.
 *   sort(all(qs), [](auto& p, auto& q) { return p.r < q.r; });
 *   PrefixXorBasis pb;
 *   int done = 0;
 *   for (auto& q : qs) {
 *     while (done <= q.r) pb.insert(a[done], done), done++;
 *     ans[q.id] = pb.max_xor(q.l);
 *   }
 *   // Online: keep a copy of pb after each insertion.
 * Notes: `n` numbers of rank $r$: $2^r$ distinct subset XORs, and
 *   each one (0 included, via the empty set) is the XOR of exactly
 *   $2^(n - r)$ subsets. Non-empty subsets with XOR 0 exist iff
 *   $n > r$; the smallest non-zero span element is `kth(1)`. The
 *   nonzero slots `p[b]` with `pos[b] >= l` of a `PrefixXorBasis`
 *   form a basis of `a[l..r]` (empty slots have `pos[b] = 0` too),
 *   so their count is its rank. For any other query `insert` them
 *   into an `XorBasis` (copying `p` would leave `rank` at 0). As `l`
 *   decreases the span grows only at the `pos[b]`.
 *   Cues: "maximum XOR of a subset", "can x be formed", "number of
 *   distinct XORs", "XOR of a path where edges may repeat" (insert
 *   `d[u] ^ d[v] ^ w` for non-tree edges, answer
 *   `min_xor(d[s] ^ d[t])`), "XOR of chosen segments" (basis of
 *   prefix XORs). See the linear algebra facts.
 * Source: Folklore linear basis; prefix (timestamp) basis from
 *   Codeforces 1100F.
 * Status: stress-tested against span enumeration
 */
#pragma once

struct XorBasis {
  static const int B = 60;
  ll p[B] = {}; // p[b]: vector with top bit b, or 0
  int rank = 0;
  bool insert(ll x) {
    for (int b; x; x ^= p[b])
      if (!p[b = __lg(x)]) return p[b] = x, ++rank;
    return false;
  }
  ll max_xor(ll x = 0) {
    for (int b = B; b--;) x = max(x, x ^ p[b]);
    return x;
  }
  ll min_xor(ll x) {
    for (int b = B; b--;) x = min(x, x ^ p[b]);
    return x;
  }
  bool contains(ll x) { return !min_xor(x); }
  ll kth(ll k) { // pivot bits of the answer spell out k
    ll x = 0;
    for (int b = B, i = rank; b--;)
      if (p[b] && (x >> b & 1) != (k >> --i & 1)) x ^= p[b];
    return x;
  }
};

struct PrefixXorBasis {
  static const int B = 60;
  ll p[B] = {};
  int pos[B] = {}; // each slot keeps the newest vector
  void insert(ll x, int i) {
    for (int b; x; x ^= p[b]) {
      if (!p[b = __lg(x)] || pos[b] < i)
        swap(p[b], x), swap(pos[b], i);
      if (!x) return; // took an empty slot
    }
  }
  ll max_xor(int l) {
    ll x = 0;
    for (int b = B; b--;)
      if (pos[b] >= l) x = max(x, x ^ p[b]);
    return x;
  }
};
