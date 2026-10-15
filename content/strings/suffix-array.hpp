/**
 * Title: Suffix array
 * Description: `sa` lists the start positions of the suffixes of
 *   `s` in sorted order and `rank` is its inverse. `lcp[i]` is the
 *   length of the longest common prefix of the suffixes starting at
 *   `sa[i]` and `sa[i + 1]`, for `i < n - 1`, so the LCP of the
 *   suffixes ranked $i < j$ is the minimum of `lcp[i, j)`.
 * Precondition: Values of `s` in $[1, "lim")$: ASCII strings work
 *   as is.
 * Time: $O(n log n)$.
 * Usage:
 *   SuffixArray S(s); // string or vector<int>
 *   // S.sa[r]: start of the suffix of rank r
 *   ll d = (ll)n * (n + 1) / 2; // distinct substrings:
 *   for (int h : S.lcp) d -= h;
 * Tweak: For large values compress them to $[1, k]$ and pass `lim =
 *   k + 1`. Sorted cyclic shifts: build on `s + s` and keep the
 *   starts `< n`.
 * Example:
 *   // LCP of the suffixes at a and b in O(1), sparse table
 *   // over lcp (needs n >= 2):
 *   SparseTable<int> st(S.lcp);
 *   auto lcp = [&](int a, int b) {
 *     if (a == b) return n - a;
 *     auto [i, j] = minmax(S.rank[a], S.rank[b]);
 *     return st.query(i, j);
 *   };
 *   // The suffixes starting with pat are those in [lo, hi):
 *   auto lo = partition_point(all(S.sa), [&](int i) {
 *     return s.compare(i, sz(pat), pat) < 0;
 *   });
 *   auto hi = partition_point(lo, end(S.sa), [&](int i) {
 *     return s.compare(i, sz(pat), pat) == 0;
 *   });
 *   int occurrences = hi - lo;
 * Notes: The LCP of two suffixes is the minimum of `lcp` between
 *   their ranks (Example). The suffixes sharing a prefix of length
 *   `>= k` with a given one form a block of consecutive ranks: find
 *   it by binary search on the sparse table, or process all `k` at
 *   once with a monotonic stack over `lcp` (the LCP array walks the
 *   suffix tree). Distinct substrings: $n(n+1)/2 - sum_i "lcp"_i$;
 *   the suffix of rank `i` brings $n - "sa"_i - "lcp"_(i-1)$ new
 *   ones, which also finds the `k`-th smallest. The longest
 *   repeated substring has length $max_i "lcp"_i$. Several strings:
 *   join them with separators that are distinct and occur nowhere
 *   (`vector<int>`, letters shifted above them); the longest common
 *   substring of two strings is the maximum `lcp[i]` over adjacent
 *   ranks from different strings, and of `k` strings a sliding
 *   window over ranks that covers all `k`. Cues: "number of
 *   distinct substrings", "longest substring occurring twice", many
 *   LCP or comparison queries between substrings. For one pattern
 *   or one comparison, hashing or Z is shorter; for a string that
 *   grows online, use the suffix automaton.
 * Source: Prefix doubling with radix sort after kactl SuffixArray
 *   (CC0); LCP by Kasai et al. (2001).
 * Status: stress-tested against brute force
 */
#pragma once

struct SuffixArray {
  vi sa, rank, lcp;
  template<class S> SuffixArray(const S& s, int lim = 256) {
    int n = sz(s) + 1; // sorts the cyclic shifts of s + [0]
    vi x(all(s)), y(n), ws(max(n, lim));
    x.push_back(0), sa = y, iota(all(sa), 0);
    for (int j = 0, p = 0; p < n; j = max(1, j * 2), lim = p) {
      p = j, iota(all(y), n - j);
      rep(i, 0, n)
        if (sa[i] >= j) y[p++] = sa[i] - j;
      fill(all(ws), 0);
      rep(i, 0, n) ws[x[i]]++;
      rep(i, 1, lim) ws[i] += ws[i - 1];
      for (int i = n; i--;) sa[--ws[x[y[i]]]] = y[i];
      swap(x, y), p = 1, x[sa[0]] = 0;
      rep(i, 1, n) {
        int a = sa[i - 1], b = sa[i];
        p += y[a] != y[b] || y[a + j] != y[b + j];
        x[b] = p - 1;
      }
    }
    sa.erase(sa.begin()), n--, rank.resize(n);
    rep(i, 0, n) rank[sa[i]] = i;
    lcp.resize(max(n - 1, 0));
    for (int i = 0, k = 0; i < n; i++, k -= k > 0)
      if (rank[i] + 1 < n) {
        int j = sa[rank[i] + 1];
        while (max(i, j) + k < n && s[i + k] == s[j + k]) k++;
        lcp[rank[i]] = k;
      }
  }
};
