/**
 * Title: Eertree
 * Description: Palindromic tree, built letter by letter. Its nodes
 *   are the distinct non-empty palindromic substrings of `s` plus
 *   two roots: node 0 of length $-1$ and node 1, the empty string.
 *   `next[v][c]` is the palindrome `c + v + c` (from node 0: just
 *   `c`) and `link[v]` the longest proper palindromic suffix of
 *   `v`. `suf[i]` is the node of the longest palindromic suffix of
 *   the prefix ending at `s[i]`.
 * Precondition: Letters `'a'` to `'z'`.
 * Time: $O(n)$ amortized; at most $n + 2$ nodes.
 * Usage:
 *   Eertree et(s); // or et.add(c), c in [0, 26)
 *   int distinct = sz(et.len) - 2; // palindromic substrings
 *   int longest = et.len[et.suf[i]]; // ending at s[i]
 *   et.count(); // then et.cnt[v]: occurrences of node v >= 2
 * Example:
 *   // Palindromic substrings with repeats: dep[v] is the number
 *   // of palindromic suffixes of v, the nodes >= 2 on its link
 *   // path (links point to smaller nodes).
 *   vi dep(sz(et.len));
 *   rep(v, 2, sz(et.len)) dep[v] = dep[et.link[v]] + 1;
 *   ll total = 0;
 *   for (int v : et.suf) total += dep[v]; // ending at each i
 *   // Max of occurrences times length, in ll: it exceeds int
 *   // already for n = 1e5.
 *   et.count();
 *   ll best = 0;
 *   rep(v, 2, sz(et.len))
 *     best = max(best, (ll)et.cnt[v] * et.len[v]);
 * Notes: A string of length `n` has at most `n` distinct non-empty
 *   palindromic substrings, since each letter adds at most one new
 *   one: its longest palindromic suffix. Nodes are created in order
 *   of first occurrence, and `link[v] < v`, so loops over `v` in
 *   increasing order are top-down along links and decreasing order
 *   is bottom-up (as in `count()`). The lengths of the palindromic
 *   suffixes of any string form $O(log n)$ arithmetic progressions;
 *   with `diff[v] = len[v] - len[link[v]]` and the series link (the
 *   first node up the link path with a different `diff`) this gives
 *   $O(n log n)$ dynamic programming over partitions into
 *   palindromes (minimum number of pieces, number of ways). For "is
 *   `s[l, r)` a palindrome" queries use Manacher instead. Cues:
 *   "distinct palindromic substrings", "occurrences of each
 *   palindrome" (the number of occurrences times the length,
 *   maximized: Example, in `ll`), palindromes ending at each
 *   position, palindromic factorization.
 * Source: Rubinchik, Shur, "EERTREE" (2015).
 * Status: stress-tested against brute force
 */
#pragma once

struct Eertree {
  vi len{-1, 0}, link{0, 0}, suf, cnt;
  vi s{-1};                    // -1, then the letters
  vector<array<int, 26>> next; // 0 means no edge
  Eertree(const string& t = "") : next(2) {
    for (char ch : t) add(ch - 'a');
  }
  int up(int v) { // c+up(v)+c is a palindrome, c = s.back()
    int i = sz(s) - 1;
    while (s[i - len[v] - 1] != s[i]) v = link[v];
    return v;
  }
  void add(int c) {
    s.push_back(c);
    int v = up(suf.empty() ? 1 : suf.back());
    if (!next[v][c]) {
      link.push_back(len[v] < 0 ? 1 : next[up(link[v])][c]);
      len.push_back(len[v] + 2), next.emplace_back();
      next[v][c] = sz(len) - 1;
    }
    suf.push_back(next[v][c]);
  }
  void count() { // cnt[v]: occurrences of palindrome v
    cnt.assign(sz(len), 0);
    for (int v : suf) cnt[v]++;
    for (int v = sz(len); --v > 1;) cnt[link[v]] += cnt[v];
  }
};
