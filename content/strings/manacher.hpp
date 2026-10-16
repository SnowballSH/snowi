/**
 * Title: Manacher
 * Description: All maximal palindromes. A palindrome `s[l, r)` has
 *   center `c = l + r - 1` in $[0, 2n - 1)$: even `c` is the
 *   character `c/2`, odd `c` the gap after character `(c-1)/2`.
 *   `plen[c]` is the length of the longest palindrome with center
 *   `c`.
 * Time: $O(n)$.
 * Usage:
 *   vi plen = manacher(s); // size 2n - 1
 *   bool pal = plen[l + r - 1] >= r - l; // s[l, r), l < r
 *   int c = max_element(all(plen)) - begin(plen); // n >= 1
 *   string best = s.substr((c + 1 - plen[c]) / 2, plen[c]);
 * Example:
 *   // Number of palindromic substrings (with repeats):
 *   ll total = 0;
 *   for (int x : plen) total += (x + 1) / 2;
 *   // Longest palindromic prefix of s: largest r with
 *   // plen[r - 1] >= r.
 * Notes: The palindromes with center `c` have lengths `p, p - 2,
 *   ...` down to 1 or 2, where `p = plen[c]`. They start at every
 *   position of $[(c + 1 - p)/2, floor(c/2)]$ and end at every
 *   position of $[ceil(c/2), (c + p - 1)/2]$, so the number of
 *   palindromes starting or ending at each position comes from a
 *   difference array (pairs of disjoint palindromes, cuts into
 *   palindromes). A string has at most `n` distinct palindromic
 *   substrings: each position ends at most one new one, the longest
 *   palindromic suffix so far; to count or list them use the
 *   eertree. Cues: many "is `s[l, r)` a palindrome" queries, "split
 *   into two palindromes" (a prefix and a suffix test), "longest
 *   palindromic substring".
 * Source: Manacher (1975); one-array form after the classic
 *   `#`-interleaved version.
 * Status: stress-tested against brute force
 */
#pragma once

template<class S> vi manacher(const S& s) {
  int n = sz(s), lo = 0, hi = -1; // rightmost found: s[lo..hi]
  vi plen(max(2 * n - 1, 0));
  rep(c, 0, sz(plen)) {
    int k = 1 - c % 2, mirror = 2 * (lo + hi) - c;
    if (c <= 2 * hi) k = min(plen[mirror], 2 * hi - c + 1);
    int i = (c - k + 1) / 2, j = (c + k - 1) / 2;
    while (i > 0 && j + 1 < n && s[i - 1] == s[j + 1]) i--, j++;
    plen[c] = j - i + 1;
    if (j > hi) lo = i, hi = j;
  }
  return plen;
}
