/**
 * Title: Lyndon factorization
 * Description: A Lyndon word is strictly smaller than all its
 *   proper suffixes. `duval(s)` splits `s` into Lyndon words
 *   $w_1 >= w_2 >= dots >= w_k$ (the split is unique) and returns
 *   their start indices. `min_rotation(s)` returns the smallest
 *   `i` for which the rotation `s[i, n) + s[0, i)` is
 *   lexicographically smallest.
 * Time: $O(n)$.
 * Usage:
 *   vi st = duval(s); // start of each factor
 *   int suf = st.back(); // smallest suffix: s[suf, n), n >= 1
 *   st.push_back(n); // now factor j is s[st[j], st[j + 1])
 *   int r = min_rotation(s);
 *   rotate(s.begin(), s.begin() + r, s.end());
 * Tweak: Largest rotation: the smallest one with the order flipped
 *   (`-x` for numbers, `'z' + 'a' - c` for letters).
 * Example:
 *   // duval("abbabbaba") = {0, 3, 6, 8}: abb, abb, ab, a.
 *   // min_rotation("baba") = 1 (abab, the first of 1 and 3).
 *   // Is s a Lyndon word? sz(duval(s)) == 1.
 * Notes: Duval keeps `s[i, j)` as $w^m w'$ with $w$ Lyndon of
 *   length `j - k` and $w'$ a proper prefix of $w$: `s[j] > s[k]`
 *   merges everything into one Lyndon word, `s[j] == s[k]` extends
 *   the period, and `s[j] < s[k]` or the end emits the $m$ copies
 *   of $w$. Facts: a Lyndon word has no non-empty border, so it
 *   is primitive; if $u < v$ are Lyndon, $u v$ is Lyndon; `s`
 *   is Lyndon iff it is strictly smaller than all its other
 *   rotations. The minimal rotation is the canonical form of a
 *   cyclic string: two are equal up to rotation iff their minimal
 *   rotations are, and its start is unique iff `s` is primitive.
 *   Pitfall: the smallest suffix is not always the best start when
 *   more text follows: `a` is the smallest suffix of `aba`, yet
 *   `abac < ac`. Cues: "smallest rotation", "canonical form of a
 *   necklace", "lexicographically smallest result of a cut",
 *   "count substrings smaller than all their suffixes".
 * Source: Duval (1983); UCup Packet 19.
 * Status: stress-tested against brute force
 */
#pragma once

template<class S> vi duval(const S& s) {
  int n = sz(s), i = 0;
  vi res;
  while (i < n) {
    int j = i + 1, k = i;
    for (; j < n && s[k] <= s[j]; j++)
      k = s[k] < s[j] ? i : k + 1;
    for (; i <= k; i += j - k) res.push_back(i);
  }
  return res;
}
template<class S> int min_rotation(const S& s) {
  int n = sz(s), i = 0, ans = 0;
  while (i < n) {
    ans = i;
    int j = i + 1, k = i;
    for (; j < 2 * n && s[k % n] <= s[j % n]; j++)
      k = s[k % n] < s[j % n] ? i : k + 1;
    while (i <= k) i += j - k;
  }
  return ans;
}
