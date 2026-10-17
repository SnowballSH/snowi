/**
 * Title: Prefix function (KMP)
 * Description: `pi[i]` is the length of the longest proper border
 *   of `s[0, i]`.
 *   `kmp_match(s, pat)` lists the start of every occurrence of `pat`
 *   in `s`, overlapping ones included.
 * Precondition: `pat` is non-empty.
 * Time: $O(n)$; `kmp_match` $O(|s| + |"pat"|)$.
 * Usage:
 *   vi pi = prefix_function(s);
 *   int per = n - pi[n - 1]; // smallest period, n >= 1
 *   vi occ = kmp_match(s, pat); // start positions
 * Example:
 *   // Prefix automaton: nxt[k][c] is the matched length after
 *   // reading c with k characters of s matched (k <= n).
 *   vector<array<int, 26>> nxt(n + 1);
 *   rep(k, 0, n + 1) rep(c, 0, 26) {
 *     if (k < n && s[k] - 'a' == c) nxt[k][c] = k + 1;
 *     else if (k) nxt[k][c] = nxt[pi[k - 1]][c];
 *   }
 *   // cnt[L] = number of occurrences of s[0, L) in s, L >= 1:
 *   vi cnt(n + 1);
 *   rep(i, 0, n) cnt[pi[i]]++;
 *   for (int L = n; L > 0; L--) cnt[pi[L - 1]] += cnt[L];
 *   rep(L, 1, n + 1) cnt[L]++;
 * Notes: The borders of `s[0, L)` form the chain `L -> pi[L - 1] ->
 *   pi[pi[L - 1] - 1] -> ... -> 0`; the step is `k = pi[k - 1]`,
 *   not `pi[k]`. They are the ancestors of `L` in the failure tree
 *   whose nodes are `0..n` and where the parent of `L` is
 *   `pi[L - 1]`; parents are smaller, so tree DP runs as a plain
 *   loop over `L` (as in the counting example) and no recursion is
 *   needed. Occurrences of every prefix of `s` in another text `t`:
 *   run `kmp_match`'s loop, add one to `cnt[k]` after each character
 *   (before falling back when `k = n`), then push counts to parents
 *   as above, without the final `+1`. The prefix automaton turns
 *   "the text must (not) contain `s`" into a DP over the matched
 *   length: count strings of length `m` avoiding `s` in
 *   $O(m n sigma)$, or with a matrix power for huge `m`. Some
 *   statements count the whole string as a border and some do not;
 *   check the samples. Cues: "prefix equals suffix", "shortest
 *   string whose repetition gives `s`", "the same piece at both
 *   ends".
 * Source: Knuth, Morris, Pratt (1977); kactl KMP (CC0).
 * Status: stress-tested against brute force
 */
#pragma once

template<class S> vi prefix_function(const S& s) {
  vi pi(sz(s));
  rep(i, 1, sz(s)) {
    int k = pi[i - 1];
    while (k && s[i] != s[k]) k = pi[k - 1];
    pi[i] = k + (s[i] == s[k]);
  }
  return pi;
}
template<class S> vi kmp_match(const S& s, const S& pat) {
  vi pi = prefix_function(pat), res;
  for (int i = 0, k = 0; i < sz(s); i++) {
    while (k && s[i] != pat[k]) k = pi[k - 1];
    if (s[i] == pat[k] && ++k == sz(pat))
      res.push_back(i - k + 1), k = pi[k - 1];
  }
  return res;
}
