/**
 * Title: Z-function
 * Description: `z[i]` is the length of the longest common prefix of
 *   `s` and its suffix `s[i, n)`. By convention `z[0] = n`.
 * Time: $O(n)$.
 * Usage:
 *   vi z = z_function(s);
 *   // pat occurs at i in s iff zt[m + 1 + i] == m, m = sz(pat):
 *   vi zt = z_function(pat + '#' + s); // '#' in neither
 * Tweak: Works on any sequence with `==`; for `vector<int>` use a
 *   separator value that occurs nowhere, like -1. Many
 *   implementations use `z[0] = 0`; check before reusing a formula.
 * Example:
 *   // p is a period of s iff p == n or z[p] == n - p:
 *   vi per;
 *   rep(p, 1, n + 1) if (p == n || z[p] == n - p) per.push_back(p);
 *   // LCP of t with every suffix of s: z of t + '#' + s.
 *   vi zt = z_function(t + '#' + s);
 *   int lcp_i = zt[sz(t) + 1 + i]; // lcp(t, s[i, n))
 * Notes: `z[i]` compares one suffix with the whole string; the
 *   prefix function `pi` compares prefixes with their borders. Both
 *   carry the same information: `p < n` is a period iff
 *   `z[p] = n - p`, iff `n - p` is on the border chain of `pi`;
 *   `p = n` is always a period and has no `z[p]`. To get `pi`
 *   from `z`: for each `i >= 1` set
 *   `pi[i + z[i] - 1] = max(.., z[i])`, then from right to left
 *   `pi[i] = max(pi[i], pi[i + 1] - 1)`. Reach for Z when you
 *   need the LCP of `s`, or of a pattern `t`, with every suffix at
 *   once: matching, the longest prefix of `t` starting at each
 *   position, comparing `s[i, n)` with `t` in $O(1)$ (LCP, then
 *   one character), and on reversed strings the longest common
 *   suffix. Distinct substrings in $O(n^2)$: after appending `c`
 *   to `s`, the new substrings are the prefixes of `u`, the
 *   reverse of `s + c`, longer than $max_(i >= 1) z_u [i]$.
 * Source: Gusfield, Algorithms on Strings, Trees and Sequences
 *   (1997); kactl Zfunc (CC0).
 * Status: stress-tested against brute force
 */
#pragma once

template<class S> vi z_function(const S& s) {
  int n = sz(s);
  vi z(n);
  if (n) z[0] = n;
  for (int i = 1, l = 0, r = 0; i < n; i++) {
    if (i < r) z[i] = min(r - i, z[i - l]);
    while (i + z[i] < n && s[z[i]] == s[i + z[i]]) z[i]++;
    if (i + z[i] > r) l = i, r = i + z[i];
  }
  return z;
}
