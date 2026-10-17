/**
 * Title: Suffix automaton
 * Description: Online automaton of the substrings of `s`. Each
 *   state `v` is a class of substrings with the same end positions:
 *   the suffixes of its longest member with lengths
 *   `len[link[v]] + 1` to `len[v]`. After `count()`, `cnt[v]` is how
 *   often each string of `v` occurs in `s`. Each `extend` adds
 *   `len[last] - len[link[last]]` new distinct substrings.
 * Precondition: Letters `'a'` to `'z'`.
 * Time: $O(n)$ amortized, at most $2n$ states and $3n$ edges;
 *   memory $2n dot 26$ ints.
 * Usage:
 *   SuffixAutomaton sam(s); // or sam.extend(c), c in [0, 26)
 *   int v = sam.find(t); // -1: t is not a substring of s
 *   sam.count(); // once, after the last extend
 *   int occ = v < 0 ? 0 : sam.cnt[v]; // t non-empty
 *   ll d = sam.distinct(); // distinct non-empty, O(n)
 * Tweak: Large alphabet or tight memory:
 *   `SuffixAutomaton<map<int, int>>`, passing any `int` to
 *   `extend`. Missing edges then read as 0 (`next[v][c]` inserts them): test edges by value, never
 *   with `count(c)`.
 * Example:
 *   // Longest common substring of s and t: walk t, keeping the
 *   // longest suffix of the prefix read so far that occurs in s.
 *   SuffixAutomaton sam(s);
 *   int v = 0, l = 0, best = 0;
 *   for (char ch : t) {
 *     int c = ch - 'a';
 *     while (v && !sam.next[v][c]) v = sam.link[v], l = sam.len[v];
 *     if (sam.next[v][c]) v = sam.next[v][c], l++;
 *     best = max(best, l);
 *   }
 * Notes: The links form a tree, the suffix tree of the reversed
 *   string: moving to `link[v]` drops letters from the front.
 *   Sorting the states by `len` orders both the links and the edges
 *   topologically, as in `count()`. Distinct substrings starting
 *   from each state: `ways[v] = 1 + sum ways[next[v][c]]` in
 *   decreasing `len`; walk down with it for the `k`-th smallest
 *   substring. First occurrence: record `len[cur] - 1` when a state
 *   is created and copy it into clones; `t` then starts at
 *   `first[v] - |t| + 1`. The end positions of `v` are those of the
 *   non-clone states in its link subtree (small-to-large merging
 *   for position queries). Matching statistics (Example) take
 *   $O(|t|)$ for any text `t`. Several strings: join them with
 *   distinct separators (a larger alphabet). Cues: distinct
 *   substrings of a growing string (sum the increments of `extend`;
 *   `distinct()` after every letter is $O(n^2)$), occurrences of
 *   many query strings, longest common substring, `k`-th distinct
 *   substring. Prefer the suffix array when memory is tight or
 *   when LCPs of arbitrary suffixes are needed.
 * Source: Blumer et al. (1985); construction as in cp-algorithms.
 * Status: stress-tested against brute force
 */
#pragma once

template<class E = array<int, 26>> struct SuffixAutomaton {
  vi len, link, cnt; // cnt: see count()
  vector<E> next;    // 0 means no edge
  int last = 0;
  SuffixAutomaton(const string& s = "") {
    add(0, -1, 0, {});
    for (char ch : s) extend(ch - 'a');
  }
  int add(int l, int lk, int c, E nx) {
    len.push_back(l), link.push_back(lk), cnt.push_back(c);
    next.push_back(nx);
    return sz(len) - 1;
  }
  void extend(int c) {
    int p = last, cur = last = add(len[p] + 1, 0, 1, {});
    for (; p != -1 && !next[p][c]; p = link[p])
      next[p][c] = cur;
    if (p == -1) return;
    int q = next[p][c];
    if (len[p] + 1 == len[q]) {
      link[cur] = q;
      return;
    }
    int cl = add(len[p] + 1, link[q], 0, next[q]);
    for (; p != -1 && next[p][c] == q; p = link[p])
      next[p][c] = cl;
    link[q] = link[cur] = cl;
  }
  void count() { // cnt[v] = number of end positions of v
    vi b(sz(len)), ord(sz(len)); // ord: states by len
    for (int l : len) b[l]++;
    rep(i, 1, sz(b)) b[i] += b[i - 1];
    for (int v = sz(len); v--;) ord[--b[len[v]]] = v;
    for (int i = sz(ord); --i;)
      cnt[link[ord[i]]] += cnt[ord[i]];
  }
  int find(const string& t) {
    int v = 0;
    for (char ch : t)
      if (!(v = next[v][ch - 'a'])) return -1;
    return v;
  }
  ll distinct() {
    ll r = 0;
    rep(v, 1, sz(len)) r += len[v] - len[link[v]];
    return r;
  }
};
