/**
 * Title: Aho-Corasick
 * Description: Automaton of a set of patterns. State `v` is a node
 *   of their trie, that is a prefix of some pattern; after reading
 *   a text, the state is its longest suffix that is such a prefix.
 *   `next[v][c]` is the move on letter `c`, `link[v]` the longest
 *   proper suffix of `v` in the trie, and `words[v]` the number of
 *   patterns that end at the current position when in state `v`.
 *   `count(t)` returns how often each pattern occurs in `t`.
 * Precondition: Non-empty patterns of letters `'a'` to `'z'`;
 *   duplicates are fine.
 * Time: $O(26 L)$ to build, $L$ the total pattern length; $O(|t| +
 *   L)$ for `count(t)`.
 * Usage:
 *   AhoCorasick ac(pats); // vector<string>
 *   vi occ = ac.count(t); // occ[i]: occurrences of pats[i]
 *   ll total = 0; // matches of all patterns in t:
 *   for (int v = 0; char ch : t)
 *     total += ac.words[v = ac.next[v][ch - 'a']];
 * Tweak: Other alphabets: change `26` and the offset `'a'`.
 *   Weighted patterns: add the weight instead of 1 to `words`,
 *   which is `ll` since its sums along links overflow `int`.
 * Example:
 *   // Every match as (pattern, end position). out[v]: nearest
 *   // state on the link path from v, v included, that ends a
 *   // pattern. Duplicate patterns report one index.
 *   int m = sz(ac.next);
 *   vi id(m, -1), out(m, -1);
 *   rep(i, 0, sz(pats)) id[ac.node[i]] = i;
 *   for (int u : ac.q) out[u] = id[u] >= 0 ? u : out[ac.link[u]];
 *   int v = 0;
 *   rep(i, 0, sz(t)) {
 *     v = ac.next[v][t[i] - 'a'];
 *     for (int u = out[v]; u >= 0; u = out[ac.link[u]])
 *       report(id[u], i); // pats[id[u]] ends at t[i]
 *   }
 * Notes: The links form the failure tree: `words[v]` sums the
 *   pattern ends on the path from `v` to the root (BFS order `q`),
 *   and `count` sums the hits over subtrees (reverse BFS order).
 *   Listing every match (Example) can output $Theta(|t| sqrt(L))$
 *   pairs; count when you can. Dynamic programming over the
 *   automaton handles "strings of length `m` that contain no
 *   pattern" (states with `words[v] = 0` only; matrix power for
 *   huge `m`), "maximum total score of matches", and "shortest
 *   string containing all patterns" (BFS over state and a mask of
 *   patterns seen; with `words` replaced by masks OR-ed along
 *   links). Queries "how often does pattern `i` occur in text `j`"
 *   offline: mark the text states and sum over the subtree of the
 *   pattern node in the failure tree (Euler tour and a Fenwick
 *   tree). Cues: many patterns against one text, "forbidden
 *   substrings", counting strings that avoid or contain given
 *   words. One pattern: use the prefix function.
 * Source: Aho, Corasick (1975); goto automaton as in kactl
 *   AhoCorasick (CC0).
 * Status: stress-tested against brute force
 */
#pragma once

struct AhoCorasick {
  vector<array<int, 26>> next; // trie edges, then all moves
  vi link, node, q; // node[i]: state of pats[i]; q: BFS order
  vector<ll> words;
  AhoCorasick(const vector<string>& pats) : next(1) {
    for (auto& p : pats) {
      int v = 0;
      for (char ch : p) {
        int c = ch - 'a';
        if (!next[v][c])
          next[v][c] = sz(next), next.emplace_back();
        v = next[v][c];
      }
      node.push_back(v);
    }
    link.resize(sz(next)), words.resize(sz(next));
    for (int v : node) words[v]++;
    q = {0};
    rep(i, 0, sz(q)) {
      int v = q[i];
      words[v] += v ? words[link[v]] : 0;
      rep(c, 0, 26) {
        int u = next[v][c], w = v ? next[link[v]][c] : 0;
        if (!u) next[v][c] = w;
        else link[u] = w, q.push_back(u);
      }
    }
  }
  vi count(const string& t) {
    vi hits(sz(next)), res;
    int v = 0;
    for (char ch : t) hits[v = next[v][ch - 'a']]++;
    for (int i = sz(q); --i;) hits[link[q[i]]] += hits[q[i]];
    for (int u : node) res.push_back(hits[u]);
    return res;
  }
};
