#include "../testlib.hpp"
#include "../../content/strings/aho-corasick.hpp"

string rnd_str(int n, char hi) {
  string s(n, 'a');
  for (char& c : s) c = (char)rnd('a', hi);
  return s;
}

bool ends_with(const string& s, const string& t) {
  return sz(t) <= sz(s) && s.compare(sz(s) - sz(t), sz(t), t) == 0;
}

// Checks the automaton against the trie of prefixes: moves, links
// and words, from the definitions.
void check_automaton(const vector<string>& pats, const AhoCorasick& ac) {
  map<string, int> state; // trie node strings, found by walking
  for (auto& p : pats)
    rep(k, 0, sz(p) + 1) {
      int v = 0;
      rep(i, 0, k) v = ac.next[v][p[i] - 'a'];
      string x = p.substr(0, k);
      if (state.count(x)) CHECK_EQ(state[x], v);
      state[x] = v;
    }
  CHECK_EQ(sz(state), sz(ac.next));
  CHECK_EQ(sz(ac.q), sz(ac.next));
  auto longest_suffix = [&](const string& x, int from) {
    rep(k, from, sz(x) + 1) {
      auto it = state.find(x.substr(k));
      if (it != state.end()) return it->second;
    }
    return -1;
  };
  for (auto& [x, v] : state) {
    rep(c, 0, 26) {
      CHECK_EQ(ac.next[v][c], longest_suffix(x + char('a' + c), 0));
    }
    if (v) CHECK_EQ(ac.link[v], longest_suffix(x, 1));
    int w = 0;
    for (auto& p : pats) w += ends_with(x, p);
    CHECK_EQ(ac.words[v], w);
  }
  rep(i, 0, sz(pats)) CHECK_EQ(ac.node[i], state[pats[i]]);
}

int occurrences(const string& s, const string& t) {
  int r = 0;
  for (size_t i = s.find(t); i != string::npos; i = s.find(t, i + 1)) r++;
  return r;
}

int main() {
  rep(it, 0, iters(1500)) {
    char hi = "abcz"[rnd(0, 3)];
    vector<string> pats((int)rnd(1, 8));
    for (auto& p : pats) p = rnd_str((int)rnd(1, 5), hi);
    if (rnd(0, 3) == 0) pats.push_back(pats[0]); // duplicate
    AhoCorasick ac(pats);
    check_automaton(pats, ac);
    string t = rnd_str((int)rnd(0, 40), hi);
    vi occ = ac.count(t);
    ll want_total = 0;
    rep(i, 0, sz(pats)) {
      CHECK_EQ(occ[i], occurrences(t, pats[i]));
      want_total += occ[i];
    }
    // Usage: total matches.
    ll total = 0;
    for (int v = 0; char ch : t) total += ac.words[v = ac.next[v][ch - 'a']];
    CHECK_EQ(total, want_total);
    // Example: list every (pattern, end) match.
    int m = sz(ac.next);
    vi id(m, -1), out(m, -1);
    rep(i, 0, sz(pats)) id[ac.node[i]] = i;
    for (int u : ac.q) out[u] = id[u] >= 0 ? u : out[ac.link[u]];
    set<pii> got, want;
    auto report = [&](int i, int end) {
      CHECK(got.insert({ac.node[i], end}).second);
    };
    int v = 0;
    rep(i, 0, sz(t)) {
      v = ac.next[v][t[i] - 'a'];
      for (int u = out[v]; u >= 0; u = out[ac.link[u]])
        report(id[u], i); // pats[id[u]] ends at t[i]
    }
    rep(i, 0, sz(pats))
      rep(e, 0, sz(t)) {
        if (ends_with(t.substr(0, e + 1), pats[i]))
          want.insert({ac.node[i], e});
      }
    CHECK(got == want);
  }
  // words is ll: totals over a long text exceed int.
  static_assert(is_same_v<decltype(AhoCorasick::words), vector<ll>>);
  // Edge cases: one letter, empty text, every letter.
  CHECK(AhoCorasick({"a"}).count("") == vi{0});
  CHECK(AhoCorasick({"a", "a"}).count("aaa") == vi({3, 3}));
  vector<string> letters;
  rep(c, 0, 26) letters.push_back(string(1, char('a' + c)));
  string abc = "thequickbrownfoxjumpsoverthelazydog";
  vi occ = AhoCorasick(letters).count(abc);
  rep(c, 0, 26) CHECK_EQ(occ[c], (int)count(all(abc), char('a' + c)));
  // Large: random patterns of total length 1e5 (a big trie) on a
  // text of length 1e5, and a, aa, aaa, ... on a^n.
  elapsed();
  vector<string> pats;
  int L = 0;
  while (L < 100'000) {
    pats.push_back(rnd_str((int)rnd(1, 20), 'b'));
    L += sz(pats.back());
  }
  string t = rnd_str(100'000, 'b');
  occ = AhoCorasick(pats).count(t);
  rep(i, 0, 300) CHECK_EQ(occ[i], occurrences(t, pats[i]));
  pats.clear();
  rep(k, 1, 450) pats.push_back(string(k, 'a'));
  occ = AhoCorasick(pats).count(string(100'000, 'a'));
  rep(k, 1, 450) CHECK_EQ(occ[k - 1], 100'000 - k + 1);
  fprintf(stderr, "  large inputs: %.2fs\n", elapsed());
  puts("OK strings/aho-corasick");
}
