#include "../testlib.hpp"
#include "../../content/strings/suffix-automaton.hpp"
#include "../../content/strings/suffix-array.hpp"

string rnd_str(int n, char hi) {
  string s(n, 'a');
  for (char& c : s) c = (char)rnd('a', hi);
  return s;
}

int occurrences(const string& s, const string& t) {
  int r = 0;
  for (size_t i = s.find(t); i != string::npos; i = s.find(t, i + 1)) r++;
  return r;
}

// Every substring maps to the state of its end-position class, the
// classes are exactly the states (minimality), and counts match.
template<class SAM> void check_brute(const string& s, SAM& sam) {
  int n = sz(s);
  map<string, vi> endpos;
  rep(l, 0, n)
    rep(r, l + 1, n + 1) endpos[s.substr(l, r - l)].push_back(r);
  CHECK_EQ(sam.distinct(), (ll)sz(endpos));
  map<vi, int> state_of;
  map<int, vi> class_of;
  for (auto& [t, pos] : endpos) {
    int v = sam.find(t);
    CHECK(v > 0);
    CHECK(sam.len[sam.link[v]] < sz(t) && sz(t) <= sam.len[v]);
    if (state_of.count(pos)) CHECK_EQ(state_of[pos], v);
    if (class_of.count(v)) CHECK(class_of[v] == pos);
    state_of[pos] = v, class_of[v] = pos;
  }
  CHECK_EQ(sz(sam.len), sz(state_of) + 1);
  CHECK(sz(sam.len) <= max(n + 1, 2 * n - 1));
  CHECK_EQ(sam.find(""), 0);
  sam.count();
  for (auto& [t, pos] : endpos) CHECK_EQ(sam.cnt[sam.find(t)], sz(pos));
  rep(it, 0, 10) {
    string t = rnd_str((int)rnd(1, 4), 'c');
    int v = sam.find(t), occ = v < 0 ? 0 : sam.cnt[v];
    CHECK_EQ(v >= 0, endpos.count(t) > 0);
    CHECK_EQ(occ, occurrences(s, t));
  }
}

int longest_common(const string& s, const string& t) {
  int best = 0;
  rep(i, 0, sz(t))
    rep(j, i + 1, sz(t) + 1) {
      if (s.find(t.substr(i, j - i)) == string::npos) break;
      best = max(best, j - i);
    }
  return best;
}

// Online: distinct substrings after each letter, by distinct() and
// by the increments of extend, with finds of letters not seen yet
// in between (in the map variant they insert 0 edges).
template<class SAM> void check_online(const string& s) {
  SAM on;
  set<string> subs;
  ll d = 0;
  rep(r, 0, sz(s)) {
    rep(it, 0, 3) {
      string t = rnd_str((int)rnd(1, 3), 'd');
      CHECK_EQ(on.find(t) >= 0, subs.count(t) > 0);
    }
    on.extend(s[r] - 'a');
    d += on.len[on.last] - on.len[on.link[on.last]];
    rep(l, 0, r + 1) subs.insert(s.substr(l, r + 1 - l));
    CHECK_EQ(on.distinct(), (ll)sz(subs));
    CHECK_EQ(d, (ll)sz(subs));
  }
  check_brute(s, on);
}

template<class SAM> void check_all(const string& s, char hi) {
  SAM sam(s);
  check_brute(s, sam);
  check_online<SAM>(s);
  // Example: longest common substring.
  string t = rnd_str((int)rnd(0, 14), hi);
  int v = 0, l = 0, best = 0;
  for (char ch : t) {
    int c = ch - 'a';
    while (v && !sam.next[v][c]) v = sam.link[v], l = sam.len[v];
    if (sam.next[v][c]) v = sam.next[v][c], l++;
    best = max(best, l);
  }
  CHECK_EQ(best, longest_common(s, t));
}

using MapSAM = SuffixAutomaton<map<int, int>>; // the Tweak

int main() {
  rep(it, 0, iters(1500)) {
    int n = (int)rnd(0, 14);
    char hi = "abcz"[rnd(0, 3)];
    string s = rnd_str(n, hi);
    check_all<SuffixAutomaton<>>(s, hi);
    if (it % 3 == 0) check_all<MapSAM>(s, hi);
  }
  // Map variant: find() of a missing letter, then extend() by it.
  MapSAM M("ab");
  CHECK_EQ(M.find("bc"), -1);
  M.extend('c' - 'a');
  CHECK_EQ(M.distinct(), 6LL);
  CHECK(M.find("abc") > 0 && M.find("bc") > 0);
  // Map variant with arbitrary int letters.
  rep(it, 0, iters(300)) {
    vi a((int)rnd(0, 12));
    for (int& x : a) x = (int)rnd(0, 2) * 1'000'000'000 - 1;
    MapSAM sam;
    set<vi> subs;
    rep(r, 0, sz(a)) {
      sam.extend(a[r]);
      rep(l, 0, r + 1) subs.insert(vi(a.begin() + l, a.begin() + r + 1));
      CHECK_EQ(sam.distinct(), (ll)sz(subs));
    }
  }
  // Large: a^n, a Fibonacci word and random texts; distinct counts
  // against the suffix array, occurrence counts against find.
  elapsed();
  int n = 100'000;
  SuffixAutomaton A(string(n, 'a'));
  CHECK_EQ(A.distinct(), (ll)n);
  A.count();
  CHECK_EQ(A.cnt[A.find(string(7, 'a'))], n - 6);
  string f = "a", g = "ab";
  while (sz(g) < n) tie(f, g) = pair(g, g + f);
  for (const string& s : {g, rnd_str(n, 'b'), rnd_str(n, 'z')}) {
    SuffixAutomaton sam;
    SuffixArray S(s);
    ll d = (ll)sz(s) * (sz(s) + 1) / 2, online = 0;
    for (int h : S.lcp) d -= h;
    for (char ch : s) {
      sam.extend(ch - 'a');
      online += sam.len[sam.last] - sam.len[sam.link[sam.last]];
    }
    CHECK_EQ(sam.distinct(), d);
    CHECK_EQ(online, d);
    CHECK_EQ(MapSAM(s).distinct(), d);
    CHECK(sz(sam.len) <= 2 * sz(s) - 1);
    sam.count();
    rep(it, 0, 30) {
      int len = (int)rnd(1, 12), at = (int)rnd(0, sz(s) - len);
      string t = it % 2 ? s.substr(at, len) : rnd_str(len, 'c');
      int v = sam.find(t);
      CHECK_EQ(v < 0 ? 0 : sam.cnt[v], occurrences(s, t));
    }
  }
  fprintf(stderr, "  large inputs: %.2fs\n", elapsed());
  puts("OK strings/suffix-automaton");
}
