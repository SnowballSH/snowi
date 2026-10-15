#include "../testlib.hpp"
#include "../../content/strings/eertree.hpp"

string rnd_str(int n, char hi) {
  string s(n, 'a');
  for (char& c : s) c = (char)rnd('a', hi);
  return s;
}

bool is_pal(const string& s) { return equal(all(s), s.rbegin()); }

int occurrences(const string& s, const string& t) {
  int r = 0;
  for (size_t i = s.find(t); i != string::npos; i = s.find(t, i + 1)) r++;
  return r;
}

// Node strings from the edges (c + v + c), checked against every
// palindromic substring, the links, suf and the counts.
void check_brute(const string& s, Eertree& et) {
  int n = sz(s), m = sz(et.len);
  vector<string> str(m);
  set<string> pals, seen;
  rep(l, 0, n)
    rep(r, l + 1, n + 1) {
      string t = s.substr(l, r - l);
      if (is_pal(t)) pals.insert(t);
    }
  CHECK_EQ(m - 2, sz(pals));
  rep(v, 0, m)
    rep(c, 0, 26) {
      int w = et.next[v][c];
      if (!w) continue;
      CHECK(w >= 2 && w > v);
      string mid = v ? str[v] : "";
      str[w] = char('a' + c) + (et.len[v] < 0 ? "" : mid + char('a' + c));
      CHECK(seen.insert(str[w]).second);
    }
  CHECK(seen == pals);
  CHECK_EQ(et.len[0], -1);
  CHECK_EQ(et.len[1], 0);
  rep(v, 2, m) {
    CHECK_EQ(et.len[v], sz(str[v]));
    int lk = et.link[v], want = 1;
    rep(k, 1, sz(str[v]))
      if (is_pal(str[v].substr(k))) {
        want = (int)(find(all(str), str[v].substr(k)) - str.begin());
        break;
      }
    CHECK_EQ(lk, want);
    CHECK(lk < v);
  }
  CHECK_EQ(sz(et.suf), n);
  rep(i, 0, n) {
    int best = 0;
    rep(l, 0, i + 1)
      if (is_pal(s.substr(l, i + 1 - l))) {
        best = i + 1 - l;
        break;
      }
    CHECK_EQ(et.len[et.suf[i]], best);
    CHECK(str[et.suf[i]] == s.substr(i + 1 - best, best));
  }
  et.count();
  rep(v, 2, m) CHECK_EQ(et.cnt[v], occurrences(s, str[v]));
  // Example: max of occurrences times length.
  ll best = 0, want_best = 0;
  rep(v, 2, sz(et.len)) best = max(best, (ll)et.cnt[v] * et.len[v]);
  for (auto& t : pals)
    want_best = max(want_best, (ll)occurrences(s, t) * sz(t));
  CHECK_EQ(best, want_best);
  // Example: palindromic substrings with repeats.
  vi dep(sz(et.len));
  rep(v, 2, sz(et.len)) dep[v] = dep[et.link[v]] + 1;
  ll total = 0, want = 0;
  for (int v : et.suf) total += dep[v];
  rep(l, 0, n)
    rep(r, l + 1, n + 1) want += is_pal(s.substr(l, r - l));
  CHECK_EQ(total, want);
}

int main() {
  rep(it, 0, iters(2000)) {
    int n = (int)rnd(0, 20);
    string s = rnd_str(n, "abcz"[rnd(0, 3)]);
    Eertree et(s);
    check_brute(s, et);
    Eertree on; // online: distinct palindromes of each prefix
    rep(i, 0, n) {
      on.add(s[i] - 'a');
      set<string> pals;
      rep(l, 0, i + 1)
        rep(r, l + 1, i + 2) {
          string t = s.substr(l, r - l);
          if (is_pal(t)) pals.insert(t);
        }
      CHECK_EQ(sz(on.len) - 2, sz(pals));
    }
  }
  // Large: a^n, a Fibonacci word (many palindromes) and random
  // texts; totals and longest suffixes by center expansion.
  elapsed();
  int n = 100'000;
  Eertree A(string(n, 'a'));
  CHECK_EQ(sz(A.len) - 2, n);
  A.count();
  rep(v, 2, n + 2) CHECK_EQ(A.cnt[v], n - A.len[v] + 1);
  ll best = 0; // Example: L * (n - L + 1) beyond int
  rep(v, 2, sz(A.len)) best = max(best, (ll)A.cnt[v] * A.len[v]);
  CHECK_EQ(best, 2'500'050'000LL);
  string f = "a", g = "ab";
  while (sz(g) < n) tie(f, g) = pair(g, g + f);
  for (const string& s : {g, rnd_str(n, 'b'), rnd_str(n, 'z')}) {
    int m = sz(s);
    Eertree et(s);
    CHECK(sz(et.len) - 2 <= m);
    et.count();
    ll total = 0, want = 0;
    rep(v, 2, sz(et.len)) total += et.cnt[v];
    vi longest(m);
    rep(c, 0, 2 * m - 1) { // all palindromes by center expansion
      int l = c / 2, r = (c + 1) / 2;
      for (; l >= 0 && r < m && s[l] == s[r]; l--, r++)
        want++, longest[r] = max(longest[r], r - l + 1);
    }
    CHECK_EQ(total, want);
    rep(i, 0, m) CHECK_EQ(et.len[et.suf[i]], longest[i]);
  }
  fprintf(stderr, "  large inputs: %.2fs\n", elapsed());
  puts("OK strings/eertree");
}
