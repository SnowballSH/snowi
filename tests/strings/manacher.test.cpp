#include "../testlib.hpp"
#include "../../content/strings/manacher.hpp"

string rnd_str(int n, char hi) {
  string s(n, 'a');
  for (char& c : s) c = (char)rnd('a', hi);
  return s;
}

bool is_pal(const string& s, int l, int r) {
  rep(i, l, r)
    if (s[i] != s[l + r - 1 - i]) return false;
  return true;
}

int main() {
  rep(it, 0, iters(3000)) {
    int n = (int)rnd(0, 24);
    string s = rnd_str(n, "abc"[rnd(0, 2)]);
    vi plen = manacher(s);
    CHECK_EQ(sz(plen), n ? 2 * n - 1 : 0);
    // Brute force: longest palindrome per center, every substring.
    rep(c, 0, 2 * n - 1) {
      int best = 0;
      rep(l, 0, n) {
        int r = c + 1 - l;
        if (l < r && r <= n && is_pal(s, l, r)) best = max(best, r - l);
      }
      CHECK_EQ(plen[c], best);
    }
    rep(l, 0, n)
      rep(r, l + 1, n + 1) CHECK_EQ(plen[l + r - 1] >= r - l, is_pal(s, l, r));
    // Usage: longest palindromic substring.
    if (n) {
      int c = max_element(all(plen)) - begin(plen);
      string best = s.substr((c + 1 - plen[c]) / 2, plen[c]);
      CHECK(is_pal(best, 0, sz(best)));
      CHECK(s.find(best) != string::npos);
      rep(l, 0, n)
        rep(r, l + 1, n + 1) CHECK(!is_pal(s, l, r) || r - l <= sz(best));
    }
    // Example and Notes: count, start and end ranges.
    ll total = 0;
    for (int x : plen) total += (x + 1) / 2;
    vi starts(n + 1), ends(n + 1), bs(n), be(n);
    rep(c, 0, sz(plen)) {
      int p = plen[c];
      starts[(c + 1 - p) / 2]++, starts[c / 2 + 1]--;
      ends[(c + 1) / 2]++, ends[(c + p - 1) / 2 + 1]--;
    }
    ll want = 0;
    rep(l, 0, n)
      rep(r, l + 1, n + 1)
        if (is_pal(s, l, r)) { want++, bs[l]++, be[r - 1]++; }
    CHECK_EQ(total, want);
    rep(i, 0, n) {
      starts[i + 1] += starts[i], ends[i + 1] += ends[i];
      CHECK_EQ(starts[i], bs[i]);
      CHECK_EQ(ends[i], be[i]);
    }
    int pre = 0;
    rep(r, 1, n + 1)
      if (plen[r - 1] >= r) pre = r;
    rep(r, 1, n + 1)
      if (is_pal(s, 0, r)) CHECK(r <= pre);
    CHECK(n == 0 || is_pal(s, 0, pre));
  }
  // Sequences of big ints.
  rep(it, 0, iters(500)) {
    vector<ll> a = rnd_vec((int)rnd(0, 20), 0, 1);
    string sa;
    for (ll& x : a) sa += char('a' + x), x = x * (ll)4e18 - (ll)2e18;
    CHECK(manacher(a) == manacher(sa));
  }
  CHECK(manacher(string()).empty());
  CHECK(manacher(string("a")) == vi{1});
  // Large: a^n, and a random text checked against hashing-free
  // properties (a maximal palindrome cannot grow).
  int n = 1'000'000;
  vi plen = manacher(string(n, 'a'));
  rep(c, 0, 2 * n - 1) CHECK_EQ(plen[c], c < n ? c + 1 : 2 * n - 1 - c);
  string s = rnd_str(n, 'b');
  plen = manacher(s);
  rep(c, 0, 2 * n - 1) {
    int l = (c + 1 - plen[c]) / 2, r = l + plen[c];
    CHECK(l + r - 1 == c);
    CHECK(l == 0 || r == n || s[l - 1] != s[r]);
    if (plen[c] <= 64) CHECK(is_pal(s, l, r));
  }
  puts("OK strings/manacher");
}
