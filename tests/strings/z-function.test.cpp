#include "../testlib.hpp"
#include "../../content/strings/z-function.hpp"

string rnd_str(int n, char hi) {
  string s(n, 'a');
  for (char& c : s) c = (char)rnd('a', hi);
  return s;
}

int lcp(const string& a, const string& b) {
  int k = 0;
  while (k < sz(a) && k < sz(b) && a[k] == b[k]) k++;
  return k;
}

vi brute_z(const string& s) {
  vi z(sz(s));
  rep(i, 0, sz(s)) z[i] = lcp(s, s.substr(i));
  return z;
}

vi brute_pi(const string& s) {
  vi pi(sz(s));
  rep(i, 0, sz(s))
    rep(k, 1, i + 1)
      if (s.compare(0, k, s, i + 1 - k, k) == 0) pi[i] = k;
  return pi;
}

int distinct_substrings(const string& s) {
  set<string> subs;
  rep(i, 0, sz(s))
    rep(j, i + 1, sz(s) + 1) subs.insert(s.substr(i, j - i));
  return sz(subs);
}

int main() {
  rep(it, 0, iters(3000)) {
    int n = (int)rnd(0, 24);
    string s = rnd_str(n, "abc"[rnd(0, 2)]);
    vi z = z_function(s);
    CHECK(z == brute_z(s));
    // Usage: matching through pat + '#' + s.
    string pat = rnd_str((int)rnd(1, 4), 'b');
    int m = sz(pat);
    vi zt = z_function(pat + '#' + s);
    rep(i, 0, n) CHECK_EQ(zt[m + 1 + i] == m, s.compare(i, m, pat) == 0);
    rep(i, 0, n) CHECK_EQ(zt[m + 1 + i], lcp(pat, s.substr(i)));
    // Notes and Example: periods (p = n has no z[p]), pi from z.
    rep(p, 1, n + 1) {
      bool period = s.compare(0, n - p, s, p, n - p) == 0;
      CHECK_EQ(p == n || z[p] == n - p, period);
    }
    vi pi(n);
    rep(i, 1, n)
      if (z[i]) pi[i + z[i] - 1] = max(pi[i + z[i] - 1], z[i]);
    for (int i = n - 2; i >= 0; i--) pi[i] = max(pi[i], pi[i + 1] - 1);
    CHECK(pi == brute_pi(s));
    // Notes: distinct substrings, one appended character at a time.
    if (n <= 12) {
      int total = 0;
      rep(len, 1, n + 1) {
        string u = s.substr(0, len);
        reverse(all(u));
        vi zu = z_function(u);
        total += len - (len > 1 ? *max_element(zu.begin() + 1, zu.end()) : 0);
      }
      CHECK_EQ(total, distinct_substrings(s));
    }
  }
  // Any sequence: big values, with -1 as the separator.
  rep(it, 0, iters(500)) {
    vector<ll> a = rnd_vec((int)rnd(0, 20), 0, 1), b = rnd_vec(2, 0, 1);
    string sa, sb;
    for (ll& x : a) sa += char('a' + x), x *= (ll)1e18;
    for (ll& x : b) sb += char('a' + x), x *= (ll)1e18;
    CHECK(z_function(a) == brute_z(sa));
    vector<ll> c = b;
    c.push_back(-1), c.insert(c.end(), all(a));
    vi zc = z_function(c);
    rep(i, 0, sz(a)) CHECK_EQ(zc[3 + i] == 2, sa.compare(i, 2, sb) == 0);
  }
  CHECK(z_function(string()).empty());
  CHECK(z_function(string("a")) == vi{1});
  { // Example: all periods of abaab, including p = n.
    string s = "abaab";
    int n = sz(s);
    vi z = z_function(s), per;
    rep(p, 1, n + 1)
      if (p == n || z[p] == n - p) per.push_back(p);
    CHECK(per == vi({3, 5}));
  }
  // Large: a^n and a random binary string against the prefix check.
  int n = 1'000'000;
  vi z = z_function(string(n, 'a'));
  rep(i, 0, n) CHECK_EQ(z[i], n - i);
  string s = rnd_str(n, 'b');
  z = z_function(s);
  rep(i, 1, n) {
    CHECK(i + z[i] <= n);
    CHECK(s.compare(0, z[i], s, i, z[i]) == 0);
    CHECK(i + z[i] == n || s[z[i]] != s[i + z[i]]);
  }
  puts("OK strings/z-function");
}
