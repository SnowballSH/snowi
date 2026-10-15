#include "../testlib.hpp"
#include "../../content/strings/prefix-function.hpp"

string rnd_str(int n, char hi) {
  string s(n, 'a');
  for (char& c : s) c = (char)rnd('a', hi);
  return s;
}

vi brute_pi(const string& s) {
  vi pi(sz(s));
  rep(i, 0, sz(s))
    rep(k, 1, i + 1)
      if (s.compare(0, k, s, i + 1 - k, k) == 0) pi[i] = k;
  return pi;
}

vi brute_match(const string& s, const string& pat) {
  vi res;
  rep(i, 0, sz(s) - sz(pat) + 1)
    if (s.compare(i, sz(pat), pat) == 0) res.push_back(i);
  return res;
}

int brute_period(const string& s) {
  rep(p, 1, sz(s) + 1) {
    bool ok = true;
    rep(i, 0, sz(s) - p) ok &= s[i] == s[i + p];
    if (ok) return p;
  }
  return -1;
}

// Notes: occurrences of every prefix of s in another text t.
void check_text_counts(const string& s, const vi& pi, const string& t) {
  int n = sz(s), k = 0;
  vi cnt(n + 1);
  for (char c : t) {
    while (k && c != s[k]) k = pi[k - 1];
    if (c == s[k]) k++;
    cnt[k]++;
    if (k == n) k = pi[k - 1];
  }
  for (int L = n; L > 0; L--) cnt[pi[L - 1]] += cnt[L];
  rep(L, 1, n + 1) CHECK_EQ(cnt[L], sz(brute_match(t, s.substr(0, L))));
}

// The code of the header's Example.
void check_example(const string& s, const vi& pi) {
  int n = sz(s);
  vector<array<int, 26>> nxt(n + 1);
  rep(k, 0, n + 1)
    rep(c, 0, 26) {
      if (k < n && s[k] - 'a' == c) nxt[k][c] = k + 1;
      else if (k) nxt[k][c] = nxt[pi[k - 1]][c];
    }
  vi cnt(n + 1);
  rep(i, 0, n) cnt[pi[i]]++;
  for (int L = n; L > 0; L--) cnt[pi[L - 1]] += cnt[L];
  rep(L, 1, n + 1) cnt[L]++;
  // nxt[k][c]: longest prefix of s that is a suffix of s[0, k) + c
  rep(k, 0, n + 1)
    rep(c, 0, 3) {
      string t = s.substr(0, k) + char('a' + c);
      int want = 0;
      rep(L, 1, min(n, sz(t)) + 1)
        if (t.compare(sz(t) - L, L, s, 0, L) == 0) want = L;
      CHECK_EQ(nxt[k][c], want);
    }
  rep(L, 1, n + 1) CHECK_EQ(cnt[L], sz(brute_match(s, s.substr(0, L))));
}

int main() {
  rep(it, 0, iters(3000)) {
    int n = (int)rnd(0, 24);
    char hi = "abc"[rnd(0, 2)];
    string s = rnd_str(n, hi);
    vi pi = prefix_function(s);
    CHECK(pi == brute_pi(s));
    if (n) CHECK_EQ(n - pi[n - 1], brute_period(s));
    check_example(s, pi);
    if (n) check_text_counts(s, pi, rnd_str((int)rnd(0, 30), hi));
    string pat = rnd_str((int)rnd(1, 5), hi);
    if (rnd(0, 3) == 0 && n) {
      int l = (int)rnd(0, n - 1);
      pat = s.substr(l, rnd(1, n - l));
    }
    CHECK(kmp_match(s, pat) == brute_match(s, pat));
  }
  // Sequences of ints, values far outside the char range.
  rep(it, 0, iters(500)) {
    vector<ll> a = rnd_vec((int)rnd(0, 20), -1, 1);
    vector<ll> b = rnd_vec((int)rnd(1, 3), -1, 1);
    string sa, sb;
    for (ll x : a) sa += char('b' + x);
    for (ll x : b) sb += char('b' + x);
    for (ll& x : a) x *= (ll)1e18;
    for (ll& x : b) x *= (ll)1e18;
    CHECK(prefix_function(a) == brute_pi(sa));
    CHECK(kmp_match(a, b) == brute_match(sa, sb));
  }
  CHECK(kmp_match(string("ab"), string("abc")).empty());
  CHECK(kmp_match(string(""), string("a")).empty());
  CHECK(kmp_match(string("a"), string("a")) == vi{0});
  // Large: a^n (longest chains) and a random text.
  int n = 1'000'000;
  string s(n, 'a');
  vi pi = prefix_function(s);
  rep(i, 0, n) CHECK_EQ(pi[i], i);
  CHECK_EQ(sz(kmp_match(s, string(1000, 'a'))), n - 999);
  s = rnd_str(n, 'b');
  string pat = "abbabaab";
  vi occ = kmp_match(s, pat), want;
  for (size_t p = s.find(pat); p != string::npos; p = s.find(pat, p + 1))
    want.push_back((int)p);
  CHECK(occ == want);
  puts("OK strings/prefix-function");
}
