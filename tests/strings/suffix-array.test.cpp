#include "../testlib.hpp"
#include "../../content/ds/sparse-table.hpp"
#include "../../content/strings/suffix-array.hpp"

string rnd_str(int n, char hi) {
  string s(n, 'a');
  for (char& c : s) c = (char)rnd('a', hi);
  return s;
}

// Sorts suffixes by comparing them directly; LCPs by scanning.
template<class S> void check_brute(const S& s, const SuffixArray& S_) {
  int n = sz(s);
  vi sa(n), lcp(max(n - 1, 0));
  iota(all(sa), 0);
  sort(all(sa), [&](int a, int b) {
    return lexicographical_compare(begin(s) + a, end(s), begin(s) + b, end(s));
  });
  rep(i, 0, n - 1) {
    int a = sa[i], b = sa[i + 1];
    while (max(a, b) + lcp[i] < n && s[a + lcp[i]] == s[b + lcp[i]]) lcp[i]++;
  }
  CHECK(S_.sa == sa);
  CHECK(S_.lcp == lcp);
  CHECK_EQ(sz(S_.rank), n);
  rep(i, 0, n) CHECK_EQ(S_.rank[sa[i]], i);
}

// The Example: O(1) LCP of two suffixes and pattern occurrences.
void check_example(const string& s, const SuffixArray& S) {
  int n = sz(s);
  if (n < 2) return;
  SparseTable<int> st(S.lcp);
  auto lcp = [&](int a, int b) {
    if (a == b) return n - a;
    auto [i, j] = minmax(S.rank[a], S.rank[b]);
    return st.query(i, j);
  };
  rep(it, 0, 20) {
    int a = (int)rnd(0, n - 1), b = (int)rnd(0, n - 1), want = 0;
    while (max(a, b) + want < n && s[a + want] == s[b + want]) want++;
    CHECK_EQ(lcp(a, b), want);
  }
  rep(it, 0, 5) {
    int l = (int)rnd(0, n), r = (int)rnd(l, min(n, l + 4));
    string pat = it ? s.substr(l, r - l) : rnd_str((int)rnd(1, 3), 'c');
    auto lo = partition_point(
        all(S.sa), [&](int i) { return s.compare(i, sz(pat), pat) < 0; });
    auto hi = partition_point(
        lo, end(S.sa), [&](int i) { return s.compare(i, sz(pat), pat) == 0; });
    vi got(lo, hi), want;
    sort(all(got));
    rep(i, 0, n)
      if (s.compare(i, sz(pat), pat) == 0) want.push_back(i);
    CHECK(got == want);
  }
}

int main() {
  rep(it, 0, iters(3000)) {
    int n = (int)rnd(0, 30);
    string s = rnd_str(n, "abcz"[rnd(0, 3)]);
    SuffixArray S(s);
    check_brute(s, S);
    check_example(s, S);
    // Usage: distinct substrings.
    set<string> subs;
    rep(l, 0, n)
      rep(r, l + 1, n + 1) subs.insert(s.substr(l, r - l));
    ll d = (ll)n * (n + 1) / 2;
    for (int h : S.lcp) d -= h;
    CHECK_EQ(d, (ll)sz(subs));
    // Tweak: sorted cyclic shifts from s + s.
    vector<string> rot, want;
    rep(i, 0, n) want.push_back(s.substr(i) + s.substr(0, i));
    sort(all(want));
    for (int i : SuffixArray(s + s).sa)
      if (i < n) rot.push_back(s.substr(i) + s.substr(0, i));
    CHECK(rot == want);
  }
  // vector<int> with small and large lim, and every lim edge value.
  rep(it, 0, iters(1000)) {
    int n = (int)rnd(0, 25), lim = (int)rnd(2, it % 2 ? 5 : 100'000);
    vi a(n);
    for (int& x : a) x = (int)rnd(max(1, lim - 3), lim - 1);
    check_brute(a, SuffixArray(a, lim));
  }
  // All of ASCII except 0.
  string bytes;
  rep(c, 1, 128) bytes += (char)c, bytes += (char)(128 - c);
  check_brute(bytes, SuffixArray(bytes));
  CHECK(SuffixArray(string()).sa.empty());
  CHECK(SuffixArray(string("z")).sa == vi{0});
  // Large: a^n (most doubling rounds) and Fibonacci words, then a
  // random text checked by ordering and LCP properties.
  elapsed();
  int n = 100'000;
  SuffixArray A(string(n, 'a'));
  rep(i, 0, n) CHECK_EQ(A.sa[i], n - 1 - i);
  rep(i, 0, n - 1) CHECK_EQ(A.lcp[i], i + 1);
  string f = "a", g = "ab";
  while (sz(g) < n) tie(f, g) = pair(g, g + f);
  string s = rnd_str(n, 'b');
  for (const string& t : {g, s}) {
    SuffixArray S(t);
    int m = sz(t);
    rep(i, 0, m - 1) {
      int a = S.sa[i], b = S.sa[i + 1], h = S.lcp[i];
      CHECK(b + h < m && (a + h == m || t[a + h] < t[b + h]));
      if (h <= 40) CHECK(t.compare(a, h, t, b, h) == 0);
    }
  }
  fprintf(stderr, "  large inputs: %.2fs\n", elapsed());
  puts("OK strings/suffix-array");
}
