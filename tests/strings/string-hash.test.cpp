#include "../testlib.hpp"
#include "../../content/strings/string-hash.hpp"

typedef uint64_t ull;
const ull M = (1ULL << 61) - 1;

// Independent evaluation: Horner with a full 128-bit remainder.
template<class S> ull eval(const S& s, int l, int r) {
  ull h = 0;
  rep(i, l, r) h = (ull)(((__uint128_t)h * StrHash::B + s[i]) % M);
  return h;
}

string rnd_str(int n, char hi) {
  string s(n, 'a');
  for (char& c : s) c = (char)rnd('a', hi);
  return s;
}

// The header's Example: LCP of s[a, n) and s[b, n), then compare.
pair<int, bool> lcp_less(StrHash& hs, const string& s, int a, int b) {
  int n = sz(s);
  int lo = 0, hi = n - max(a, b);
  while (lo < hi) {
    int mid = (lo + hi + 1) / 2;
    if (hs.get(a, a + mid) == hs.get(b, b + mid)) lo = mid;
    else hi = mid - 1;
  }
  return {lo, lo < n - b && (lo == n - a || s[a + lo] < s[b + lo])};
}

int main() {
  CHECK(StrHash::B < M);
  rep(it, 0, iters(1500)) {
    int n = (int)rnd(0, 16);
    string s = rnd_str(n, "abc"[rnd(0, 2)]);
    StrHash hs(s);
    rep(l, 0, n + 1)
      rep(r, l, n + 1) {
        CHECK_EQ(hs.get(l, r), eval(s, l, r));
        // equal hashes exactly for equal substrings
        rep(l2, 0, n - (r - l) + 1)
          CHECK_EQ(hs.get(l, r) == hs.get(l2, l2 + r - l),
                   s.compare(l, r - l, s, l2, r - l) == 0);
      }
    rep(a, 0, n)
      rep(b, 0, n) {
        auto [lcp, less] = lcp_less(hs, s, a, b);
        int want = 0;
        while (max(a, b) + want < n && s[a + want] == s[b + want]) want++;
        CHECK_EQ(lcp, want);
        CHECK_EQ(less, s.substr(a) < s.substr(b));
      }
    // Tweak: hash of a concatenation; another string, same base.
    string t = rnd_str((int)rnd(0, 8), 'c');
    StrHash ht(t), hst(s + t);
    rep(l, 0, n + 1)
      rep(r, 0, sz(t) + 1) {
        ull x = hs.get(l, n), y = ht.get(0, r);
        CHECK_EQ(ht.norm(ht.mul(x, ht.pw[r]) + y), hst.get(l, n + r));
      }
  }
  // Values up to M - 1, zeros, and sums that hit M exactly.
  rep(it, 0, iters(1500)) {
    int n = (int)rnd(1, 12);
    vector<ll> a(n);
    for (ll& x : a) x = rnd(0, 3) ? rnd(0, 2) : (ll)rnd(M - 3, M - 1);
    ull h = 0;
    rep(i, 0, n) {
      if (h && rnd(0, 2) == 0) a[i] = (ll)(M - StrHash::mul(h, StrHash::B));
      h = eval(a, 0, i + 1);
    }
    StrHash ha(a);
    rep(l, 0, n + 1)
      rep(r, l, n + 1) CHECK_EQ(ha.get(l, r), eval(a, l, r));
  }
  StrHash empty(string{});
  CHECK_EQ(empty.get(0, 0), 0ULL);
  // Notes: Thue-Morse breaks overflow hashing mod 2^64, not this one.
  int len = 1 << 11;
  string tm(len, 'a'), co(len, 'b');
  rep(i, 0, len)
    if (__builtin_popcount(i) & 1) tm[i] = 'b', co[i] = 'a';
  ull base = rnd(0, LLONG_MAX) | 1, h1 = 0, h2 = 0; // any odd base
  rep(i, 0, len) h1 = h1 * base + tm[i], h2 = h2 * base + co[i];
  CHECK_EQ(h1, h2);
  CHECK(StrHash(tm).get(0, len) != StrHash(co).get(0, len));
  // Large: 1e5 random equal-length comparisons in a periodic string
  // with a few changed characters, so many pairs are equal.
  int n = 200'000;
  string s(n, 'a');
  rep(i, 0, n) s[i] = "abcab"[i % 5];
  rep(k, 0, 20) s[rnd(0, n - 1)] = 'c';
  StrHash hs(s);
  int equal = 0;
  rep(it, 0, 100'000) {
    int len2 = (int)rnd(1, 2000);
    int a = (int)rnd(0, n - len2), b = (int)rnd(0, (n - len2) / 5) * 5;
    if (rnd(0, 1)) b = min(n - len2, b + a % 5);
    bool same = s.compare(a, len2, s, b, len2) == 0;
    CHECK_EQ(hs.get(a, a + len2) == hs.get(b, b + len2), same);
    equal += same;
  }
  CHECK(equal > 10'000);
  puts("OK strings/string-hash");
}
