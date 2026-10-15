#include "../testlib.hpp"
#include "../../content/strings/lyndon.hpp"

string rnd_str(int n, char hi) {
  string s(n, 'a');
  for (char& c : s) c = (char)rnd('a', hi);
  return s;
}

bool is_lyndon(const string& w) {
  rep(i, 1, sz(w))
    if (w.substr(i) <= w) return false;
  return !w.empty();
}

// Factorization properties; with uniqueness they pin it down.
void check_duval(const string& s, const vi& st) {
  int n = sz(s);
  CHECK(n == 0 ? st.empty() : !st.empty() && st[0] == 0);
  string prev;
  rep(j, 0, sz(st)) {
    int end = j + 1 < sz(st) ? st[j + 1] : n;
    CHECK(st[j] < end);
    string w = s.substr(st[j], end - st[j]);
    if (sz(w) <= 64) CHECK(is_lyndon(w));
    if (j) CHECK(w <= prev);
    prev = w;
  }
}

int brute_rotation(const string& s) {
  int best = 0;
  rep(i, 1, sz(s))
    if (s.substr(i) + s.substr(0, i) < s.substr(best) + s.substr(0, best))
      best = i;
  return best;
}

int main() {
  rep(it, 0, iters(4000)) {
    int n = (int)rnd(0, 20);
    string s = rnd_str(n, "abc"[rnd(0, 2)]);
    if (rnd(0, 3) == 0 && n) { // periodic: w^m w'
      int p = (int)rnd(1, 4);
      rep(i, p, n) s[i] = s[i - p];
    }
    vi st = duval(s);
    check_duval(s, st);
    // Usage: starts only; append n to get the end of the last one.
    vi ends = st;
    ends.push_back(n);
    string joined;
    rep(j, 0, sz(st)) joined += s.substr(ends[j], ends[j + 1] - ends[j]);
    CHECK(joined == s);
    if (n) {
      string smallest = s;
      rep(i, 0, n) smallest = min(smallest, s.substr(i));
      CHECK(s.substr(st.back()) == smallest);
    }
    CHECK_EQ(sz(st) == 1, is_lyndon(s));
    CHECK_EQ(min_rotation(s), brute_rotation(s));
    vector<ll> a(n); // big values with the same order
    rep(i, 0, n) a[i] = (s[i] - 'b') * (ll)4e18;
    CHECK(duval(a) == st);
    CHECK_EQ(min_rotation(a), min_rotation(s));
    // Tweak: largest rotation from the flipped order.
    string f = s;
    for (char& c : f) c = char('z' + 'a' - c);
    int r = min_rotation(f), best = 0;
    rep(i, 0, n)
      if (s.substr(i) + s.substr(0, i) > s.substr(best) + s.substr(0, best))
        best = i;
    CHECK(n == 0 ||
          s.substr(r) + s.substr(0, r) == s.substr(best) + s.substr(0, best));
  }
  // Example values.
  CHECK(duval(string("abbabbaba")) == vi({0, 3, 6, 8})); // no 9
  CHECK_EQ(min_rotation(string("baba")), 1);
  CHECK_EQ(min_rotation(string("")), 0);
  // Large: 1e6 characters.
  int n = 1'000'000;
  vi st = duval(string(n, 'a'));
  CHECK_EQ(sz(st), n);
  string s = rnd_str(n, 'b');
  check_duval(s, duval(s));
  rep(i, 0, n) s[i] = "bc"[rnd(0, 1)];
  int pos = (int)rnd(0, n - 1);
  rep(i, 0, 20) s[(pos + i) % n] = 'a'; // unique run of 20 a's
  CHECK_EQ(min_rotation(s), pos);
  string per;
  while (sz(per) < n) per += "aab";
  CHECK_EQ(min_rotation(per), 0);
  CHECK_EQ(min_rotation(per.substr(1) + per[0]), 2); // first of many
  puts("OK strings/lyndon");
}
