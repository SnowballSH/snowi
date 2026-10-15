#include "../testlib.hpp"
#include "../../content/ds/line-container.hpp"

ll brute(const vector<pair<ll, ll>>& lines, ll x) {
  ll y = LLONG_MIN;
  for (auto [k, m] : lines) y = max(y, k * x + m);
  return y;
}

void stress(int rounds, ll ck, ll cm, ll cx) {
  rep(it, 0, iters(rounds)) {
    LineContainer lc;
    vector<pair<ll, ll>> lines;
    int n = (int)rnd(1, 30);
    rep(i, 0, n) {
      ll k = rnd(-ck, ck), m = rnd(-cm, cm);
      lines.push_back({k, m}), lc.add(k, m);
      rep(q, 0, 3) {
        ll x = rnd(-cx, cx);
        CHECK_EQ(lc.query(x), brute(lines, x));
      }
    }
    // Breakpoints strictly increase (no line is kept for an empty
    // range) and each line is a maximum at its own p.
    ll last = LLONG_MIN;
    for (auto& l : lc) {
      CHECK(last < l.p);
      last = l.p;
      if (abs(l.p) > cx) continue;
      CHECK_EQ(l.k * l.p + l.m, brute(lines, l.p));
      for (ll x : {l.p - 1, l.p + 1})
        if (abs(x) <= cx) CHECK_EQ(lc.query(x), brute(lines, x));
    }
    CHECK_EQ(last, LineContainer::inf);
  }
}

int main() {
  stress(3000, 4, 10, 12); // equal slopes, ties, negative floors
  stress(1500, 50, 1000, 100);
  stress(800, 1'000'000'000, 1'000'000'000'000'000'000, 1'000'000'000);
  stress(300, (1LL << 62) - 1, (1LL << 62) - 1, 1);
  { // minimum by negation (the Tweak)
    LineContainer lc;
    for (auto [k, m] : {pair<ll, ll>{1, 0}, {-1, 0}, {0, 2}}) lc.add(-k, -m);
    rep(x, -5, 6) CHECK_EQ(-lc.query(x), min({(ll)x, (ll)-x, 2LL}));
  }
  // Large: tangents 2 t x - t^2 of x^2 in random order give x^2.
  elapsed();
  int n = 100'000;
  vector<ll> ts(n);
  iota(all(ts), 0);
  shuffle(all(ts), test_rng);
  LineContainer lc;
  for (ll t : ts) lc.add(2 * t, -t * t);
  CHECK_EQ(sz(lc), n);
  rep(x, 0, n) CHECK_EQ(lc.query(x), (ll)x * x);
  LineContainer rnd_lc;
  vector<pair<ll, ll>> lines;
  rep(i, 0, n) {
    lines.push_back({rnd(-1e9, 1e9), rnd(-1e18, 1e18)});
    rnd_lc.add(lines[i].first, lines[i].second);
  }
  rep(i, 0, 300) {
    ll x = rnd(-1e9, 1e9);
    CHECK_EQ(rnd_lc.query(x), brute(lines, x));
  }
  fprintf(stderr, "  large inputs: %.2fs\n", elapsed());
  puts("OK ds/line-container");
}
