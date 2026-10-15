#include "../testlib.hpp"
#include "../../content/ds/li-chao.hpp"

struct Seg {
  ll k, m, l, r; // k x + m for x in [l, r)
};

ll brute(const vector<Seg>& segs, ll x) {
  ll y = LiChao::INF;
  for (auto& s : segs)
    if (s.l <= x && x < s.r) y = min(y, s.k * x + s.m);
  return y;
}

// xs sorted, possibly with duplicates; values drawn from [-c, c].
vector<ll> rnd_xs(int n, ll c, bool distinct) {
  vector<ll> xs = rnd_vec(n, -c, c);
  sort(all(xs));
  if (distinct) xs.erase(unique(all(xs)), xs.end());
  return xs;
}

void stress(int rounds, ll cx, ll ck, ll cm) {
  rep(it, 0, iters(rounds)) {
    vector<ll> xs = rnd_xs((int)rnd(1, 14), cx, rnd(0, 3) > 0);
    LiChao lc(xs);
    vector<Seg> segs;
    rep(op, 0, 40) {
      ll k = rnd(-ck, ck), m = rnd(-cm, cm);
      if (rnd(0, 1)) {
        segs.push_back({k, m, LLONG_MIN, LLONG_MAX});
        lc.add_line(k, m);
      } else {
        ll l = rnd(-cx - 2, cx + 2), r = rnd(-cx - 2, cx + 2);
        if (l > r && rnd(0, 3)) swap(l, r); // else empty
        segs.push_back({k, m, l, r});
        lc.add_segment(k, m, l, r);
      }
      ll x = xs[rnd(0, sz(xs) - 1)];
      CHECK_EQ(lc.query(x), brute(segs, x));
    }
    for (ll x : xs) CHECK_EQ(lc.query(x), brute(segs, x));
  }
}

// The Notes' rollback on trees relies on add_line changing at most
// one node per level: count the changed (K, M) entries.
void check_changes_per_line() {
  rep(it, 0, iters(300)) {
    vector<ll> xs = rnd_xs((int)rnd(1, 40), 50, true);
    LiChao lc(xs);
    int levels = __lg(lc.size) + 1;
    rep(op, 0, 30) {
      vector<ll> k0 = lc.K, m0 = lc.M;
      lc.add_line(rnd(-5, 5), rnd(-100, 100));
      int changed = 0;
      rep(v, 0, sz(k0)) changed += k0[v] != lc.K[v] || m0[v] != lc.M[v];
      CHECK(changed <= levels);
    }
  }
}

int main() {
  check_changes_per_line();
  // Small values: many ties, equal slopes, empty points.
  stress(3000, 8, 3, 10);
  stress(1000, 30, 10, 100);
  // Extremes: |k x| <= 1e18 and |m| <= 1e18 (UBSan checks overflow).
  stress(1000, 1'000'000'000, 1'000'000'000, 1'000'000'000'000'000'000);
  { // one point, no line yet
    LiChao lc({7});
    CHECK_EQ(lc.query(7), LiChao::INF);
    lc.add_segment(1, 1, 8, 9);
    CHECK_EQ(lc.query(7), LiChao::INF);
    lc.add_segment(2, -3, 7, 8);
    CHECK_EQ(lc.query(7), 11LL);
  }
  // Large: tangents -2 t x + t^2 of -x^2 in random order give -x^2
  // at every x in the set; then random segments against brute force.
  elapsed();
  int n = 100'000;
  vector<ll> xs(n), ts(n);
  iota(all(xs), 0), iota(all(ts), 0);
  shuffle(all(ts), test_rng);
  LiChao lc(xs);
  for (ll t : ts) lc.add_line(-2 * t, t * t);
  rep(x, 0, n) CHECK_EQ(lc.query(x), -(ll)x * x);
  LiChao seg(xs);
  vector<Seg> segs;
  rep(i, 0, n) {
    ll l = rnd(0, n - 1), r = rnd(l + 1, n);
    segs.push_back({rnd(-1e6, 1e6), rnd(-1e12, 1e12), l, r});
    seg.add_segment(segs[i].k, segs[i].m, l, r);
  }
  rep(i, 0, 200) {
    ll x = rnd(0, n - 1);
    CHECK_EQ(seg.query(x), brute(segs, x));
  }
  fprintf(stderr, "  large inputs: %.2fs\n", elapsed());
  puts("OK ds/li-chao");
}
