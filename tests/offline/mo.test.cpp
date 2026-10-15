#include "../testlib.hpp"
#include "../../content/offline/mo.hpp"

struct Result {
  vi distinct;
  vector<ll> sumsq;
  ll moves = 0;
};

// Runs mo with callbacks that check the window: every add or
// remove happens at an end, and answer(q) sees exactly qs[q]. An
// empty window has no position: removing its last element from
// either end looks the same.
Result run(const vi& a, const vector<pii>& qs, int V, bool hil) {
  int n = sz(a), lo = 0, hi = 0, distinct = 0;
  vi cnt(V);
  ll sq = 0;
  Result res{vi(sz(qs), -1), vector<ll>(sz(qs)), 0};
  auto add = [&](int i) {
    CHECK(0 <= i && i < n);
    if (lo == hi) lo = i, hi = i + 1;
    else {
      CHECK(i == lo - 1 || i == hi);
      i == hi ? hi++ : lo--;
    }
    distinct += !cnt[a[i]]++;
    sq += 2 * cnt[a[i]] - 1;
    res.moves++;
  };
  auto remove = [&](int i) {
    CHECK(lo < hi && (i == lo || i == hi - 1));
    i == lo ? lo++ : hi--;
    distinct -= !--cnt[a[i]];
    sq -= 2 * cnt[a[i]] + 1;
    res.moves++;
  };
  auto answer = [&](int q) {
    auto [l, r] = qs[q];
    CHECK(lo == hi ? l == r : lo == l && hi == r);
    CHECK_EQ(res.distinct[q], -1); // answered once
    res.distinct[q] = distinct, res.sumsq[q] = sq;
  };
  mo(n, qs, add, remove, answer, hil);
  return res;
}

void stress(int rounds, int maxn, int maxv) {
  rep(it, 0, iters(rounds)) {
    int n = (int)rnd(0, maxn), q = (int)rnd(0, 30);
    vi a(n);
    for (int& x : a) x = (int)rnd(0, maxv - 1);
    vector<pii> qs(q);
    for (auto& [l, r] : qs) {
      l = (int)rnd(0, n), r = (int)rnd(0, n);
      if (l > r) swap(l, r);
    }
    for (bool hil : {false, true}) {
      Result res = run(a, qs, maxv, hil);
      rep(j, 0, q) {
        map<int, ll> c;
        rep(i, qs[j].first, qs[j].second) c[a[i]]++;
        ll sq = 0;
        for (auto [v, k] : c) sq += k * k;
        CHECK_EQ(res.distinct[j], sz(c));
        CHECK_EQ(res.sumsq[j], sq);
      }
    }
  }
}

// The first 4^k cells of the curve fill the 2^k square at the
// origin, and consecutive cells are neighbours.
void check_hilbert(int k) {
  int side = 1 << k;
  vector<pii> at(side * side, {-1, -1});
  rep(x, 0, side)
    rep(y, 0, side) {
      ll d = hilbert(x, y);
      CHECK(d < side * side);
      CHECK(at[d].first < 0);
      at[d] = {x, y};
    }
  rep(d, 1, side * side) {
    auto [x0, y0] = at[d - 1];
    auto [x1, y1] = at[d];
    CHECK_EQ(abs(x0 - x1) + abs(y0 - y1), 1);
  }
}

int main() {
  check_hilbert(6);
  CHECK_EQ(hilbert((1 << 21) - 1, 0), (1LL << 42) - 1); // last cell
  stress(2000, 12, 3);
  stress(1000, 40, 10);
  stress(300, 300, 1000);
  { // n = 0, empty queries, no queries
    Result res = run({}, {{0, 0}, {0, 0}}, 1, false);
    CHECK(res.distinct == vi({0, 0}));
    CHECK_EQ(res.moves, 0LL);
    CHECK(run({5}, {}, 6, true).distinct.empty());
  }
  // Large: n = q = 1e5, both orders; the number of moves must stay
  // near n sqrt(q) (6.3e7 if the sort were wrong, e.g. by l only).
  elapsed();
  int n = 100'000, q = 100'000, V = 1000;
  vi a(n);
  for (int& x : a) x = (int)rnd(0, V - 1);
  vector<pii> qs(q);
  for (auto& [l, r] : qs) l = (int)rnd(0, n - 1), r = (int)rnd(l + 1, n);
  Result block = run(a, qs, V, false), hil = run(a, qs, V, true);
  CHECK(block.distinct == hil.distinct);
  CHECK(block.sumsq == hil.sumsq);
  double bound = 1.2 * n * sqrt(q);
  fprintf(stderr, "  moves / (n sqrt q): blocks %.2f, Hilbert %.2f\n",
          block.moves / (n * sqrt(q)), hil.moves / (n * sqrt(q)));
  CHECK(block.moves < bound);
  CHECK(hil.moves < bound);
  for (int j = 0; j < q; j += 997) {
    auto [l, r] = qs[j];
    vi seen(V);
    int d = 0;
    rep(i, l, r) d += !seen[a[i]]++;
    CHECK_EQ(block.distinct[j], d);
  }
  fprintf(stderr, "  large inputs: %.2fs\n", elapsed());
  puts("OK offline/mo");
}
