#include "../testlib.hpp"
#include "../../content/ds2/hash-map.hpp"

// Random operations on both containers against map.
template<class M> void stress(ll lo, ll hi) {
  rep(it, 0, iters(300)) {
    M h;
    map<ll, int> ref;
    rep(op, 0, 200) {
      ll x = rnd(lo, hi);
      int kind = (int)rnd(0, 2);
      if (kind == 0) h[x] += 1, ref[x] += 1;
      else if (kind == 1) CHECK_EQ(h.erase(x) > 0, ref.erase(x) > 0);
      else {
        auto pos = h.find(x);
        CHECK_EQ(pos != h.end(), ref.count(x) > 0);
        if (pos != h.end()) CHECK_EQ(pos->second, ref[x]);
      }
      CHECK_EQ(sz(h), sz(ref));
    }
    map<ll, int> back(h.begin(), h.end());
    CHECK(back == ref);
  }
}

// Distinct values of the low (or high) 16 bits over 2^16 keys.
// Random hashes give about 1 - 1/e = 63% distinct.
int spread(const vector<uint64_t>& keys, int shift) {
  vector<bool> seen(1 << 16);
  int d = 0;
  for (auto k : keys) {
    size_t b = SafeHash{}(k) >> shift & 0xffff;
    d += !seen[b], seen[b] = 1;
  }
  return d;
}

// Pairs of ll (and of int) against map, extreme coordinates.
template<class M> void stress_pairs() {
  rep(it, 0, iters(100)) {
    M h;
    map<pair<ll, ll>, int> ref;
    rep(op, 0, 300) {
      ll a = rnd(0, 3), b = rnd(-2, 2);
      if (rnd(0, 1)) a = a << 40 | rnd(0, 1); // high bits differ
      else a = rnd(LLONG_MIN, LLONG_MAX);
      if (rnd(0, 1)) b = rnd(LLONG_MIN, LLONG_MAX);
      h[{a, b}]++, ref[{a, b}]++;
      ll x = rnd(0, 3) << 40, y = rnd(-2, 2);
      auto pos = h.find({x, y});
      CHECK_EQ(pos != h.end(), ref.count({x, y}) > 0);
    }
    CHECK_EQ(sz(h), sz(ref));
    map<pair<ll, ll>, int> back(h.begin(), h.end());
    CHECK(back == ref);
  }
}

// Regression: keys built to defeat a pair hash that mixes the
// first component without the seed, as (*this)(mix(a) ^ b) does:
// (a, mix(a) ^ C) all reach the seeded hash as C. And pair<ll, ll>
// must not be truncated to pii (keys differing in high bits).
void check_pair_attacks() {
  SafeHash h;
  int n = 20000;
  set<size_t> crafted, high;
  rep(a, 0, n) {
    ll b = ll(SafeHash::mix(a) ^ 12345);
    crafted.insert(h(pair<ll, ll>(a, b)));
    high.insert(h(pair<ll, ll>((ll)a << 32, 7)));
  }
  CHECK(sz(crafted) > n - 10);
  CHECK(sz(high) > n - 10);
  CHECK(h(pair<ll, ll>(1LL << 40, 3)) != h(pii(0, 3)));
  // Nor a linear combiner: (a, C - K * a) all collide under
  // a * K + b.
  for (ll K : {31LL, 1'000'003LL}) {
    set<size_t> lin;
    rep(a, 0, n) lin.insert(h(pair<ll, ll>(a, 12345 - K * a)));
    CHECK(sz(lin) > n - 10);
  }
}

int main() {
  stress_pairs<unordered_map<pair<ll, ll>, int, SafeHash>>();
  stress_pairs<gp_hash_table<pair<ll, ll>, int, SafeHash>>();
  check_pair_attacks();
  stress<unordered_map<ll, int, SafeHash>>(-20, 20);
  stress<gp_hash_table<ll, int, SafeHash>>(-20, 20);
  stress<unordered_map<ll, int, SafeHash>>(LLONG_MIN, LLONG_MAX);
  stress<gp_hash_table<ll, int, SafeHash>>(LLONG_MIN, LLONG_MAX);
  // Reference outputs of the splitmix64 generator from state 0.
  CHECK_EQ(SafeHash::mix(0), 0xe220a8397b1dcdafULL);
  CHECK_EQ(SafeHash::mix(0x9e3779b97f4a7c15), 0x6e789e6aa1b965f4ULL);
  // Same key, same hash; the seed is applied; pair keys work.
  SafeHash h;
  rep(i, 0, 100) {
    ll x = rnd(LLONG_MIN, LLONG_MAX);
    CHECK_EQ(h(x), h(x));
    CHECK(h(x) != SafeHash::mix(x));
  }
  unordered_map<pii, int, SafeHash> grid;
  map<pii, int> ref;
  rep(i, 0, 20000) {
    pii p(rnd(-30, 30), rnd(INT_MIN, INT_MIN + 30));
    grid[p]++, ref[p]++;
  }
  map<pii, int> back(all(grid));
  CHECK(back == ref);
  // Structured keys spread over the low and the high bits.
  for (int s : {0, 20, 32, 47}) {
    vector<uint64_t> keys(1 << 16);
    rep(i, 0, 1 << 16) keys[i] = (uint64_t)i << s;
    CHECK(spread(keys, 0) > 38000);
    CHECK(spread(keys, 48) > 38000);
  }
  // Large: keys x << 20, which make a masked hash collide.
  elapsed();
  gp_hash_table<ll, int, SafeHash> g({}, {}, {}, {}, {1 << 18});
  rep(i, 0, 200'000) g[(ll)i << 20] = i;
  rep(i, 0, 200'000) CHECK_EQ(g[(ll)i << 20], i);
  CHECK_EQ(sz(g), 200'000);
  fprintf(stderr, "  large inputs: %.2fs\n", elapsed());
  puts("OK ds2/hash-map");
}
