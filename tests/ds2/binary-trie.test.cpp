#include "../testlib.hpp"
#include "../../content/ds2/binary-trie.hpp"

const ll FULL = 1LL << BinaryTrie::B;

// Every query against the stored multiset ys, for one x.
void check_all(BinaryTrie& tr, const vector<ll>& ys, ll x) {
  int n = sz(ys);
  CHECK_EQ(tr.size(), n);
  vector<ll> xs, sorted_ys = ys;
  for (ll y : ys) xs.push_back(y ^ x);
  sort(all(xs)), sort(all(sorted_ys));
  rep(k, 0, n) {
    CHECK_EQ(tr.kth(k, x), xs[k]);
    CHECK_EQ(tr.kth(k), sorted_ys[k]);
  }
  if (n) {
    CHECK_EQ(tr.min_xor(x), xs[0]);
    CHECK_EQ(tr.max_xor(x), xs[n - 1]);
  }
  for (ll k : {0LL, 1LL, x, FULL - 1, FULL, rnd(0, FULL)}) {
    int want = (int)count_if(all(xs), [&](ll v) { return v < k; });
    CHECK_EQ(tr.count_less(x, k), want);
  }
  CHECK_EQ(tr.count_less(x, 1), (int)count(all(ys), x));
}

void stress(int rounds, ll hi) {
  rep(it, 0, iters(rounds)) {
    BinaryTrie tr;
    vector<ll> ys;
    rep(op, 0, 30) {
      int kind = (int)rnd(0, 3);
      if (kind <= 1 || ys.empty()) {
        ll y = rnd(0, hi);
        int c = kind == 0 ? 1 : (int)rnd(1, 3);
        tr.add(y, c), ys.insert(ys.end(), c, y);
      } else {
        int i = (int)rnd(0, sz(ys) - 1);
        tr.add(ys[i], -1), ys.erase(ys.begin() + i);
      }
      check_all(tr, ys, rnd(0, hi));
    }
    for (ll y : ys) tr.add(y, -1); // empty again, nodes stay
    check_all(tr, {}, rnd(0, hi));
    tr.add(hi, 1);
    check_all(tr, {hi}, rnd(0, hi));
  }
}

int main() {
  stress(1500, 15);
  stress(1000, FULL - 1);
  // Extremes: 0 and 2^B - 1 against x = 0 and x = 2^B - 1.
  BinaryTrie tr;
  tr.add(0, 1), tr.add(FULL - 1, 1);
  CHECK_EQ(tr.max_xor(0), FULL - 1);
  CHECK_EQ(tr.min_xor(FULL - 1), 0LL);
  CHECK_EQ(tr.count_less(0, FULL), 2);
  CHECK_EQ(tr.count_less(0, FULL - 1), 1);
  // Large: 1e5 random values, sampled queries against brute force.
  elapsed();
  int n = 100'000;
  BinaryTrie big;
  vector<ll> ys = rnd_vec(n, 0, FULL - 1);
  for (ll y : ys) big.add(y, 1);
  rep(it, 0, 100) {
    ll x = rnd(0, FULL - 1), k = rnd(0, FULL), mn = FULL, mx = -1;
    int less = 0;
    for (ll y : ys)
      mn = min(mn, x ^ y), mx = max(mx, x ^ y), less += (x ^ y) < k;
    CHECK_EQ(big.min_xor(x), mn);
    CHECK_EQ(big.max_xor(x), mx);
    CHECK_EQ(big.count_less(x, k), less);
  }
  sort(all(ys));
  rep(k, 0, n) CHECK_EQ(big.kth(k), ys[k]);
  fprintf(stderr, "  large inputs: %.2fs\n", elapsed());
  puts("OK ds2/binary-trie");
}
