#include "../testlib.hpp"
#include "../../content/ds/dsu-rollback.hpp"

int main() {
  rep(it, 0, iters(2000)) {
    int n = (int)rnd(1, 12);
    RollbackDSU d(n);
    vi label(n);
    iota(all(label), 0);
    vector<pair<int, vi>> saved; // (checkpoint, labels)
    rep(op, 0, 40) {
      int kind = (int)rnd(0, 4);
      int a = (int)rnd(0, n - 1), b = (int)rnd(0, n - 1);
      if (kind <= 1) {
        CHECK_EQ(d.join(a, b), label[a] != label[b]);
        int from = label[b];
        for (int& x : label)
          if (x == from) x = label[a];
      } else if (kind == 2) {
        saved.push_back({d.time(), label});
      } else if (kind == 3 && !saved.empty()) {
        int k = (int)rnd(0, sz(saved) - 1);
        d.rollback(saved[k].first), label = saved[k].second;
        saved.resize(k + 1);
      }
      CHECK_EQ(d.find(a) == d.find(b), label[a] == label[b]);
      CHECK_EQ(label[d.find(a)], label[a]);
      CHECK_EQ(d.size(a), (int)count(all(label), label[a]));
      set<int> comps(all(label)); // the Tweak's component count
      CHECK_EQ(n - d.time() / 2, sz(comps));
    }
  }
  // Large: a DFS-like pattern of joins and rollbacks.
  int n = 1 << 16;
  RollbackDSU d(n);
  for (int len = 1; len < n; len *= 2) {
    int t = d.time();
    for (int i = 0; i + len < n; i += 2 * len) CHECK(d.join(i, i + len));
    CHECK_EQ(d.size(0), 2 * len);
    if (len < n / 2) {
      d.rollback(t);
      CHECK_EQ(d.size(0), len);
      for (int i = 0; i + len < n; i += 2 * len) CHECK(d.join(i, i + len));
    }
  }
  CHECK_EQ(d.size(n - 1), n);
  d.rollback(0);
  rep(i, 0, n) CHECK_EQ(d.size(i), 1);
  { // Union by size: depth <= log2 n. Joining the big set under the
    // singleton (a flipped comparison) would build a path here.
    RollbackDSU u(1 << 12);
    rep(i, 1, 1 << 12) {
      if (i % 2) u.join(i, 0);
      else u.join(0, i);
    }
    rep(x, 0, 1 << 12) {
      int depth = 0;
      for (int y = x; u.e[y] >= 0; y = u.e[y]) depth++;
      CHECK(depth <= 12);
    }
  }
  puts("OK ds/dsu-rollback");
}
