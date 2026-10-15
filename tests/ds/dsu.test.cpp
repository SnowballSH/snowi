#include "../testlib.hpp"
#include "../../content/ds/dsu.hpp"

int main() {
  rep(it, 0, iters(3000)) {
    int n = (int)rnd(1, 15);
    DSU d(n);
    vi label(n);
    iota(all(label), 0);
    rep(op, 0, 30) {
      int a = (int)rnd(0, n - 1), b = (int)rnd(0, n - 1);
      if (rnd(0, 1)) {
        bool merged = label[a] != label[b];
        CHECK_EQ(d.join(a, b), merged);
        int from = label[b];
        for (int& x : label)
          if (x == from) x = label[a];
      }
      CHECK_EQ(d.same(a, b), label[a] == label[b]);
      CHECK_EQ(label[d.find(a)], label[a]);
      CHECK_EQ(d.size(a), (int)count(all(label), label[a]));
    }
  }
  // Large: merge by doubling (deepest trees), then a long chain.
  int n = 1 << 16;
  DSU d(n);
  for (int len = 1; len < n; len *= 2)
    for (int i = 0; i + len < n; i += 2 * len) CHECK(d.join(i + len, i));
  CHECK_EQ(d.size(n - 1), n);
  rep(i, 0, n) CHECK_EQ(d.find(i), d.find(0));
  DSU c(n);
  rep(i, 1, n) CHECK(c.join(i - 1, i));
  rep(i, 0, n) CHECK(!c.join(0, i));
  CHECK_EQ(c.size(0), n);
  puts("OK ds/dsu");
}
