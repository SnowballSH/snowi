#include "../testlib.hpp"
#include "../../content/ds/fenwick-2d.hpp"

int main() {
  rep(it, 0, iters(2000)) {
    int n = (int)rnd(0, 8), m = (int)rnd(0, 8);
    vector<vector<ll>> a(n, vector<ll>(m));
    Fenwick2D fw(n, m);
    ll big = rnd(0, 3) ? 10 : ll(1e15);
    rep(op, 0, 30) {
      if (n && m) {
        int x = (int)rnd(0, n - 1), y = (int)rnd(0, m - 1);
        ll v = rnd(-big, big);
        a[x][y] += v, fw.add(x, y, v);
      }
      int x1 = (int)rnd(0, n), x2 = (int)rnd(x1, n);
      int y1 = (int)rnd(0, m), y2 = (int)rnd(y1, m);
      ll want = 0;
      rep(i, x1, x2)
        rep(j, y1, y2) want += a[i][j];
      CHECK_EQ(fw.rect(x1, y1, x2, y2), want);
      want = 0;
      rep(i, 0, x2)
        rep(j, 0, y2) want += a[i][j];
      CHECK_EQ(fw.query(x2, y2), want);
    }
  }
  // The Tweak's rectangle add: corners at x2 = n or y2 = m are no-ops.
  rep(it, 0, iters(1000)) {
    int n = (int)rnd(0, 8), m = (int)rnd(0, 8);
    vector<vector<ll>> a(n, vector<ll>(m));
    Fenwick2D fw(n, m);
    rep(op, 0, 20) {
      int x1 = (int)rnd(0, n), x2 = (int)rnd(x1, n);
      int y1 = (int)rnd(0, m), y2 = (int)rnd(y1, m);
      ll v = rnd(-10, 10);
      fw.add(x1, y1, v), fw.add(x2, y2, v);
      fw.add(x1, y2, -v), fw.add(x2, y1, -v);
      rep(i, x1, x2)
        rep(j, y1, y2) a[i][j] += v;
      rep(i, 0, n)
        rep(j, 0, m) CHECK_EQ(fw.query(i + 1, j + 1), a[i][j]);
    }
  }
  // Large: 1e5 adds on a 512 x 300 grid, then 1e5 rectangle sums
  // checked against 2D prefix sums.
  int n = 512, m = 300, q = 100'000;
  Fenwick2D fw(n, m);
  vector<vector<ll>> pre(n + 1, vector<ll>(m + 1));
  rep(it, 0, q) {
    int x = (int)rnd(0, n - 1), y = (int)rnd(0, m - 1);
    ll v = rnd(-1e9, 1e9);
    fw.add(x, y, v), pre[x + 1][y + 1] += v;
  }
  rep(i, 0, n)
    rep(j, 0, m) pre[i + 1][j + 1] += pre[i][j + 1] + pre[i + 1][j] - pre[i][j];
  rep(it, 0, q) {
    int x1 = (int)rnd(0, n), x2 = (int)rnd(x1, n);
    int y1 = (int)rnd(0, m), y2 = (int)rnd(y1, m);
    ll want = pre[x2][y2] - pre[x1][y2] - pre[x2][y1] + pre[x1][y1];
    CHECK_EQ(fw.rect(x1, y1, x2, y2), want);
  }
  puts("OK ds/fenwick-2d");
}
