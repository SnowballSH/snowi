#include "../testlib.hpp"
#include "../../content/graph/floyd-warshall.hpp"

typedef vector<vector<ll>> Mat;
const ll NONE = LLONG_MAX;

// walk[k][x][v]: lightest walk with exactly k edges from x to v.
Mat brute(const Mat& w) {
  int n = sz(w);
  vector<Mat> walk(n + 1, Mat(n, vector<ll>(n, NONE)));
  rep(x, 0, n) walk[0][x][x] = 0;
  rep(k, 0, n)
    rep(x, 0, n)
      rep(a, 0, n)
        rep(b, 0, n)
          if (walk[k][x][a] != NONE && w[a][b] != NONE)
            walk[k + 1][x][b] = min(walk[k + 1][x][b], walk[k][x][a] + w[a][b]);
  auto reach = [&](int x, int v) {
    rep(k, 0, n)
      if (walk[k][x][v] != NONE) return true;
    return false;
  };
  vector<bool> on_neg(n);
  rep(x, 0, n)
    rep(k, 1, n + 1)
      if (walk[k][x][x] != NONE) on_neg[x] = on_neg[x] || walk[k][x][x] < 0;
  Mat d(n, vector<ll>(n, NONE));
  rep(i, 0, n)
    rep(j, 0, n) {
      rep(k, 0, n) d[i][j] = min(d[i][j], walk[k][i][j]);
      rep(y, 0, n)
        if (on_neg[y] && reach(i, y) && reach(y, j)) d[i][j] = LLONG_MIN;
    }
  return d;
}

// The Example's O(n^2) update for a new edge a -> b of weight w.
void add_edge(Mat& d, int a, int b, ll w) {
  int n = sz(d);
  rep(i, 0, n)
    rep(j, 0, n)
      if (d[i][a] != NONE && d[b][j] != NONE)
        d[i][j] = min(d[i][j], d[i][a] + w + d[b][j]);
}

int main() {
  int finite = 0, neg = 0, added = 0;
  rep(it, 0, iters(2000)) {
    int n = (int)rnd(1, 6);
    ll lo = rnd(0, 2) ? -4 : 0, hi = 10;
    if (rnd(0, 4) == 0) lo = -1'000'000'000'000, hi = -lo;
    Mat w(n, vector<ll>(n, NONE));
    int m = (int)rnd(0, 12);
    rep(e, 0, m) {
      ll& x = w[rnd(0, n - 1)][rnd(0, n - 1)];
      x = min(x, rnd(lo, hi));
    }
    Mat d = w;
    floyd_warshall(d);
    Mat want = brute(w);
    CHECK(d == want);
    bool cyc = false;
    for (auto& row : want)
      for (ll x : row)
        finite += x != NONE && x != LLONG_MIN, neg += x == LLONG_MIN,
            cyc |= x == LLONG_MIN;
    // A new edge that closes no negative cycle: O(n^2) update.
    int a = (int)rnd(0, n - 1), b = (int)rnd(0, n - 1);
    ll x = rnd(lo, hi);
    if (cyc || (d[b][a] != NONE && d[b][a] + x < 0)) continue;
    add_edge(d, a, b, x), added++;
    w[a][b] = min(w[a][b], x);
    floyd_warshall(w);
    CHECK(d == w);
  }
  CHECK(min(finite, neg) > (finite + neg) / 10);
  CHECK(added > iters(2000) / 4);
  // Unreachable pairs stay unreachable (no overflow, no INF - 1).
  for (ll x : {4LL, -1LL}) {
    Mat d(3, vector<ll>(3, NONE));
    d[0][1] = 5;
    floyd_warshall(d);
    add_edge(d, 1, 2, x);
    Mat want = {{0, 5, 5 + x}, {NONE, 0, x}, {NONE, NONE, 0}};
    CHECK(d == want);
  }
  // Many heavy negative cycles: values would overflow without the
  // clamp (UBSan aborts on signed overflow).
  int n = 60;
  Mat d(n, vector<ll>(n, -1'000'000'000'000'000LL));
  floyd_warshall(d);
  CHECK(d == Mat(n, vector<ll>(n, LLONG_MIN)));
  // Large, positive weights: certify against the edges directly.
  n = 200;
  Mat w(n, vector<ll>(n, NONE));
  rep(e, 0, 3000) w[rnd(0, n - 1)][rnd(0, n - 1)] = rnd(1, 1'000'000'000);
  d = w;
  floyd_warshall(d);
  rep(i, 0, n)
    rep(j, 0, n) {
      if (i == j) {
        CHECK_EQ(d[i][j], 0LL);
        continue;
      }
      bool tight = false;
      rep(k, 0, n)
        if (w[k][j] != NONE && d[i][k] != NONE) {
          CHECK(d[i][j] <= d[i][k] + w[k][j]);
          tight |= d[i][j] == d[i][k] + w[k][j];
        }
      CHECK(tight || d[i][j] == NONE);
    }
  // Weights at the bound n * max|w| < 2^62: restricted paths can be
  // much longer than every distance, but no sum may overflow (UBSan).
  {
    int nb = 6;
    ll big = ((1LL << 62) - 1) / nb;
    Mat wb(nb, vector<ll>(nb, NONE));
    for (auto [a, b] : vector<pii>{{3, 0}, {0, 1}, {1, 2}, {2, 0}, {1, 4}})
      wb[a][b] = big;
    rep(v, 0, 5) wb[v][5] = wb[5][v] = 0;
    Mat db = wb;
    floyd_warshall(db);
    CHECK(db == brute(wb));
  }
  rep(it, 0, iters(500)) {
    int nb = (int)rnd(2, 6);
    ll big = ((1LL << 62) - 1) / nb;
    Mat wb(nb, vector<ll>(nb, NONE));
    rep(e, 0, (int)rnd(1, 12)) {
      ll x = rnd(0, 1) ? rnd(big - 9, big) : rnd(-big, big);
      ll& y = wb[rnd(0, nb - 1)][rnd(0, nb - 1)];
      y = min(y, x);
    }
    Mat db = wb;
    floyd_warshall(db);
    CHECK(db == brute(wb));
  }
  puts("OK graph/floyd-warshall");
}
