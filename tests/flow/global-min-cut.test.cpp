#include "../testlib.hpp"
#include "../../content/flow/global-min-cut.hpp"
#include "../../content/flow/dinic.hpp"

typedef vector<vector<ll>> Mat;

ll cut_weight(const Mat& w, const vector<bool>& in) {
  ll c = 0;
  rep(u, 0, sz(w))
    rep(v, 0, sz(w))
      if (in[u] && !in[v]) c += w[u][v];
  return c;
}

// Minimum over every non-empty proper subset.
ll brute(const Mat& w) {
  int n = sz(w);
  ll best = LLONG_MAX;
  rep(mask, 1, (1 << n) - 1) {
    vector<bool> in(n);
    rep(v, 0, n) in[v] = mask >> v & 1;
    best = min(best, cut_weight(w, in));
  }
  return best;
}

// Minimum over t of the max flow from 0 to t.
ll by_flows(const Mat& w) {
  int n = sz(w);
  ll best = LLONG_MAX;
  rep(t, 1, n) {
    Dinic d(n);
    rep(u, 0, n)
      rep(v, 0, u)
        if (w[u][v]) d.add_edge(u, v, w[u][v], w[u][v]);
    best = min(best, d.max_flow(0, t));
  }
  return best;
}

// Returns the cut after checking that side is a proper part
// whose weight is the cut.
ll solve(const Mat& w) {
  int n = sz(w);
  auto [cut, side] = global_min_cut(w);
  CHECK(0 < sz(side) && sz(side) < n);
  vector<bool> in(n);
  for (int v : side) {
    CHECK(0 <= v && v < n && !in[v]);
    in[v] = true;
  }
  CHECK_EQ(cut_weight(w, in), cut);
  return cut;
}

Mat random_graph(int n, int density, ll maxw) {
  Mat w(n, vector<ll>(n));
  rep(u, 0, n) {
    w[u][u] = rnd(0, 1) ? rnd(0, maxw) : LLONG_MAX; // ignored
    rep(v, 0, u)
      if (rnd(1, 100) <= density) w[u][v] = w[v][u] = rnd(0, maxw);
  }
  return w;
}

int main() {
  CHECK_EQ(solve({{0, 0}, {0, 0}}), 0LL);
  CHECK_EQ(solve({{5, 7}, {7, 9}}), 7LL);
  // A sentinel on the diagonal must not overflow when merged.
  CHECK_EQ(solve({{LLONG_MAX, 1}, {1, LLONG_MAX}}), 1LL);
  rep(it, 0, iters(3000)) {
    int n = (int)rnd(2, 10), density = (int)rnd(10, 100);
    ll maxw = rnd(0, 3) ? rnd(1, 10) : (ll)1e15;
    Mat w = random_graph(n, density, maxw);
    CHECK_EQ(solve(w), brute(w));
  }
  rep(it, 0, iters(60)) {
    int n = (int)rnd(10, 40);
    Mat w = random_graph(n, (int)rnd(5, 100), rnd(1, 1000));
    CHECK_EQ(solve(w), by_flows(w));
  }
  // Extremes: total weight close to 2^63.
  {
    ll big = 1LL << 59;
    Mat w(4, vector<ll>(4, big));
    rep(v, 0, 4) w[v][v] = 0;
    CHECK_EQ(solve(w), 3 * big);
  }
  // Large: two dense halves joined by a few light edges.
  elapsed();
  {
    int n = 500;
    Mat w = random_graph(n, 50, 1000);
    vi half(n);
    rep(v, 0, n) half[v] = (int)rnd(0, 1);
    rep(u, 0, n)
      rep(v, 0, n)
        if (half[u] != half[v]) w[u][v] = 0;
    ll planted = 0;
    rep(it, 0, 10) {
      int u = (int)rnd(0, n - 1), v = (int)rnd(0, n - 1);
      if (half[u] != half[v] && !w[u][v]) w[u][v] = w[v][u] = 1, planted++;
    }
    CHECK_EQ(solve(w), planted);
  }
  fprintf(stderr, "  large case: %.2fs\n", elapsed());
  puts("OK flow/global-min-cut");
}
