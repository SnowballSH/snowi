#include "../testlib.hpp"
#include "../../content/graph/bellman-ford.hpp"

typedef vector<tuple<int, int, ll>> Edges;
const ll NONE = LLONG_MAX;

// walk[k][x][v]: lightest walk with exactly k edges from x to v.
vector<ll> brute(int n, const Edges& edges, int s) {
  vector<vector<vector<ll>>> walk(n + 1,
                                  vector<vector<ll>>(n, vector<ll>(n, NONE)));
  rep(x, 0, n) walk[0][x][x] = 0;
  rep(k, 0, n)
    rep(x, 0, n)
      for (auto [a, b, w] : edges)
        if (walk[k][x][a] != NONE)
          walk[k + 1][x][b] = min(walk[k + 1][x][b], walk[k][x][a] + w);
  auto reach = [&](int x, int v) {
    rep(k, 0, n)
      if (walk[k][x][v] != NONE) return true;
    return false;
  };
  vector<bool> on_neg(n);
  rep(x, 0, n)
    rep(k, 1, n + 1)
      on_neg[x] = on_neg[x] || (walk[k][x][x] != NONE && walk[k][x][x] < 0);
  vector<ll> d(n, NONE);
  rep(v, 0, n) {
    rep(k, 0, n) d[v] = min(d[v], walk[k][s][v]);
    rep(y, 0, n)
      if (on_neg[y] && reach(s, y) && reach(y, v)) d[v] = LLONG_MIN;
  }
  return d;
}

// Finite vertices: par edges are tight and lead back to s.
void check_par(int n, const Edges& edges, int s, BellmanFord& bf) {
  auto finite = [&](int v) {
    return bf.dist[v] != LLONG_MAX && bf.dist[v] != LLONG_MIN;
  };
  vector<bool> tight(n);
  for (auto [a, b, w] : edges)
    if (finite(b) && a == bf.par[b] && bf.dist[a] + w == bf.dist[b])
      tight[b] = true;
  rep(v, 0, n)
    if (finite(v)) {
      int u = v;
      for (int step = 0; u != s; step++) {
        CHECK(step < n && tight[u]);
        u = bf.par[u];
      }
    }
  if (finite(s)) CHECK_EQ(bf.par[s], -1);
}

// Real weights must not be truncated into ll silently (see Tweak).
static_assert(!is_constructible_v<BellmanFord, int,
                                  vector<tuple<int, int, double>>, int>);
static_assert(
    is_constructible_v<BellmanFord, int, vector<tuple<int, int, ll>>, int>);

int main() {
  rep(it, 0, iters(3000)) {
    int n = (int)rnd(1, 7), m = (int)rnd(0, 14);
    ll lo = rnd(0, 2) ? -5 : 0, hi = 10;
    if (rnd(0, 4) == 0) lo = -1'000'000'000'000, hi = -lo;
    Edges edges;
    rep(i, 0, m)
      edges.push_back({(int)rnd(0, n - 1), (int)rnd(0, n - 1), rnd(lo, hi)});
    int s = (int)rnd(0, n - 1);
    BellmanFord bf(n, edges, s);
    CHECK(bf.dist == brute(n, edges, s));
    check_par(n, edges, s, bf);
  }
  // A cycle listed backwards with total weight -1: reaching every
  // vertex takes n - 1 rounds, marking all of them n - 1 more.
  int n = 300;
  Edges edges;
  for (int v = n - 1; v > 0; v--) edges.push_back({v - 1, v, -1});
  edges.push_back({n - 1, 0, n - 2});
  CHECK(BellmanFord(n, edges, 0).dist == vector<ll>(n, LLONG_MIN));
  // Speed: a chain listed backwards needs n - 1 full rounds.
  n = 1000, edges.clear();
  for (int v = n - 1; v > 0; v--) edges.push_back({v - 1, v, -1});
  while (sz(edges) < 5000) {
    int a = (int)rnd(0, n - 1), b = (int)rnd(0, a);
    edges.push_back({a, b, 1'000'000'000});
  }
  BellmanFord bf(n, edges, 0);
  rep(v, 0, n) CHECK_EQ(bf.dist[v], (ll)-v);
  check_par(n, edges, 0, bf);
  puts("OK graph/bellman-ford");
}
