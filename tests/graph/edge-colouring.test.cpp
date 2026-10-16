#include "../testlib.hpp"
#include "../../content/graph/edge-colouring.hpp"

// Proper (distinct colours at every vertex) with exactly the colours
// [0, Delta), Delta the maximum degree: Konig's bound is tight.
void check(int L, int R, const vector<pii>& edges) {
  vi col = edge_colouring(L, R, edges);
  CHECK_EQ(sz(col), sz(edges));
  vi deg(L + R);
  for (auto [u, v] : edges) deg[u]++, deg[L + v]++;
  int D = sz(deg) ? *max_element(all(deg)) : 0;
  set<pii> seen; // (vertex, colour)
  rep(e, 0, sz(edges)) {
    CHECK(0 <= col[e] && col[e] < D);
    CHECK(seen.insert({edges[e].first, col[e]}).second);
    CHECK(seen.insert({L + edges[e].second, col[e]}).second);
  }
}

// vizing: proper, with colours in [0, Delta].
void check_vizing(int n, const vector<pii>& edges) {
  vi col = vizing(n, edges);
  CHECK_EQ(sz(col), sz(edges));
  vi deg(n);
  for (auto [u, v] : edges) deg[u]++, deg[v]++;
  int D = *max_element(all(deg));
  set<pii> seen;
  rep(e, 0, sz(edges)) {
    CHECK(0 <= col[e] && col[e] <= D);
    CHECK(seen.insert({edges[e].first, col[e]}).second);
    CHECK(seen.insert({edges[e].second, col[e]}).second);
  }
}

int main() {
  check(0, 0, {});
  check(3, 2, {});
  check(1, 1, {{0, 0}, {0, 0}, {0, 0}});
  rep(it, 0, iters(3000)) {
    int L = (int)rnd(1, 8), R = (int)rnd(1, 8), m = (int)rnd(0, 30);
    vector<pii> edges(m);
    for (auto& [u, v] : edges) u = (int)rnd(0, L - 1), v = (int)rnd(0, R - 1);
    check(L, R, edges);
  }
  // k-regular: k perfect matchings (union of k random permutations).
  rep(it, 0, iters(300)) {
    int n = (int)rnd(1, 30), k = (int)rnd(1, 8);
    vector<pii> edges;
    vi p(n);
    iota(all(p), 0);
    rep(j, 0, k) {
      shuffle(all(p), test_rng);
      rep(i, 0, n) edges.push_back({i, p[i]});
    }
    check(n, n, edges);
  }
  check_vizing(1, {});
  check_vizing(3, {{0, 1}, {1, 2}, {2, 0}});
  rep(n, 1, 14) { // complete graphs: odd n needs Delta + 1
    vector<pii> edges;
    rep(a, 0, n)
      rep(b, 0, a) edges.push_back({a, b});
    check_vizing(n, edges);
  }
  rep(it, 0, iters(3000)) {
    int n = (int)rnd(1, 12);
    double p = rnd_real(0, 1);
    vector<pii> edges;
    rep(a, 0, n)
      rep(b, 0, a)
        if (rnd_real(0, 1) < p) {
          if (rnd(0, 1)) edges.push_back({a, b});
          else edges.push_back({b, a});
        }
    shuffle(all(edges), test_rng);
    check_vizing(n, edges);
  }
  int L = 2000, R = 2000, m = 200'000;
  vector<pii> edges(m);
  for (auto& [u, v] : edges) u = (int)rnd(0, L - 1), v = (int)rnd(0, R - 1);
  double start = elapsed();
  vi col = edge_colouring(L, R, edges);
  fprintf(stderr, "  2e5 edges on 2000 + 2000 vertices: %.2fs\n",
          elapsed() - start);
  check(L, R, edges);
  set<pii> simple;
  for (auto [u, v] : edges)
    if (u != v) simple.insert(minmax(u, v));
  vector<pii> es(all(simple));
  start = elapsed();
  col = vizing(L, es);
  fprintf(stderr, "  vizing, %d edges on 2000 vertices: %.2fs\n", sz(es),
          elapsed() - start);
  check_vizing(L, es);
  puts("OK graph/edge-colouring");
}
