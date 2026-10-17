#include "../testlib.hpp"
#include "../../content/geometry/manhattan-mst.hpp"

typedef P<ll> Pl;

ll dist(Pl a, Pl b) { return abs(a.x - b.x) + abs(a.y - b.y); }
// Prim on the complete graph, O(n^2).
ll brute(const vector<Pl>& ps) {
  int n = sz(ps);
  vector<ll> d(n, LLONG_MAX);
  vector<bool> in(n);
  ll total = 0;
  d[0] = 0;
  rep(it, 0, n) {
    int u = -1;
    rep(v, 0, n)
      if (!in[v] && (u < 0 || d[v] < d[u])) u = v;
    in[u] = true, total += d[u];
    rep(v, 0, n) d[v] = min(d[v], dist(ps[u], ps[v]));
  }
  return total;
}
ll kruskal(int n, vector<array<ll, 3>> e) {
  vi par(n);
  iota(all(par), 0);
  auto find = [&](int x) {
    while (par[x] != x) x = par[x] = par[par[x]];
    return x;
  };
  sort(all(e));
  ll total = 0;
  int joined = 0;
  for (auto [w, i, j] : e) {
    int a = find((int)i), b = find((int)j);
    if (a != b) par[a] = b, total += w, joined++;
  }
  CHECK(joined == n - 1);
  return total;
}
void check(const vector<Pl>& ps) {
  auto e = manhattan_mst(ps);
  CHECK(sz(e) <= 4 * sz(ps));
  for (auto [w, i, j] : e) CHECK(i != j && w == dist(ps[i], ps[j]));
  CHECK_EQ(kruskal(sz(ps), e), brute(ps));
}

int main() {
  check({{5, 5}});
  check({{0, 0}, {0, 0}});
  check({{0, 0}, {1, 1}, {2, 2}, {3, 3}});
  check({{0, 0}, {1, -1}, {2, -2}, {3, -3}});
  rep(it, 0, iters(4000)) {
    int n = (int)rnd(1, 30);
    ll c = vector<ll>{1, 3, 10, 1000, (ll)1e9}[it % 5];
    vector<Pl> ps(n);
    for (Pl& p : ps) p = {rnd(-c, c), rnd(-c, c)};
    if (it % 7 == 0) // points on diagonals and axes: octant borders
      for (Pl& p : ps) p.y = rnd(0, 1) ? p.x : (rnd(0, 1) ? -p.x : 0);
    check(ps);
  }
  int n = 200'000;
  vector<Pl> ps(n);
  for (Pl& p : ps) p = {rnd(-1e9, 1e9), rnd(-1e9, 1e9)};
  double start = elapsed();
  auto e = manhattan_mst(ps);
  ll total = kruskal(n, e);
  fprintf(stderr, "  2e5 points: %.2fs (%lld)\n", elapsed() - start, total);
  puts("OK geometry/manhattan-mst");
}
