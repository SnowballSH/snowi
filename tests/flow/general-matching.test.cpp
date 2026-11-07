#include "../testlib.hpp"
#include "../../content/flow/general-matching.hpp"

typedef vector<pii> Edges;

vector<vi> adjacency(int n, const Edges& edges) {
  vector<vi> g(n);
  for (auto [u, v] : edges) {
    g[u].push_back(v);
    if (u != v) g[v].push_back(u);
  }
  rep(u, 0, n) shuffle(all(g[u]), test_rng);
  return g;
}

// Maximum matching by DP over the set of vertices still free.
int brute(int n, const Edges& edges) {
  vi nb(n);
  for (auto [u, v] : edges)
    if (u != v) nb[u] |= 1 << v, nb[v] |= 1 << u;
  vi f(1 << n);
  rep(mask, 1, 1 << n) {
    int i = __builtin_ctz(mask), rest = mask & ~(1 << i);
    f[mask] = f[rest];
    for (int m = nb[i] & rest; m; m &= m - 1)
      f[mask] = max(f[mask], 1 + f[rest & ~(m & -m)]);
  }
  return f[(1 << n) - 1];
}

// Rank of a random Tutte matrix mod p: twice the maximum matching
// size, except with probability below n / p.
int tutte(int n, const Edges& edges) {
  const ll p = 1'000'000'007;
  vector<vector<ll>> a(n, vector<ll>(n));
  for (auto [u, v] : edges)
    if (u != v) {
      ll x = rnd(1, p - 1);
      a[u][v] = (a[u][v] + x) % p, a[v][u] = (a[v][u] + p - x) % p;
    }
  auto power = [&](ll b, ll e) {
    ll r = 1;
    for (; e; e >>= 1, b = b * b % p)
      if (e & 1) r = r * b % p;
    return r;
  };
  int rank = 0;
  rep(c, 0, n) {
    int r = rank;
    while (r < n && !a[r][c]) r++;
    if (r == n) continue;
    swap(a[r], a[rank]);
    ll inv = power(a[rank][c], p - 2);
    rep(i, rank + 1, n) {
      ll f = a[i][c] * inv % p;
      if (f)
        rep(j, c, n) a[i][j] = (a[i][j] + (p - f) * a[rank][j]) % p;
    }
    rank++;
  }
  return rank;
}

// mate is a matching of k pairs that uses edges of the graph.
void certify(int n, const Edges& edges, int k, const vi& mate) {
  CHECK_EQ(sz(mate), n);
  set<pii> es;
  for (auto [u, v] : edges) es.insert({u, v}), es.insert({v, u});
  int ends = 0;
  rep(v, 0, n)
    if (mate[v] >= 0) {
      CHECK(mate[v] < n && mate[v] != v);
      CHECK_EQ(mate[mate[v]], v);
      CHECK(es.count({v, mate[v]}));
      ends++;
    } else CHECK_EQ(mate[v], -1);
  CHECK_EQ(ends, 2 * k);
}

int solve(int n, const Edges& edges) {
  auto [k, mate] = general_matching(adjacency(n, edges));
  certify(n, edges, k, mate);
  return k;
}

Edges random_edges(int n, int m) {
  Edges e;
  if (n)
    rep(i, 0, m) e.push_back({(int)rnd(0, n - 1), (int)rnd(0, n - 1)});
  return e;
}

// Odd cycles a, v, ..., end, a through an earlier vertex a:
// nested blossoms.
Edges flowers(int n) {
  Edges e;
  for (int v = 1; v < n;) {
    int a = (int)rnd(0, v - 1), len = 2 * (int)rnd(1, 3);
    int end = min(n, v + len) - 1;
    e.push_back({a, v});
    rep(x, v, end) e.push_back({x, x + 1});
    e.push_back({end, a});
    v = end + 1;
  }
  return e;
}

int main() {
  CHECK_EQ(solve(0, {}), 0);
  CHECK_EQ(solve(1, {{0, 0}}), 0);
  CHECK_EQ(solve(3, {{0, 1}, {1, 2}, {2, 0}}), 1);
  // Two triangles joined by an edge: the blossom must be entered.
  CHECK_EQ(solve(6, {{0, 1}, {1, 2}, {2, 0}, {2, 3}, {3, 4}, {4, 5}, {5, 3}}),
           3);
  { // Both sides of a blossom must be shrunk: shrinking only the
    // scanned vertex's side finds 2 pairs with these exact lists.
    vector<vi> g = {{5, 3}, {2, 5}, {1, 4, 3}, {2, 5, 0}, {2}, {1, 0, 3}};
    CHECK_EQ(general_matching(g).first, 3);
  }
  rep(it, 0, iters(4000)) {
    int n = (int)rnd(0, 13);
    int m = (int)rnd(0, n * (n + 1) / 2 + 2);
    Edges e = rnd(0, 4) ? random_edges(n, m) : flowers(n);
    CHECK_EQ(solve(n, e), brute(n, e));
  }
  rep(it, 0, iters(300)) {
    int n = (int)rnd(10, 80);
    Edges e = rnd(0, 1) ? random_edges(n, (int)rnd(n / 2, 3 * n)) : flowers(n);
    CHECK_EQ(2 * solve(n, e), tutte(n, e));
  }
  rep(n, 0, 30) { // complete graphs
    Edges e;
    rep(i, 0, n)
      rep(j, 0, i) e.push_back({i, j});
    CHECK_EQ(solve(n, e), n / 2);
  }
  // Large: a planted perfect matching hidden among random edges
  // and odd cycles; then disjoint triangles (no perfect matching).
  elapsed();
  {
    int n = 100'000;
    vi p(n);
    iota(all(p), 0), shuffle(all(p), test_rng);
    Edges e = random_edges(n, 100'000), f = flowers(n);
    e.insert(e.end(), all(f));
    for (int i = 0; i < n; i += 2) e.push_back({p[i], p[i + 1]});
    CHECK_EQ(solve(n, e), n / 2);
  }
  {
    int n = 60'000;
    Edges e;
    for (int i = 0; i < n; i += 3)
      e.push_back({i, i + 1}), e.push_back({i + 1, i + 2}),
          e.push_back({i + 2, i});
    CHECK_EQ(solve(n, e), n / 3);
  }
  {
    int n = 1000;
    CHECK_EQ(solve(n, random_edges(n, 100'000)), n / 2);
  }
  { // Broom: a path 0..L-1, then n - L leaves on L - 1. Each leaf's
    // search fails after walking the path, unless the trees of
    // failed searches are retired.
    int n = 100'000, L = n / 2;
    Edges e;
    rep(i, 0, L - 1) e.push_back({i, i + 1});
    rep(i, L, n) e.push_back({L - 1, i});
    double start = elapsed();
    CHECK_EQ(solve(n, e), (L + 1) / 2);
    CHECK(elapsed() - start < 2); // without retiring: > 2 min
  }
  fprintf(stderr, "  large cases: %.2fs\n", elapsed());
  puts("OK flow/general-matching");
}
