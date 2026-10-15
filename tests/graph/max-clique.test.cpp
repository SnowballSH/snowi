#include "../testlib.hpp"
#include "../../content/graph/max-clique.hpp"

typedef MaxClique::B B;

// Largest clique by DP over all subsets: a set is a clique iff
// its lowest vertex sees the rest and the rest is a clique.
int brute(const vector<B>& g) {
  int n = sz(g), best = 0;
  vi nb(n);
  rep(u, 0, n)
    rep(v, 0, n)
      if (g[u][v]) nb[u] |= 1 << v;
  vector<char> ok(1 << n);
  ok[0] = 1;
  rep(mask, 1, 1 << n) {
    int low = __builtin_ctz(mask), rest = mask & (mask - 1);
    ok[mask] = ok[rest] && (nb[low] & rest) == rest;
    if (ok[mask]) best = max(best, __builtin_popcount(mask));
  }
  return best;
}

// Bron-Kerbosch with a pivot: an independent check for larger n.
void bk(const vector<B>& g, B p, B x, int size, int& best) {
  if (p.none() && x.none()) best = max(best, size);
  if (p.none()) return;
  int piv = -1;
  rep(u, 0, sz(g))
    if (p[u] || x[u]) {
      if (piv < 0 || (p & g[u]).count() > (p & g[piv]).count()) piv = u;
    }
  rep(v, 0, sz(g))
    if (p[v] && !g[piv][v]) {
      bk(g, p & g[v], x & g[v], size + 1, best);
      p[v] = 0, x[v] = 1;
    }
}

vector<B> random_graph(int n, double prob) {
  vector<B> g(n);
  rep(u, 0, n)
    rep(v, 0, u)
      if (rnd_real(0, 1) < prob) g[u][v] = g[v][u] = 1;
  return g;
}

// The answer must be a sorted clique of the expected size.
void check(const vector<B>& g, const vi& c, int want) {
  CHECK_EQ(sz(c), want);
  CHECK(is_sorted(all(c)));
  rep(i, 0, sz(c)) {
    CHECK(0 <= c[i] && c[i] < sz(g));
    rep(j, 0, i) CHECK(g[c[i]][c[j]]);
  }
}

int main() {
  CHECK(MaxClique({}).solve().empty());
  CHECK(MaxClique(vector<B>(1)).solve() == vi{0});
  CHECK_EQ(sz(MaxClique(vector<B>(5)).solve()), 1);
  vector<B> full(7);
  rep(u, 0, 7)
    rep(v, 0, 7) full[u][v] = u != v;
  CHECK(MaxClique(full).solve() == (vi{0, 1, 2, 3, 4, 5, 6}));
  rep(it, 0, iters(3000)) {
    int n = (int)rnd(1, 12);
    auto g = random_graph(n, rnd_real(0, 1));
    check(g, MaxClique(g).solve(), brute(g));
  }
  rep(it, 0, iters(40)) {
    int n = (int)rnd(20, 60);
    auto g = random_graph(n, rnd_real(0.3, 0.9));
    B all_;
    rep(i, 0, n) all_[i] = 1;
    int want = 0;
    bk(g, all_, B(), 0, want);
    check(g, MaxClique(g).solve(), want);
  }
  // A planted clique of 24 in G(128, 1/2), whose own cliques
  // have about 11 vertices.
  int n = 128;
  auto g = random_graph(n, 0.5);
  vi perm(n);
  iota(all(perm), 0);
  shuffle(all(perm), test_rng);
  rep(i, 0, 24)
    rep(j, 0, i) g[perm[i]][perm[j]] = 1, g[perm[j]][perm[i]] = 1;
  vi want(perm.begin(), perm.begin() + 24);
  sort(all(want));
  CHECK(MaxClique(g).solve() == want);
  // Dense random graph: the hard case.
  g = random_graph(n, 0.9);
  double t0 = elapsed();
  vi c = MaxClique(g).solve();
  fprintf(stderr, "  n = 128, p = 0.9: %.2fs, size %d\n", elapsed() - t0,
          sz(c));
  check(g, c, sz(c));
  puts("OK graph/max-clique");
}
