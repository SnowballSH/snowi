#include "../testlib.hpp"
#include "../../content/tree/lca-rmq.hpp"

// Shapes: 0 random, 1 path, 2 star, 3 long thin; labels shuffled.
vector<vi> random_tree(int n, int shape) {
  vi perm(n);
  iota(all(perm), 0), shuffle(all(perm), test_rng);
  vector<vi> adj(n);
  rep(i, 1, n) {
    int p = shape == 1   ? i - 1
            : shape == 2 ? 0
            : shape == 3 ? (int)rnd(max(0, i - 3), i - 1)
                         : (int)rnd(0, i - 1);
    adj[perm[i]].push_back(perm[p]), adj[perm[p]].push_back(perm[i]);
  }
  for (auto& a : adj) shuffle(all(a), test_rng);
  return adj;
}

struct Brute {
  vi par, dep;
  Brute(const vector<vi>& adj, int root) : par(sz(adj), -1), dep(sz(adj)) {
    vi q{root}, seen(sz(adj));
    seen[root] = 1;
    rep(i, 0, sz(q))
      for (int v : adj[q[i]])
        if (!seen[v]) {
          seen[v] = 1, par[v] = q[i], dep[v] = dep[q[i]] + 1;
          q.push_back(v);
        }
  }
  int lca(int u, int v) {
    while (dep[u] > dep[v]) u = par[u];
    while (dep[v] > dep[u]) v = par[v];
    while (u != v) u = par[u], v = par[v];
    return u;
  }
  bool is_anc(int a, int v) {
    while (v >= 0 && v != a) v = par[v];
    return v == a;
  }
};

void check_tree(const vector<vi>& adj, int root) {
  int n = sz(adj);
  RmqLCA lc(adj, root);
  Brute b(adj, root);
  CHECK_EQ(sz(lc.order), n);
  rep(v, 0, n) {
    CHECK_EQ(lc.depth[v], b.dep[v]);
    CHECK_EQ(lc.par[v], b.par[v]);
    CHECK_EQ(lc.order[lc.tin[v]], v);
    int size = 0;
    rep(x, 0, n) size += b.is_anc(v, x);
    rep(x, 0, n) {
      bool in = lc.tin[v] <= lc.tin[x] && lc.tin[x] < lc.tin[v] + size;
      CHECK_EQ(in, b.is_anc(v, x));
    }
  }
  rep(u, 0, n)
    rep(v, 0, n) {
      int a = b.lca(u, v);
      CHECK_EQ(lc.lca(u, v), a);
      CHECK_EQ(lc.dist(u, v), b.dep[u] + b.dep[v] - 2 * b.dep[a]);
    }
  // The header's Example: Steiner tree size of a random set s.
  vi s;
  rep(v, 0, n)
    if (rnd(0, 2) == 0) s.push_back(v);
  shuffle(all(s), test_rng);
  sort(all(s), [&](int x, int y) { return lc.tin[x] < lc.tin[y]; });
  ll edges = 0;
  rep(i, 0, sz(s)) edges += lc.dist(s[i], s[(i + 1) % sz(s)]);
  edges /= 2;
  // Edge v - par[v] is needed iff subtree(v) holds some, not all, of s.
  ll want = 0;
  rep(v, 0, n)
    if (v != root) {
      int in = 0;
      for (int x : s) in += b.is_anc(v, x);
      want += 0 < in && in < sz(s);
    }
  CHECK_EQ(edges, want);
}

int main() {
  rep(it, 0, iters(1500)) {
    int n = (int)rnd(1, it < 100 ? 4 : 24);
    check_tree(random_tree(n, it % 4), (int)rnd(0, n - 1));
  }
  // Large: random recursive tree (shallow: naive walk is cheap), and
  // a 2e5 path rooted at one end (deep: the build must be iterative).
  int n = 200'000;
  auto adj = random_tree(n, 0);
  int root = (int)rnd(0, n - 1);
  RmqLCA lc(adj, root);
  Brute b(adj, root);
  rep(it, 0, 100'000) {
    int u = (int)rnd(0, n - 1), v = (int)rnd(0, n - 1);
    int a = b.lca(u, v);
    CHECK_EQ(lc.lca(u, v), a);
    CHECK_EQ(lc.dist(u, v), b.dep[u] + b.dep[v] - 2 * b.dep[a]);
  }
  vi perm(n);
  iota(all(perm), 0), shuffle(all(perm), test_rng);
  vector<vi> path(n);
  rep(i, 1, n)
    path[perm[i]].push_back(perm[i - 1]), path[perm[i - 1]].push_back(perm[i]);
  RmqLCA lp(path, perm[0]);
  rep(it, 0, 100'000) {
    int i = (int)rnd(0, n - 1), j = (int)rnd(0, n - 1);
    CHECK_EQ(lp.lca(perm[i], perm[j]), perm[min(i, j)]);
    CHECK_EQ(lp.dist(perm[i], perm[j]), abs(i - j));
  }
  puts("OK tree/lca-rmq");
}
