#include "../testlib.hpp"
#include "../../content/tree/lca.hpp"

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
  vi path(int u, int v) { // vertices u .. v
    int a = lca(u, v);
    vi left, right;
    for (; u != a; u = par[u]) left.push_back(u);
    for (; v != a; v = par[v]) right.push_back(v);
    left.push_back(a);
    left.insert(left.end(), right.rbegin(), right.rend());
    return left;
  }
};

void check_tree(const vector<vi>& adj, int root) {
  int n = sz(adj);
  LCA lc(adj, root);
  Brute b(adj, root);
  CHECK_EQ(sz(lc.order), n);
  rep(v, 0, n) {
    CHECK_EQ(lc.depth[v], b.dep[v]);
    CHECK_EQ(lc.up[0][v], v == root ? root : b.par[v]);
    CHECK_EQ(lc.order[lc.tin[v]], v);
  }
  // Preorder: each subtree is a contiguous block of tin values.
  rep(v, 0, n) {
    int size = 0;
    rep(x, 0, n) size += b.is_anc(v, x);
    rep(x, 0, n) {
      bool in = lc.tin[v] <= lc.tin[x] && lc.tin[x] < lc.tin[v] + size;
      CHECK_EQ(in, b.is_anc(v, x));
    }
  }
  rep(u, 0, n) {
    rep(k, 0, b.dep[u] + 3) {
      int want = u;
      rep(s, 0, k) want = want < 0 ? -1 : b.par[want];
      if (k > b.dep[u]) want = -1;
      CHECK_EQ(lc.kth_ancestor(u, k), want);
    }
    rep(v, 0, n) {
      vi p = b.path(u, v);
      CHECK_EQ(lc.lca(u, v), b.lca(u, v));
      CHECK_EQ(lc.dist(u, v), sz(p) - 1);
      rep(k, 0, sz(p) + 2) CHECK_EQ(lc.jump(u, v, k), k < sz(p) ? p[k] : -1);
    }
  }
}

// The header's Example: max edge weight on a path.
void check_example(const vector<vi>& adj, int root) {
  int n = sz(adj);
  LCA lc(adj, root);
  Brute b(adj, root);
  vi w(n);
  rep(v, 0, n) w[v] = (int)rnd(-1e9, 1e9);
  vector<vi> mx(sz(lc.up), w); // max of 2^k edges above v
  rep(k, 1, sz(mx))
    rep(v, 0, n) mx[k][v] = max(mx[k - 1][v], mx[k - 1][lc.up[k - 1][v]]);
  auto climb = [&](int v, int d) { // max of d edges above v
    int r = INT_MIN;
    rep(k, 0, sz(mx))
      if (d >> k & 1) r = max(r, mx[k][v]), v = lc.up[k][v];
    return r;
  };
  rep(it, 0, 3 * n) {
    int u = (int)rnd(0, n - 1), v = (int)rnd(0, n - 1);
    int a = lc.lca(u, v), d = lc.depth[a];
    int best = max(climb(u, lc.depth[u] - d), climb(v, lc.depth[v] - d));
    int want = INT_MIN;
    for (int x = u; x != a; x = b.par[x]) want = max(want, w[x]);
    for (int x = v; x != a; x = b.par[x]) want = max(want, w[x]);
    CHECK_EQ(best, want);
  }
}

int main() {
  rep(it, 0, iters(1500)) {
    int n = (int)rnd(1, it < 100 ? 4 : 24);
    auto adj = random_tree(n, it % 4);
    int root = (int)rnd(0, n - 1);
    check_tree(adj, root);
    check_example(adj, root);
  }
  // Large: a random recursive tree (shallow, so the naive walk is
  // cheap), then a 2e5 path rooted at one end (deep: no recursion).
  int n = 200'000;
  auto adj = random_tree(n, 0);
  int root = (int)rnd(0, n - 1);
  LCA lc(adj, root);
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
  LCA lp(path, perm[0]);
  rep(it, 0, 100'000) {
    int i = (int)rnd(0, n - 1), j = (int)rnd(0, n - 1);
    CHECK_EQ(lp.lca(perm[i], perm[j]), perm[min(i, j)]);
    CHECK_EQ(lp.dist(perm[i], perm[j]), abs(i - j));
    int k = (int)rnd(0, n);
    CHECK_EQ(lp.kth_ancestor(perm[i], k), k <= i ? perm[i - k] : -1);
    int s = i <= j ? i + k : i - k;
    CHECK_EQ(lp.jump(perm[i], perm[j], k), k <= abs(i - j) ? perm[s] : -1);
  }
  puts("OK tree/lca");
}
