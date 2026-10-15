#include "../testlib.hpp"
#include "../../content/graph/dominator-tree.hpp"

// Vertices reachable from s without passing through `banned`.
vector<bool> reach(const vector<vi>& adj, int s, int banned) {
  vector<bool> seen(sz(adj));
  if (s == banned) return seen;
  vi st = {s};
  seen[s] = true;
  while (sz(st)) {
    int u = st.back();
    st.pop_back();
    for (int v : adj[u])
      if (v != banned && !seen[v]) seen[v] = true, st.push_back(v);
  }
  return seen;
}

// d strictly dominates v iff removing d cuts v off from s. The
// strict dominators form a chain; idom is the one with the most
// dominators itself.
vi brute(const vector<vi>& adj, int s) {
  int n = sz(adj);
  auto all_ = reach(adj, s, -1);
  vector<vi> dom(n); // strict dominators
  rep(d, 0, n) {
    if (!all_[d]) continue;
    auto r = reach(adj, s, d);
    rep(v, 0, n)
      if (all_[v] && v != d && !r[v]) dom[v].push_back(d);
  }
  vi idom(n, -1);
  rep(v, 0, n)
    if (all_[v]) {
      if (v == s) {
        idom[v] = s;
        continue;
      }
      int best = -1;
      for (int d : dom[v])
        if (best < 0 || sz(dom[d]) > sz(dom[best])) best = d;
      idom[v] = best;
    }
  return idom;
}

vector<vi> random_graph(int n, int m) {
  vector<vi> adj(n);
  rep(i, 0, m) adj[rnd(0, n - 1)].push_back((int)rnd(0, n - 1));
  return adj;
}

int main() {
  CHECK(dominator_tree({{}}, 0) == vi{0});
  CHECK(dominator_tree({{0, 0}}, 0) == vi{0});
  CHECK(dominator_tree({{1}, {}, {1}}, 0) == (vi{0, 0, -1}));
  CHECK(dominator_tree({{1, 2}, {3}, {3}, {}}, 0) == (vi{0, 0, 0, 0}));
  rep(it, 0, iters(5000)) {
    int n = (int)rnd(1, 10), m = (int)rnd(0, 3 * n);
    auto adj = random_graph(n, m);
    int s = (int)rnd(0, n - 1);
    CHECK(dominator_tree(adj, s) == brute(adj, s));
  }
  rep(it, 0, iters(100)) {
    int n = (int)rnd(10, 300);
    auto adj = random_graph(n, (int)rnd(n, 3 * n));
    CHECK(dominator_tree(adj, 0) == brute(adj, 0));
  }
  // Random graph, n = 2000: about 1.6e7 steps of brute force.
  auto adj = random_graph(2000, 6000);
  CHECK(dominator_tree(adj, 0) == brute(adj, 0));
  // A path plus back edges to ancestors: idom is the tree parent.
  // Deep DFS and long compression paths.
  int n = 200'000;
  vi parent(n, -1);
  vector<vi> g(n);
  rep(v, 1, n) {
    parent[v] = v < n / 2 ? v - 1 : (int)rnd(0, v - 1);
    g[parent[v]].push_back(v);
  }
  rep(v, 1, n) {
    int a = v;
    rep(j, 0, 3) a = a ? parent[a] : 0;
    g[v].push_back(a);
  }
  g[n / 2 - 1].push_back(1);
  for (vi& out : g) shuffle(all(out), test_rng);
  double t0 = elapsed();
  vi idom = dominator_tree(g, 0);
  fprintf(stderr, "  n = 2e5, m = 4e5: %.2fs\n", elapsed() - t0);
  parent[0] = 0;
  CHECK(idom == parent);
  // Star: n DFS children of the root share one bucket.
  vector<vi> star(n);
  rep(v, 1, n) star[0].push_back(v);
  CHECK(dominator_tree(star, 0) == vi(n, 0));
  // Path with edges from its end to every vertex: each step
  // evaluates the whole path unless it is compressed.
  vector<vi> path(n);
  rep(v, 1, n) path[v - 1].push_back(v), path[n - 1].push_back(v);
  vi want(n);
  rep(v, 1, n) want[v] = v - 1;
  CHECK(dominator_tree(path, 0) == want);
  puts("OK graph/dominator-tree");
}
