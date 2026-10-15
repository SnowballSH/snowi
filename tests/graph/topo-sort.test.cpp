#include "../testlib.hpp"
#include "../../content/graph/topo-sort.hpp"

// A cycle exists iff some vertex reaches itself in >= 1 step.
bool brute_has_cycle(const vector<vi>& adj) {
  int n = sz(adj);
  vector<vector<bool>> r(n, vector<bool>(n));
  rep(u, 0, n)
    for (int v : adj[u]) r[u][v] = true;
  rep(k, 0, n)
    rep(i, 0, n)
      rep(j, 0, n)
        if (r[i][k] && r[k][j]) r[i][j] = true;
  rep(v, 0, n)
    if (r[v][v]) return true;
  return false;
}

void check_order(const vector<vi>& adj, const vi& order) {
  int n = sz(adj);
  CHECK_EQ(sz(order), n);
  vi pos(n, -1);
  rep(i, 0, n) {
    CHECK(pos[order[i]] == -1);
    pos[order[i]] = i;
  }
  rep(u, 0, n)
    for (int v : adj[u]) CHECK(pos[u] < pos[v]);
}

// Random graph; if dag, edges follow a random permutation.
vector<vi> random_graph(int n, int m, bool dag) {
  vi rank(n);
  iota(all(rank), 0);
  shuffle(all(rank), test_rng);
  vector<vi> adj(n);
  rep(i, 0, m) {
    int u = (int)rnd(0, n - 1), v = (int)rnd(0, n - 1);
    if (dag && rank[u] >= rank[v]) {
      if (u == v) continue;
      swap(u, v);
    }
    adj[u].push_back(v);
  }
  return adj;
}

int main() {
  int cyclic = 0;
  rep(it, 0, iters(4000)) {
    int n = (int)rnd(0, 8), m = n ? (int)rnd(0, 16) : 0;
    vector<vi> adj = random_graph(n, m, rnd(0, 1));
    vi order = topo_sort(adj);
    if (brute_has_cycle(adj)) {
      CHECK(order.empty());
      cyclic++;
    } else {
      check_order(adj, order);
    }
  }
  int total = iters(4000);
  CHECK(cyclic > total / 8 && cyclic < total * 7 / 8);
  // Large DAG, and a long chain closed into a cycle.
  int n = 200'000;
  vector<vi> adj = random_graph(n, 400'000, true);
  check_order(adj, topo_sort(adj));
  vector<vi> chain(n);
  rep(v, 0, n - 1) chain[v].push_back(v + 1);
  CHECK(topo_sort(chain) == [&] {
    vi id(n);
    iota(all(id), 0);
    return id;
  }());
  chain[n - 1].push_back(0);
  CHECK(topo_sort(chain).empty());
  puts("OK graph/topo-sort");
}
