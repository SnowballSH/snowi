#include "../testlib.hpp"
#include "../../content/graph/scc.hpp"

vector<vector<bool>> reach(const vector<vi>& adj) {
  int n = sz(adj);
  vector<vector<bool>> r(n, vector<bool>(n));
  rep(s, 0, n) {
    vi stack = {s};
    r[s][s] = true;
    while (sz(stack)) {
      int u = stack.back();
      stack.pop_back();
      for (int v : adj[u])
        if (!r[s][v]) r[s][v] = true, stack.push_back(v);
    }
  }
  return r;
}

void check_order(const vector<vi>& adj, const vi& comp, int k) {
  int n = sz(adj);
  CHECK_EQ(sz(comp), n);
  vector<bool> used(k);
  rep(v, 0, n) {
    CHECK(0 <= comp[v] && comp[v] < k);
    used[comp[v]] = true;
  }
  CHECK(count(all(used), false) == 0);
  rep(u, 0, n)
    for (int v : adj[u]) CHECK(comp[u] >= comp[v]);
}

int main() {
  rep(it, 0, iters(4000)) {
    int n = (int)rnd(0, 9), m = n ? (int)rnd(0, 2 * n + 2) : 0;
    vector<vi> adj(n);
    rep(i, 0, m) adj[rnd(0, n - 1)].push_back((int)rnd(0, n - 1));
    auto [comp, k] = scc(adj);
    auto r = reach(adj);
    check_order(adj, comp, k);
    rep(u, 0, n)
      rep(v, 0, n) CHECK_EQ(comp[u] == comp[v], r[u][v] && r[v][u]);
  }
  // Large: blocks of cycles chained into a path, labels shuffled.
  // The condensation is a path, so the numbering is forced.
  int n = 200'000, b = 1000, k = n / b;
  vi label(n);
  iota(all(label), 0);
  shuffle(all(label), test_rng);
  vector<vi> adj(n);
  rep(v, 0, n) {
    int next = v % b == b - 1 ? v - b + 1 : v + 1;
    adj[label[v]].push_back(label[next]);
    if (v + b < n) adj[label[v]].push_back(label[v + (int)rnd(1, b)]);
  }
  for (vi& out : adj) shuffle(all(out), test_rng);
  auto [comp, cnt] = scc(adj);
  CHECK_EQ(cnt, k);
  rep(v, 0, n) CHECK_EQ(comp[label[v]], k - 1 - v / b);
  // One cycle through all vertices: deep DFS, one component.
  vector<vi> cyc(n);
  rep(v, 0, n) cyc[v].push_back((v + 1) % n);
  CHECK_EQ(scc(cyc).second, 1);
  puts("OK graph/scc");
}
