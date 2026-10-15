#include "../testlib.hpp"
#include "../../content/graph/euler-path.hpp"

bool brute(int u, int left, const vector<pii>& edges, bool directed,
           vector<bool>& used) {
  if (!left) return true;
  rep(e, 0, sz(edges))
    if (!used[e]) {
      auto [a, b] = edges[e];
      for (int dir = 0; dir < (directed ? 1 : 2); dir++, swap(a, b))
        if (a == u) {
          used[e] = true;
          bool ok = brute(b, left - 1, edges, directed, used);
          used[e] = false;
          if (ok) return true;
        }
    }
  return false;
}

void check_walk(const vector<pii>& edges, bool directed, int src, const vi& vs,
                const vi& es) {
  int m = sz(edges);
  CHECK_EQ(sz(vs), m + 1);
  CHECK_EQ(sz(es), m);
  CHECK_EQ(vs[0], src);
  vector<bool> seen(m);
  rep(i, 0, m) {
    int e = es[i];
    CHECK(0 <= e && e < m && !seen[e]);
    seen[e] = true;
    auto [a, b] = edges[e];
    bool fwd = a == vs[i] && b == vs[i + 1];
    bool bwd = b == vs[i] && a == vs[i + 1];
    CHECK(fwd || (!directed && bwd));
  }
}

// Edges of a random walk from src, shuffled: an Euler path exists.
vector<pii> walk_edges(int n, int m, int src, bool directed) {
  vector<pii> edges;
  int u = src;
  rep(i, 0, m) {
    int v = (int)rnd(0, n - 1);
    if (!directed && rnd(0, 1)) edges.push_back({v, u});
    else edges.push_back({u, v});
    u = v;
  }
  shuffle(all(edges), test_rng);
  return edges;
}

// Notes: orient an undirected graph so that |in - out| <= 1, with
// one circuit per component and dummy edges pairing odd vertices.
void check_orientation(int n, const vector<pii>& edges) {
  int m = sz(edges);
  vector<vi> adj(n);
  vi comp(n, -1), bal(n);
  for (auto [u, v] : edges) adj[u].push_back(v), adj[v].push_back(u);
  rep(s, 0, n) {
    if (comp[s] >= 0 || adj[s].empty()) continue;
    vi q = {s}, id, odd;
    comp[s] = s;
    rep(i, 0, sz(q))
      for (int y : adj[q[i]])
        if (comp[y] < 0) comp[y] = s, q.push_back(y);
    vector<pii> es;
    rep(e, 0, m)
      if (comp[edges[e].first] == s) es.push_back(edges[e]), id.push_back(e);
    for (int x : q)
      if (sz(adj[x]) % 2) odd.push_back(x);
    for (int i = 0; i < sz(odd); i += 2) es.push_back({odd[i], odd[i + 1]});
    auto [vs, ep] = euler_path(n, es, false, s);
    CHECK(sz(vs) == sz(es) + 1 && vs.back() == s);
    rep(i, 0, sz(ep))
      if (ep[i] < sz(id)) bal[vs[i]]++, bal[vs[i + 1]]--;
  }
  rep(v, 0, n) CHECK(abs(bal[v]) <= 1);
}

int main() {
  int found = 0, total = iters(4000);
  rep(it, 0, total) {
    int n = (int)rnd(1, 5), m = (int)rnd(0, 7), src = (int)rnd(0, n - 1);
    bool directed = rnd(0, 1);
    vector<pii> edges;
    if (rnd(0, 1)) edges = walk_edges(n, m, (int)rnd(0, n - 1), directed);
    else
      rep(i, 0, m) edges.push_back({(int)rnd(0, n - 1), (int)rnd(0, n - 1)});
    vector<bool> used(m);
    bool want = brute(src, m, edges, directed, used);
    auto [vs, es] = euler_path(n, edges, directed, src);
    CHECK_EQ(!vs.empty(), want);
    if (want) {
      found++;
      check_walk(edges, directed, src, vs, es);
    } else {
      CHECK(es.empty());
    }
  }
  CHECK(found > total / 5 && found < total * 4 / 5);
  rep(it, 0, iters(1000)) {
    int n = (int)rnd(1, 8), m = (int)rnd(0, 10);
    vector<pii> edges;
    rep(i, 0, m) edges.push_back({(int)rnd(0, n - 1), (int)rnd(0, n - 1)});
    check_orientation(n, edges);
  }
  // Two triangles: no single walk, but orientable per component.
  vector<pii> tri2 = {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}};
  CHECK(euler_path(6, tri2, false, 0).first.empty());
  check_orientation(6, tri2);
  // Large random walks (directed and undirected), and a long cycle.
  int n = 50'000, m = 100'000;
  for (bool directed : {false, true}) {
    vector<pii> edges = walk_edges(n, m, 7, directed);
    auto [vs, es] = euler_path(n, edges, directed, 7);
    check_walk(edges, directed, 7, vs, es);
    edges.push_back({n - 1, n - 2});
    edges.push_back({n - 3, n - 4});
    CHECK(euler_path(n, edges, directed, 7).first.empty());
  }
  vector<pii> cyc;
  rep(v, 0, m) cyc.push_back({v, (v + 1) % m});
  auto [vs, es] = euler_path(m, cyc, true, 0);
  check_walk(cyc, true, 0, vs, es);
  puts("OK graph/euler-path");
}
