#include "../testlib.hpp"
#include "../../content/graph/dijkstra.hpp"

typedef vector<vector<pair<int, ll>>> Graph;

// Relax every edge n times (Bellman-Ford without early exit).
vector<ll> brute(const Graph& adj, const vi& srcs) {
  int n = sz(adj);
  vector<ll> d(n, LLONG_MAX);
  for (int s : srcs) d[s] = 0;
  rep(it, 0, n)
    rep(u, 0, n)
      if (d[u] < LLONG_MAX)
        for (auto [v, w] : adj[u]) d[v] = min(d[v], d[u] + w);
  return d;
}

// dist is feasible (no edge improves it), every reached vertex
// has a tight path back to a source, and path() follows it.
void certify(const Graph& adj, const vi& srcs, Dijkstra& dj) {
  int n = sz(adj);
  vector<bool> is_src(n);
  for (int s : srcs) {
    is_src[s] = true;
    CHECK_EQ(dj.dist[s], 0LL);
  }
  rep(u, 0, n)
    if (dj.dist[u] < LLONG_MAX)
      for (auto [v, w] : adj[u]) CHECK(dj.dist[v] <= dj.dist[u] + w);
  rep(t, 0, n) {
    vi p = dj.path(t);
    if (dj.dist[t] == LLONG_MAX) {
      CHECK(p.empty());
      CHECK_EQ(dj.par[t], -1);
      continue;
    }
    CHECK(!p.empty() && p.back() == t && is_src[p[0]]);
    CHECK_EQ(dj.par[p[0]], -1);
    rep(i, 1, sz(p)) {
      bool tight = false;
      for (auto [v, w] : adj[p[i - 1]])
        tight |= v == p[i] && dj.dist[p[i - 1]] + w == dj.dist[v];
      CHECK(tight);
      CHECK_EQ(dj.par[p[i]], p[i - 1]);
    }
  }
}

// The bottleneck Tweak (maximize the smallest edge) applied to a
// copy of the struct.
struct Widest {
  vector<ll> dist;
  vi par;
  Widest(const Graph& adj, const vi& srcs)
      : dist(sz(adj), -1), par(sz(adj), -1) {
    priority_queue<pair<ll, int>> pq;
    for (int s : srcs) pq.push({dist[s] = LLONG_MAX, s});
    while (sz(pq)) {
      auto [d, u] = pq.top();
      pq.pop();
      if (d < dist[u]) continue;
      for (auto [v, w] : adj[u])
        if (min(d, w) > dist[v]) pq.push({dist[v] = min(d, w), v}), par[v] = u;
    }
  }
  vi path(int t) {
    vi p;
    if (dist[t] >= 0)
      for (; t >= 0; t = par[t]) p.push_back(t);
    reverse(all(p));
    return p;
  }
};

void check_widest(const Graph& adj, const vi& srcs) {
  int n = sz(adj);
  vector<ll> want(n, -1);
  for (int s : srcs) want[s] = LLONG_MAX;
  rep(it, 0, n)
    rep(u, 0, n)
      for (auto [v, w] : adj[u]) want[v] = max(want[v], min(want[u], w));
  Widest wd(adj, srcs);
  CHECK(wd.dist == want);
  rep(t, 0, n) {
    vi p = wd.path(t);
    CHECK_EQ(p.empty(), want[t] < 0);
    if (p.empty()) continue;
    CHECK(count(all(srcs), p[0]) && p.back() == t);
    ll low = LLONG_MAX;
    rep(i, 1, sz(p)) {
      ll best = -1;
      for (auto [v, w] : adj[p[i - 1]])
        if (v == p[i]) best = max(best, w);
      low = min(low, best);
    }
    CHECK_EQ(low, want[t]);
  }
}

int main() {
  rep(it, 0, iters(3000)) {
    int n = (int)rnd(1, 8), m = (int)rnd(0, 20);
    ll maxw = rnd(0, 1) ? 3 : 1'000'000'000'000LL;
    Graph adj(n);
    rep(i, 0, m) {
      int u = (int)rnd(0, n - 1), v = (int)rnd(0, n - 1);
      adj[u].push_back({v, rnd(0, maxw)});
    }
    vi srcs;
    int k = (int)rnd(1, 3);
    rep(i, 0, k) srcs.push_back((int)rnd(0, n - 1));
    Dijkstra dj(adj, srcs);
    CHECK(dj.dist == brute(adj, srcs));
    certify(adj, srcs, dj);
    check_widest(adj, srcs);
  }
  // Bottleneck Tweak on 0 -> 1 (5), 1 -> 2 (3).
  Widest wd({{{1, 5}}, {{2, 3}}, {}}, {0});
  CHECK(wd.dist == vector<ll>({LLONG_MAX, 5, 3}));
  // Large: random graph, then a long path with huge weights.
  int n = 100'000;
  Graph adj(n);
  rep(i, 0, 200'000) {
    int u = (int)rnd(0, n - 1), v = (int)rnd(0, n - 1);
    adj[u].push_back({v, rnd(0, 1'000'000'000)});
  }
  Dijkstra dj(adj, {0});
  certify(adj, {0}, dj);
  Graph line(n);
  rep(i, 0, n - 1) line[i].push_back({i + 1, 1'000'000'000'000LL});
  Dijkstra dl(line, {0});
  CHECK_EQ(dl.dist[n - 1], (n - 1) * 1'000'000'000'000LL);
  CHECK_EQ(sz(dl.path(n - 1)), n);
  puts("OK graph/dijkstra");
}
