#include "../testlib.hpp"
#include "../../content/flow/min-cost-flow.hpp"

struct Arc {
  int u, v;
  ll cap, cost;
};

// Exhaustive: every integer flow; returns {value, min cost} for
// value min(limit, max value).
pair<ll, ll> brute_exhaustive(int n, const vector<Arc>& arcs, int s, int t,
                              ll limit) {
  int m = sz(arcs);
  vector<ll> f(m, 0);
  map<ll, ll> best; // value -> min cost
  while (true) {
    vector<ll> net(n);
    ll c = 0;
    rep(i, 0, m)
      net[arcs[i].u] -= f[i], net[arcs[i].v] += f[i], c += f[i] * arcs[i].cost;
    bool ok = true;
    rep(v, 0, n)
      if (v != s && v != t && net[v]) ok = false;
    if (ok && net[t] >= 0) {
      auto it = best.find(net[t]);
      if (it == best.end() || c < it->second) best[net[t]] = c;
    }
    int i = 0;
    while (i < m && f[i] == arcs[i].cap) f[i++] = 0;
    if (i == m) break;
    f[i]++;
  }
  ll value = min(limit, best.rbegin()->first);
  return {value, best[value]};
}

// Successive shortest paths with plain Bellman-Ford on an edge list.
pair<ll, ll> brute_ssp(int n, const vector<Arc>& arcs, int s, int t, ll limit) {
  vi to, from;
  vector<ll> cap, cost;
  for (auto& a : arcs) {
    from.push_back(a.u), to.push_back(a.v), cap.push_back(a.cap),
        cost.push_back(a.cost);
    from.push_back(a.v), to.push_back(a.u), cap.push_back(0),
        cost.push_back(-a.cost);
  }
  ll flow = 0, total = 0;
  while (flow < limit) {
    vector<ll> d(n, LLONG_MAX);
    vi pe(n, -1);
    d[s] = 0;
    rep(round, 0, n)
      rep(e, 0, sz(to))
        if (cap[e] > 0 && d[from[e]] < LLONG_MAX &&
            d[from[e]] + cost[e] < d[to[e]]) {
          d[to[e]] = d[from[e]] + cost[e], pe[to[e]] = e;
        }
    if (d[t] == LLONG_MAX) break;
    ll push = limit - flow;
    for (int v = t; v != s; v = from[pe[v]]) push = min(push, cap[pe[v]]);
    for (int v = t; v != s; v = from[pe[v]])
      cap[pe[v]] -= push, cap[pe[v] ^ 1] += push;
    flow += push, total += push * d[t];
  }
  return {flow, total};
}

// Conservation, bounds, value and cost; optimality: no negative
// cycle in the residual graph, and no s-t path if below the limit.
void certify(MinCostFlow& mcf, int n, const vector<Arc>& arcs, const vi& ids,
             int s, int t, ll limit, ll f, ll c) {
  vector<ll> net(n);
  ll cost = 0;
  vector<tuple<int, int, ll>> res; // residual arcs
  rep(i, 0, sz(arcs)) {
    auto& a = arcs[i];
    ll x = mcf.flow(ids[i]);
    CHECK(0 <= x && x <= a.cap);
    net[a.u] -= x, net[a.v] += x, cost += x * a.cost;
    if (x < a.cap) res.push_back({a.u, a.v, a.cost});
    if (x > 0) res.push_back({a.v, a.u, -a.cost});
  }
  rep(v, 0, n)
    if (v != s && v != t) CHECK_EQ(net[v], 0LL);
  CHECK_EQ(net[t], f);
  CHECK_EQ(cost, c);
  vector<ll> d(n, 0);
  bool changed = true;
  for (int r = 0; r < n + 1 && changed; r++) {
    changed = false;
    for (auto [u, v, w] : res)
      if (d[u] + w < d[v]) d[v] = d[u] + w, changed = true;
  }
  CHECK(!changed);
  if (f < limit) {
    vector<bool> seen(n);
    seen[s] = true;
    for (bool grew = true; grew;) {
      grew = false;
      for (auto [u, v, w] : res)
        if (seen[u] && !seen[v]) seen[v] = grew = true;
    }
    CHECK(!seen[t]);
  }
}

// Costs phi[v] - phi[u] + slack: negative edges, no negative cycle.
vector<Arc> random_arcs(int n, int m, ll maxcap, ll maxslack) {
  vector<ll> phi(n);
  for (ll& x : phi) x = rnd(-5, 5);
  vector<Arc> arcs;
  rep(i, 0, m) {
    int u = (int)rnd(0, n - 1), v = (int)rnd(0, n - 1);
    arcs.push_back({u, v, rnd(0, maxcap), phi[v] - phi[u] + rnd(0, maxslack)});
  }
  return arcs;
}

pair<MinCostFlow, vi> build(int n, const vector<Arc>& arcs) {
  MinCostFlow mcf(n);
  vi ids;
  for (auto& a : arcs) ids.push_back(mcf.add_edge(a.u, a.v, a.cap, a.cost));
  return {mcf, ids};
}

int main() {
  rep(it, 0, iters(1500)) {
    int n = (int)rnd(2, 4), m = (int)rnd(0, 6);
    auto arcs = random_arcs(n, m, 2, rnd(0, 1) ? 0 : 4);
    int s = (int)rnd(0, n - 1), t = (int)rnd(0, n - 2);
    if (t >= s) t++;
    ll limit = rnd(0, 2) ? LLONG_MAX : rnd(0, 4);
    auto [mcf, ids] = build(n, arcs);
    auto [f, c] = mcf.max_flow(s, t, limit);
    CHECK(make_pair(f, c) == brute_exhaustive(n, arcs, s, t, limit));
    certify(mcf, n, arcs, ids, s, t, limit, f, c);
  }
  rep(it, 0, iters(1500)) {
    int n = (int)rnd(2, 12), m = (int)rnd(0, 30);
    auto arcs = random_arcs(n, m, rnd(1, 6), rnd(0, 6));
    ll limit = rnd(0, 2) ? LLONG_MAX : rnd(0, 10);
    auto [mcf, ids] = build(n, arcs);
    auto [f, c] = mcf.max_flow(0, n - 1, limit);
    CHECK(make_pair(f, c) == brute_ssp(n, arcs, 0, n - 1, limit));
    certify(mcf, n, arcs, ids, 0, n - 1, limit, f, c);
  }
  rep(it, 0, iters(1500)) { // Costs at the documented bound.
    int n = (int)rnd(2, 8), m = (int)rnd(0, 7);
    ll C = LLONG_MAX / (8 * n);
    vector<ll> phi(n);
    for (ll& x : phi) x = rnd(-C / 4, C / 4);
    vector<Arc> arcs;
    rep(i, 0, m) {
      int u = (int)rnd(0, n - 1), v = (int)rnd(0, n - 1);
      arcs.push_back({u, v, rnd(0, 1), phi[v] - phi[u] + rnd(0, C / 2)});
    }
    auto [mcf, ids] = build(n, arcs);
    auto [f, c] = mcf.max_flow(0, n - 1);
    CHECK(make_pair(f, c) == brute_ssp(n, arcs, 0, n - 1, LLONG_MAX));
    certify(mcf, n, arcs, ids, 0, n - 1, LLONG_MAX, f, c);
  }
  { // Extremes: total cost 4e18, negative edge on the cheap path.
    ll C = (ll)1e9, K = (ll)2e9;
    vector<Arc> arcs = {
        {0, 1, K, C}, {1, 3, K, C}, {0, 2, K, -C}, {2, 3, K, C}, {1, 2, K, -C}};
    auto [mcf, ids] = build(4, arcs);
    auto [f, c] = mcf.max_flow(0, 3);
    CHECK_EQ(f, 2 * K);
    CHECK_EQ(c, 2 * K * C);
    certify(mcf, 4, arcs, ids, 0, 3, LLONG_MAX, f, c);
  }
  elapsed();
  { // Planted assignment: rows -> columns with tight duals.
    int k = 150, n = 2 * k + 2, s = n - 2, t = n - 1;
    vi perm(k);
    iota(all(perm), 0);
    shuffle(all(perm), test_rng);
    vector<ll> u(k), v(k);
    ll opt = 0;
    rep(i, 0, k)
      u[i] = rnd(-1e6, 1e6), v[i] = rnd(-1e6, 1e6), opt += u[i] + v[i];
    vector<Arc> arcs;
    rep(i, 0, k) arcs.push_back({s, i, 1, 0}), arcs.push_back({k + i, t, 1, 0});
    rep(i, 0, k)
      rep(j, 0, k)
        arcs.push_back(
            {i, k + j, 1, u[i] + v[j] + (perm[i] == j ? 0 : rnd(1, 1e6))});
    auto [mcf, ids] = build(n, arcs);
    auto [f, c] = mcf.max_flow(s, t);
    CHECK_EQ(f, (ll)k);
    CHECK_EQ(c, opt);
    certify(mcf, n, arcs, ids, s, t, LLONG_MAX, f, c);
  }
  { // Random network, costs with negative edges.
    int n = 400;
    auto arcs = random_arcs(n, 4000, 20, 100);
    auto [mcf, ids] = build(n, arcs);
    auto [f, c] = mcf.max_flow(0, n - 1);
    certify(mcf, n, arcs, ids, 0, n - 1, LLONG_MAX, f, c);
  }
  { // Comb s -> a_i (cost i) -> x (cost -2i) -> tail -> t. Without
    // the Bellman-Ford potentials the first Dijkstra improves x k
    // times and rescans the tail each time: k * L = 2.5e9 steps.
    int k = 50000, L = 50000, x = k + 1, t = k + L + 2;
    MinCostFlow mcf(t + 1);
    rep(i, 1, k + 1) mcf.add_edge(0, i, 1, i), mcf.add_edge(i, x, 1, -2 * i);
    rep(j, 0, L) mcf.add_edge(x + j, x + j + 1, k, 0);
    mcf.add_edge(x + L, t, 1, 0);
    CHECK(mcf.max_flow(0, t) == make_pair(1LL, (ll)-k));
  }
  { // K cheap paths s -> x_i -> t next to a costly side chain with
    // 4e5 edges. Dijkstra must stop at t, not scan the chain K times.
    int K = 600, L = 20000, s = 0, t = 1, y = K + 2;
    MinCostFlow mcf(K + L + 2);
    rep(i, 0, K) mcf.add_edge(s, 2 + i, 1, i), mcf.add_edge(2 + i, t, 1, 0);
    mcf.add_edge(s, y, 1, 1e6);
    rep(j, 0, L - 1)
      rep(p, 0, 20) mcf.add_edge(y + j, y + j + 1, 1, 1 + p);
    double start = elapsed();
    CHECK(mcf.max_flow(s, t) == make_pair((ll)K, (ll)K * (K - 1) / 2));
    CHECK(elapsed() - start < 2); // without the early stop: ~25 s
  }
  fprintf(stderr, "  large cases: %.2fs\n", elapsed());
  puts("OK flow/min-cost-flow");
}
